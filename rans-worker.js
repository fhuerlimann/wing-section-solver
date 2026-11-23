/* Steady 2D incompressible RANS solver, Spalart–Allmaras turbulence model.
   Cell-centred finite volume on a body-fitted O-grid (optionally with a blunt trailing-edge base), SIMPLE pressure–velocity
   coupling with Rhie–Chow interpolation, line-by-line (ADI) TDMA solves.
   Non-dimensional: chord = 1, freestream speed = 1, density = 1, nu = 1/Re. */
"use strict";
if (typeof importScripts === "function") importScripts("sections.js");

const CB1 = 0.1355, CB2 = 0.622, SIG = 2 / 3, KAP = 0.41, CV1 = 7.1, CW2 = 0.3, CW3 = 2;
const CW1 = CB1 / (KAP * KAP) + (1 + CB2) / SIG;
const CV1_3 = CV1 ** 3, CW3_6 = CW3 ** 6;

/* ---------- geometry ---------- */
const GRID = { wakeSharp: true, wakeW: 10 };   // steer grid lines into the wake behind sharp trailing edges too
let gridRetries = 0;   // layers rebuilt because their grid lines converged (diagnostic)
/* Number of base points in each half of a blunt trailing edge. The base spacing is
   about twice the cosine spacing of the surface at the corners, so cell sizes match. */
function baseHalf(NI, te) {
  if (!(te > 0)) return 0;
  const dc = Math.PI * Math.PI / (NI * NI);                    // surface spacing at the TE, ~ (pi/ns)^2/4
  // thick wedge bases (> 2 % chord) may use up to twice as many base nodes
  return Math.max(2, Math.min(Math.floor(NI / (te > 0.02 ? 12 : 24)), Math.round(te / (2 * dc) / 2)));
}

/* Closed ring of NI nodes, counter-clockwise: base midpoint -> upper corner -> upper surface
   -> LE (index NI/2) -> lower surface -> lower corner -> base. With te = 0 the base vanishes and
   node 0 is the sharp trailing edge. te is the total trailing-edge thickness in chords, added
   as a linear thickness ramp so the section stays smooth. */
function nacaSurface(code, NI, invert, te) {
  te = te > 0 ? te : 0;
  const sh = sectionShape(code, invert, te), pt = sh.pt;   // NACA 4-digit or a built-in section (sections.js)
  const n = NI / 2, hb = baseHalf(NI, te), ns = n - hb;
  const xk = new Float64Array(ns + 1);
  for (let k = 0; k <= ns; k++) xk[k] = 0.5 * (1 - Math.cos(Math.PI * k / ns));
  // corners (the shoulder of a double wedge) must be mesh nodes, or they get cut off
  if (sh.snap) for (const xs of sh.snap) { let kb = 1; for (let k = 1; k < ns; k++) if (Math.abs(xk[k] - xs) < Math.abs(xk[kb] - xs)) kb = k; xk[kb] = xs; }
  const xc = k => xk[k];
  const X = new Float64Array(NI), Y = new Float64Array(NI);
  for (let k = ns; k >= 0; k--) { const q = pt(xc(k), true); X[hb + ns - k] = q[0]; Y[hb + ns - k] = q[1]; }
  for (let k = 1; k <= ns; k++) { const i = n + k; if (i >= NI) break; const q = pt(xc(k), false); X[i] = q[0]; Y[i] = q[1]; }
  if (hb > 0) {
    const ux = X[hb], uy = Y[hb], lx = X[NI - hb], ly = Y[NI - hb], mx = 0.5 * (ux + lx), my = 0.5 * (uy + ly);
    for (let i = 0; i < hb; i++) { X[i] = mx + (ux - mx) * i / hb; Y[i] = my + (uy - my) * i / hb; }
    for (let q = 1; q < hb; q++) { X[NI - hb + q] = lx + (mx - lx) * q / hb; Y[NI - hb + q] = ly + (my - ly) * q / hb; }
  }
  return { X, Y, hb };
}

function buildGrid(code, invert, NI, NJ, ds1, R, te) {
  const S = nacaSurface(code, NI, invert, te);
  // geometric wall-normal spacing from ds1 to radius R
  const M = NJ - 1;
  let lo = 1.0001, hi = 2;
  for (let k = 0; k < 80; k++) { const r = (lo + hi) / 2; (ds1 * (r ** M - 1) / (r - 1) > R) ? hi = r : lo = r; }
  const ratio = (lo + hi) / 2, d = new Float64Array(NJ);
  for (let j = 1; j < NJ; j++) d[j] = d[j - 1] + ds1 * ratio ** (j - 1);
  d[NJ - 1] = R;
  // March each layer along the normals of the previous one, then slide nodes
  // tangentially to even out spacing. Smoothing grows with wall distance, so the
  // near-wall cells stay orthogonal while the trailing-edge fan spreads out.
  const X = new Float64Array(NI * NJ), Y = new Float64Array(NI * NJ);
  // Wake band: base nodes and a few surface nodes either side of a blunt trailing edge.
  // Their grid lines are steered downstream (C-grid-like wake lines) and smoothed less,
  // so the wake keeps fine cross-stream resolution instead of a fan from one point.
  // sharp trailing edges get a wider band (no base nodes to carry the wake)
  const W = Math.max(6, Math.round(NI / (te > 0 ? 16 : GRID.wakeW))), wte = new Float64Array(NI);
  const round = (sharpSection(code) || {}).round === true;   // a circle has no trailing edge: plain radial lines
  const wake = !round && (te > 0 || GRID.wakeSharp);
  if (wake) for (let i = 0; i < NI; i++) { const e = Math.min(i, NI - i) - S.hb; wte[i] = e <= 0 ? 1 : e >= W ? 0 : 0.5 * (1 + Math.cos(Math.PI * e / W)); }
  const px = new Float64Array(NI), py = new Float64Array(NI), mx = new Float64Array(NI), my = new Float64Array(NI);
  const nx0 = new Float64Array(NI), ny0 = new Float64Array(NI);
  let carry = 0;   // extra smoothing carried to the next layer after a retry, decaying
  for (let i = 0; i < NI; i++) { X[i] = S.X[i]; Y[i] = S.Y[i]; }
  for (let j = 1; j < NJ; j++) {
    const o = (j - 1) * NI, h = d[j] - d[j - 1];
    for (let i = 0; i < NI; i++) {
      const ip = o + (i + 1) % NI, im = o + (i - 1 + NI) % NI;
      const tx = X[ip] - X[im], ty = Y[ip] - Y[im], l = Math.hypot(tx, ty);
      mx[i] = ty / l; my[i] = -tx / l;
    }
    nx0.set(mx); ny0.set(my);
    // Concave regions (the cusped aft lower surface of supercritical sections) focus the normals,
    // and the layers fold where they meet. A layer whose grid lines converge or fold is redone with
    // more normal smoothing and dent filling, and the next layers keep some of it. Meshes that
    // never trigger this (all but extreme-camber NACA 4-digit sections) are built exactly as before.
    for (let extra = carry; ; extra += 6) {
    if (extra > carry) { mx.set(nx0); my.set(ny0); }
    const sn = Math.min(1, d[j] / 0.05) * 0.5;               // normal smoothing
    for (let pass = 0; pass < 2 + extra; pass++) {
      for (let i = 0; i < NI; i++) {
        const ip = (i + 1) % NI, im = (i - 1 + NI) % NI;
        const ax = mx[i] + sn * 0.5 * (mx[ip] + mx[im] - 2 * mx[i]), ay = my[i] + sn * 0.5 * (my[ip] + my[im] - 2 * my[i]);
        const l = Math.hypot(ax, ay); px[i] = ax / l; py[i] = ay / l;
      }
      mx.set(px); my.set(py);
    }
    // far from the body, lean towards radial directions so concave regions cannot fold
    const wr = Math.min(1, Math.max(0, (d[j] - 0.05) / 1.5)), wrs = wr * wr * (3 - 2 * wr);
    if (wrs > 0) for (let i = 0; i < NI; i++) {
      const rx = X[o + i] - 0.5, ry = Y[o + i], rl = Math.hypot(rx, ry);
      const ax = (1 - wrs) * px[i] + wrs * rx / rl, ay = (1 - wrs) * py[i] + wrs * ry / rl, l = Math.hypot(ax, ay);
      px[i] = ax / l; py[i] = ay / l;
    }
    // steer the wake band towards the chord direction, blending to radial with the far field
    if (wake) for (let i = 0; i < NI; i++) {
      const b = wte[i] * Math.min(1, d[j] / 0.05);
      if (b > 0) {
        const rx = X[o + i] - 0.5, ry = Y[o + i], rl = Math.hypot(rx, ry), gx = (1 - wrs) + wrs * rx / rl, gy = wrs * ry / rl;
        const ax = (1 - b) * px[i] + b * gx, ay = (1 - b) * py[i] + b * gy, l = Math.hypot(ax, ay);
        px[i] = ax / l; py[i] = ay / l;
      }
    }
    const q = j * NI;
    for (let i = 0; i < NI; i++) { X[q + i] = X[o + i] + px[i] * h; Y[q + i] = Y[o + i] + py[i] * h; }
    const s = d[j] < 0.004 ? 0 : 0.35 * Math.min(1, (d[j] - 0.004) / 0.3);   // tangential smoothing
    for (let pass = 0; pass < 4 && s > 0; pass++) {
      for (let i = 0; i < NI; i++) { px[i] = X[q + i]; py[i] = Y[q + i]; }
      for (let i = 0; i < NI; i++) {
        const ip = (i + 1) % NI, im = (i - 1 + NI) % NI;
        const tx = px[ip] - px[im], ty = py[ip] - py[im], l = Math.hypot(tx, ty);
        const dx = 0.5 * (px[ip] + px[im]) - px[i], dy = 0.5 * (py[ip] + py[im]) - py[i];
        const t = (dx * tx + dy * ty) / l;
        const sd = s * (1 - 0.8 * wte[i]); X[q + i] = px[i] + sd * t * tx / l; Y[q + i] = py[i] + sd * t * ty / l;
      }
    }
    // On a retry, also fill dents in the new layer: a node lying inside the chord of its two
    // neighbours (behind a concavity) is moved outward onto it. Convex parts are not touched,
    // so the layer does not shrink.
    for (let pass = 0; pass < extra; pass++) {
      for (let i = 0; i < NI; i++) { px[i] = X[q + i]; py[i] = Y[q + i]; }
      for (let i = 0; i < NI; i++) {
        const ip = (i + 1) % NI, im = (i - 1 + NI) % NI;
        const tx = px[ip] - px[im], ty = py[ip] - py[im], l = Math.hypot(tx, ty), nx = ty / l, ny = -tx / l;
        const ln = (0.5 * (px[ip] + px[im]) - px[i]) * nx + (0.5 * (py[ip] + py[im]) - py[i]) * ny;
        if (ln > 0) { const w = 0.5 * (1 - wte[i]) * ln; X[q + i] = px[i] + w * nx; Y[q + i] = py[i] + w * ny; }
      }
    }
    let folded = false;
    for (let i = 0; i < NI && !folded; i++) {
      const a = o + i, b = o + (i + 1) % NI, c = q + (i + 1) % NI, e = q + i;
      if ((X[c] - X[a]) * (Y[e] - Y[b]) - (X[e] - X[b]) * (Y[c] - Y[a]) >= 0) folded = true;
    }
    // Converging grid lines (behind a concavity) are caught before they cross: no segment of the
    // new layer may shrink well below its length on the surface or on the previous layer.
    if (!folded && j > 1) for (let i = 0; i < NI && !folded; i++) {
      if (wte[i] > 0 || wte[(i + 1) % NI] > 0) continue;
      const a = o + i, b = o + (i + 1) % NI, c = q + (i + 1) % NI, e = q + i, ln = Math.hypot(X[c] - X[e], Y[c] - Y[e]);
      if (ln < 0.85 * Math.hypot(X[b] - X[a], Y[b] - Y[a]) || ln < 0.6 * Math.hypot(X[(i + 1) % NI] - X[i], Y[(i + 1) % NI] - Y[i])) folded = true;
    }
    if (folded) gridRetries++;
    if (!folded || extra >= 90) { carry = Math.floor(0.85 * extra); break; }
    }
  }
  return { X, Y, ratio, hb: S.hb, d, round };
}

/* ---------- tridiagonal solvers ---------- */
let TW = null;
function work(n) {
  if (!TW || TW.n < n) TW = { n, cp: new Float64Array(n), dp: new Float64Array(n), bb: new Float64Array(n), u: new Float64Array(n), z: new Float64Array(n) };
  return TW;
}
function thomas(n, a, b, c, d, x) {
  const { cp, dp } = work(n);
  cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
  for (let i = 1; i < n; i++) {
    const m = b[i] - a[i] * cp[i - 1];
    cp[i] = c[i] / m; dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
  }
  x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
}
// cyclic tridiagonal: alpha = A[n-1][0], beta = A[0][n-1]
function cyclic(n, a, b, c, alpha, beta, r, x) {
  const W = work(n), bb = W.bb, u = W.u, z = W.z;
  const gamma = -b[0];
  for (let i = 0; i < n; i++) bb[i] = b[i];
  bb[0] = b[0] - gamma; bb[n - 1] = b[n - 1] - alpha * beta / gamma;
  thomas(n, a, bb, c, r, x);
  u.fill(0); u[0] = gamma; u[n - 1] = alpha;
  thomas(n, a, bb, c, u, z);
  const fact = (x[0] + beta * x[n - 1] / gamma) / (1 + z[0] + beta * z[n - 1] / gamma);
  for (let i = 0; i < n; i++) x[i] -= fact * z[i];
}

/* ---------- agglomeration multigrid for the pressure-correction equation ----------
   Additive-correction multigrid (Hutchinson & Raithby): each coarse cell merges 2 x 2 fine
   cells, coarse coefficients are the sums of the fine links that cross the merged block
   (Galerkin with piecewise-constant transfer), restriction sums the residuals and
   prolongation injects the coarse correction. Alternating line TDMA is the smoother,
   which copes with the 1000:1 cells next to the wall. */
function lineSweep(ni, nc, A, b, phi, n, buf) {
  const aP = A.aP, aW = A.aW, aE = A.aE, aS = A.aS, aN = A.aN, la = buf[0], lb = buf[1], lc = buf[2], lr = buf[3], lx = buf[4];
  for (let s = 0; s < n; s++) {
    for (let i = 0; i < ni; i++) {                 // wall-normal lines
      const iw = (i - 1 + ni) % ni, ie = (i + 1) % ni;
      for (let j = 0; j < nc; j++) {
        const c = j * ni + i;
        la[j] = -aS[c]; lb[j] = aP[c]; lc[j] = -aN[c];
        lr[j] = b[c] + aW[c] * phi[j * ni + iw] + aE[c] * phi[j * ni + ie];
      }
      thomas(nc, la, lb, lc, lr, lx);
      for (let j = 0; j < nc; j++) phi[j * ni + i] = lx[j];
    }
    for (let j = 0; j < nc; j++) {                 // circumferential (periodic) lines
      const o = j * ni;
      for (let i = 0; i < ni; i++) {
        const c = o + i;
        la[i] = -aW[c]; lb[i] = aP[c]; lc[i] = -aE[c];
        let r = b[c];
        if (j > 0) r += aS[c] * phi[c - ni];
        if (j < nc - 1) r += aN[c] * phi[c + ni];
        lr[i] = r;
      }
      cyclic(ni, la, lb, lc, -aE[o + ni - 1], -aW[o], lr, lx);
      for (let i = 0; i < ni; i++) phi[o + i] = lx[i];
    }
  }
}

function makeMG(ni0, nc0, opt) {
  const levels = [];
  let ni = ni0, nc = nc0;
  while (ni % 2 === 0 && ni >= 24 && nc >= 8) {
    const L = { ni: ni / 2, nc: Math.ceil(nc / 2) };
    const n = L.ni * L.nc;
    for (const k of ["aP", "aW", "aE", "aS", "aN", "b", "x", "r"]) L[k] = new Float64Array(n);
    levels.push(L); ni = L.ni; nc = L.nc;
  }
  const mx = Math.max(ni0, nc0), buf = [0, 1, 2, 3, 4].map(() => new Float64Array(mx)), rFine = new Float64Array(ni0 * nc0);
  const nPre = opt.pre ?? 1, nPost = opt.post ?? 1, nCoarse = opt.coarse ?? 10;

  function coarsen(F, fni, fnc, C) {
    C.aP.fill(0); C.aW.fill(0); C.aE.fill(0); C.aS.fill(0); C.aN.fill(0);
    const cni = C.ni;
    for (let j = 0; j < fnc; j++) for (let i = 0; i < fni; i++) {
      const c = j * fni + i, I = i >> 1, J = j >> 1, K = J * cni + I;
      C.aP[K] += F.aP[c];
      // a = coupling of c to its east neighbour, b = of the east neighbour to c (equal for a
      // symmetric matrix; the compressible pressure correction is not symmetric)
      const ie = (i + 1) % fni, IE = ie >> 1, a = F.aE[c], b = F.aW[j * fni + ie];
      if (IE === I) C.aP[K] -= a + b; else { C.aE[K] += a; C.aW[J * cni + IE] += b; }
      if (j < fnc - 1) {
        const JN = (j + 1) >> 1, an = F.aN[c], bs = F.aS[c + fni];
        if (JN === J) C.aP[K] -= an + bs; else { C.aN[K] += an; C.aS[JN * cni + I] += bs; }
      }
    }
  }
  function residual(A, ni, nc, b, x, r) {
    for (let j = 0; j < nc; j++) for (let i = 0; i < ni; i++) {
      const c = j * ni + i;
      let s = b[c] - A.aP[c] * x[c] + A.aW[c] * x[j * ni + (i - 1 + ni) % ni] + A.aE[c] * x[j * ni + (i + 1) % ni];
      if (j > 0) s += A.aS[c] * x[c - ni];
      if (j < nc - 1) s += A.aN[c] * x[c + ni];
      r[c] = s;
    }
  }
  function vcycle(k, A, ni, nc, b, x) {
    if (k === levels.length) { lineSweep(ni, nc, A, b, x, nCoarse, buf); return; }
    lineSweep(ni, nc, A, b, x, nPre, buf);
    const C = levels[k], r = k === 0 ? rFine : levels[k - 1].r;
    residual(A, ni, nc, b, x, r);
    C.b.fill(0); C.x.fill(0);
    for (let j = 0; j < nc; j++) for (let i = 0; i < ni; i++) C.b[(j >> 1) * C.ni + (i >> 1)] += r[j * ni + i];
    vcycle(k + 1, C, C.ni, C.nc, C.b, C.x);
    for (let j = 0; j < nc; j++) for (let i = 0; i < ni; i++) x[j * ni + i] += C.x[(j >> 1) * C.ni + (i >> 1)];
    lineSweep(ni, nc, A, b, x, nPost, buf);
  }
  return {
    levels: levels.length,
    solve(A, b, x, cycles) {
      let F = A, fni = ni0, fnc = nc0;
      for (const C of levels) { coarsen(F, fni, fnc, C); F = C; fni = C.ni; fnc = C.nc; }
      for (let k = 0; k < cycles; k++) vcycle(0, A, ni0, nc0, b, x);
    },
  };
}

/* ---------- solution-adaptive mesh ----------
   The base O-grid is a map from computational coordinates (xi around the section, eta
   outwards, one unit per cell) to physical space. An adapted mesh keeps the NI x NJ
   topology and moves its nodes to new computational positions, so the map stays smooth
   and cannot fold. Grid lines are pulled towards cells where |vorticity| x (base cell size)
   is large: shear layers, the wake and shed vortices. The radial budget comes from the
   outer far field, the tangential budget from the quiet parts of each ring. The wall ring
   and the near-wall rings stay exactly where they were, so y+ and the boundary layer
   resolution do not change. */
const smooth01 = (a, b, x) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };
const ADAPT = { eps: 0.04, wmax: 4, extra: 0.25, kmax: 1.6, bt0: 0.02, bt1: 0.4, win: 0.3, rsm: 4, kappa: 300, epsS: 0.01, wmaxS: 6, bt0S: 0.003, bt1S: 0.03, kappaS: 30 };

// Catmull-Rom weights
function crw(t, w) {
  const t2 = t * t, t3 = t2 * t;
  w[0] = 0.5 * (-t3 + 2 * t2 - t); w[1] = 0.5 * (3 * t3 - 5 * t2 + 2); w[2] = 0.5 * (-3 * t3 + 4 * t2 + t); w[3] = 0.5 * (t3 - t2);
}
// physical position of computational point (xi, eta) on the base node lattice (bicubic, periodic in xi)
function evalBase(B, xi, eta, out) {
  const NI = B.NI, NJ = B.NJ, X = B.X, Y = B.Y, wx = B._wx || (B._wx = new Float64Array(4)), wy = B._wy || (B._wy = new Float64Array(4));
  xi = ((xi % NI) + NI) % NI;
  const i0 = Math.floor(xi);
  let j0 = Math.floor(eta); if (j0 > NJ - 2) j0 = NJ - 2; if (j0 < 0) j0 = 0;
  crw(xi - i0, wx); crw(eta - j0, wy);
  let x = 0, y = 0;
  for (let b = 0; b < 4; b++) {
    const j = j0 - 1 + b, jj = j < 0 ? 0 : j > NJ - 1 ? NJ - 1 : j, jn = j < 0 ? 1 : NJ - 2;   // linear extrapolation outside
    let rx = 0, ry = 0;
    for (let a = 0; a < 4; a++) {
      const i = (i0 - 1 + a + NI) % NI;
      let px = X[jj * NI + i], py = Y[jj * NI + i];
      if (j !== jj) { px = 2 * px - X[jn * NI + i]; py = 2 * py - Y[jn * NI + i]; }
      rx += wx[a] * px; ry += wx[a] * py;
    }
    x += wy[b] * rx; y += wy[b] * ry;
  }
  out[0] = x; out[1] = y;
}
// computational coordinates of the cell centres of a mesh given by its node coordinates (XI unwrapped: node i ~ i)
function cellComp(NI, NJ, XI, ET) {
  const NC = NJ - 1, A = new Float64Array(NI * NC), Bv = new Float64Array(NI * NC);
  for (let j = 0; j < NC; j++) for (let i = 0; i < NI; i++) {
    const ip = (i + 1) % NI, sh = ip === 0 ? NI : 0, a = j * NI + i, b = j * NI + ip, e = (j + 1) * NI + ip, f = (j + 1) * NI + i;
    A[a] = 0.25 * (XI[a] + XI[b] + sh + XI[e] + sh + XI[f]);
    Bv[a] = 0.25 * (ET[a] + ET[b] + ET[e] + ET[f]);
  }
  return { A, B: Bv };
}
/* Bilinear interpolation on a structured, curvilinear lattice of points (A, B), periodic in i:
   finds the lattice quad containing (xi, eta) by walking from the previous one. */
function makeLocator(NI, NR, A, Bv) {
  let gi = 0, gj = 0;
  const r = { c00: 0, c10: 0, c01: 0, c11: 0, s: 0, t: 0 };
  function at(xi, eta) {
    let i = gi, j = gj, s = 0.5, t = 0.5;
    for (let it = 0; it < 3 * (NI + NR); it++) {
      const ip = (i + 1) % NI, sh = ip === 0 ? NI : 0, c00 = j * NI + i, c10 = j * NI + ip, c01 = c00 + NI, c11 = c10 + NI;
      const x00 = A[c00], x10 = A[c10] + sh, x01 = A[c01], x11 = A[c11] + sh, y00 = Bv[c00], y10 = Bv[c10], y01 = Bv[c01], y11 = Bv[c11];
      let x = xi; while (x - x00 > NI / 2) x -= NI; while (x - x00 < -NI / 2) x += NI;
      s = 0.5; t = 0.5;
      for (let k = 0; k < 8; k++) {
        const px = (1 - s) * (1 - t) * x00 + s * (1 - t) * x10 + (1 - s) * t * x01 + s * t * x11 - x;
        const py = (1 - s) * (1 - t) * y00 + s * (1 - t) * y10 + (1 - s) * t * y01 + s * t * y11 - eta;
        const xs = (1 - t) * (x10 - x00) + t * (x11 - x01), xt = (1 - s) * (x01 - x00) + s * (x11 - x10);
        const ys = (1 - t) * (y10 - y00) + t * (y11 - y01), yt = (1 - s) * (y01 - y00) + s * (y11 - y10);
        const det = xs * yt - xt * ys;
        if (!(Math.abs(det) > 1e-14)) break;
        const ds = (px * yt - py * xt) / det, dt = (py * xs - px * ys) / det;
        s -= ds; t -= dt;
        if (Math.abs(ds) + Math.abs(dt) < 1e-10) break;
      }
      let moved = false;
      if (s < -1e-7) { i = (i - 1 + NI) % NI; moved = true; } else if (s > 1 + 1e-7) { i = ip; moved = true; }
      if (t < -1e-7 && j > 0) { j--; moved = true; } else if (t > 1 + 1e-7 && j < NR - 2) { j++; moved = true; }
      if (!moved) break;
    }
    gi = i; gj = j;
    s = s < 0 ? 0 : s > 1 ? 1 : s; t = t < 0 ? 0 : t > 1 ? 1 : t;
    const ip = (i + 1) % NI;
    r.c00 = j * NI + i; r.c10 = j * NI + ip; r.c01 = r.c00 + NI; r.c11 = r.c10 + NI; r.s = s; r.t = t;
    return r;
  }
  at.val = (f, q) => (1 - q.s) * (1 - q.t) * f[q.c00] + q.s * (1 - q.t) * f[q.c10] + (1 - q.s) * q.t * f[q.c01] + q.s * q.t * f[q.c11];
  return at;
}
const identityComp = (NI, NJ) => {
  const XI = new Float64Array(NI * NJ), ET = new Float64Array(NI * NJ);
  for (let j = 0; j < NJ; j++) for (let i = 0; i < NI; i++) { XI[j * NI + i] = i; ET[j * NI + i] = j; }
  return { XI, ET, NJ };
};

/* New node positions. B: base mesh {NI, NJ, X, Y, d}; cur: computational node coordinates of the
   mesh the solution {u, v} lives on. Returns {X, Y, comp, bad, stats}.
   The indicator is the change of total pressure per base cell in each grid direction, measured
   on the current mesh. Total pressure only changes where the flow has lost energy (boundary layer,
   shear layers, wake, vortex cores), not where the outer flow merely accelerates; around the rings
   it marks wakes and shed vortices that cross them, but not the boundary layer, which changes
   only slowly along the wall. */
function adaptGrid(B, cur, sol, strength) {
  const { NI, NJ, d, X, Y } = B, NC = NJ - 1, N = NI * NC, o = ADAPT;
  const KC = cur.NJ - 1, NK = NI * KC, p0 = new Float64Array(NK);   // current mesh: KC rings of cells
  let shocks = false;   // any supersonic flow: then shocks are refined (sensor, whole ring, thinner wall band)
  const isComp = !!sol.comp, GM = 1.4, M2 = isComp ? sol.mach * sol.mach : 0, pa = isComp ? new Float64Array(NK) : null;
  for (let c = 0; c < NK; c++) {
    const q2 = sol.u[c] * sol.u[c] + sol.v[c] * sol.v[c];
    if (!isComp) { p0[c] = 2 * sol.p[c] + sol.u[c] * sol.u[c] + sol.v[c] * sol.v[c]; continue; }   // total pressure / q
    // compressible: isentropic total pressure / q (it drops across shocks as well as in viscous layers)
    const T = Math.max(0.05, 1 + (GM - 1) * M2 * (sol.hh[c] - 0.5 * (q2 - 1))), P = Math.max(0.02, 1 + GM * M2 * sol.p[c]);
    pa[c] = P; p0[c] = 2 * P * Math.pow(1 + 0.5 * (GM - 1) * M2 * q2 / T, GM / (GM - 1)) / (GM * M2);
    if (M2 * q2 / T > 1) shocks = true;
  }
  // 1. total-pressure change per unit of base computational coordinate, on the current cells
  const cc = cellComp(NI, cur.NJ, cur.XI, cur.ET), A = cc.A, Bc = cc.B, gT = new Float64Array(NK), gR = new Float64Array(NK);
  const sTc = new Float64Array(NK), sRc = new Float64Array(NK);
  for (let j = 0; j < KC; j++) for (let i = 0; i < NI; i++) {
    const c = j * NI + i, ip = j * NI + (i + 1) % NI, im = j * NI + (i - 1 + NI) % NI;
    let da = A[ip] - A[im]; da = ((da % NI) + NI) % NI;
    gT[c] = Math.abs(p0[ip] - p0[im]) / Math.max(da, 1e-6);
    const jp = j < KC - 1 ? c + NI : c, jm = j > 0 ? c - NI : c;
    gR[c] = Math.abs(p0[jp] - p0[jm]) / Math.max(Bc[jp] - Bc[jm], 1e-6);
    if (shocks) {
      // Shock sensor: normalised second difference of pressure along each grid direction. A captured
      // shock is a jump spread over a few cells, so this stays large however fine the cells get, while
      // smooth compressions and expansions fade as the mesh is refined. Scaled so that eps -> epsS.
      const sT = Math.abs(pa[ip] - 2 * pa[c] + pa[im]) / (pa[ip] + 2 * pa[c] + pa[im]);
      const sR = j > 0 && j < KC - 1 ? Math.abs(pa[jp] - 2 * pa[c] + pa[jm]) / (pa[jp] + 2 * pa[c] + pa[jm]) : 0;
      sTc[c] = sT / o.epsS; sRc[c] = sR / o.epsS;
    }
  }
  // 2. weights on the base cells: wanted refinement radially (wR) and around the rings (wT)
  const loc = makeLocator(NI, KC, A, Bc), wR = new Float64Array(N), wT = new Float64Array(N);
  const wRS = shocks ? new Float64Array(N) : null, wTS = shocks ? new Float64Array(N) : null;   // shock weights
  for (let j = 0; j < NC; j++) {
    const dc = 0.5 * (d[j] + d[j + 1]), mR = smooth01(0.01, 0.04, dc) * (1 - smooth01(1, 2, dc));   // not in the boundary layer, nor where the far field gives up cells
    for (let i = 0; i < NI; i++) {
      const c = j * NI + i, q = loc(i + 0.5, j + 0.5);
      wR[c] = 1 + mR * (Math.min(o.wmax, Math.max(1, loc.val(gR, q) / o.eps)) - 1);
      wT[c] = Math.min(o.wmax, Math.max(1, loc.val(gT, q) / o.eps));
      if (shocks) {
        wRS[c] = 1 + mR * (Math.min(o.wmaxS, Math.max(1, loc.val(sRc, q))) - 1);
        wTS[c] = Math.min(o.wmaxS, Math.max(1, loc.val(sTc, q)));
      }
    }
  }
  // widen the refined zones a little, then smooth: along the rings for wR, along the grid lines for wT
  const tmp = new Float64Array(N);
  const dilate = w => {
    tmp.set(w);
    for (let j = 0; j < NC; j++) for (let i = 0; i < NI; i++) {
      const c = j * NI + i; let m = tmp[c];
      m = Math.max(m, tmp[j * NI + (i + 1) % NI], tmp[j * NI + (i - 1 + NI) % NI]);
      if (j > 0) m = Math.max(m, tmp[c - NI]); if (j < NC - 1) m = Math.max(m, tmp[c + NI]);
      w[c] = m;
    }
  };
  const smoothI = (w, n) => { for (let k = 0; k < n; k++) { tmp.set(w); for (let j = 0; j < NC; j++) for (let i = 0; i < NI; i++) w[j * NI + i] = 0.25 * (tmp[j * NI + (i + 1) % NI] + tmp[j * NI + (i - 1 + NI) % NI]) + 0.5 * tmp[j * NI + i]; } };
  const smoothJ = (w, n) => { for (let k = 0; k < n; k++) { tmp.set(w); for (let j = 0; j < NC; j++) for (let i = 0; i < NI; i++) { const c = j * NI + i, jm = j > 0 ? c - NI : c, jp = j < NC - 1 ? c + NI : c; w[c] = 0.25 * (tmp[jm] + tmp[jp]) + 0.5 * tmp[c]; } } };
  dilate(wR); dilate(wR); dilate(wT); dilate(wT);
  // The wake band behind the trailing edge packs many grid lines into a thin strip. Rings must move
  // together across it, or its cells shear (and URANS decouples there): one radial weight per ring.
  const nb = B.hb + Math.max(6, Math.round(NI / 16));
  if (!B.round && (B.hb > 0 || GRID.wakeSharp)) for (let j = 0; j < NC; j++) {
    let m = 0;
    for (let q = -nb; q <= nb; q++) m = Math.max(m, wR[j * NI + (q + NI) % NI]);
    for (let q = -nb; q <= nb; q++) wR[j * NI + (q + NI) % NI] = m;
  }
  smoothI(wR, 12); smoothJ(wR, 4); smoothJ(wT, 12); smoothI(wT, 8);
  // Shocks are thin: their weights are widened once and smoothed mostly along the shock, so the
  // refinement stays concentrated across it. The larger of the two weights wins.
  if (shocks) {
    dilate(wRS); dilate(wTS);
    smoothI(wRS, 6); smoothJ(wRS, 2); smoothJ(wTS, 6); smoothI(wTS, 2);
    for (let c = 0; c < N; c++) { wR[c] = Math.max(wR[c], wRS[c]); wT[c] = Math.max(wT[c], wTS[c]); }
  }
  for (let c = 0; c < N; c++) { wR[c] = 1 + strength * (wR[c] - 1); wT[c] = 1 + strength * (wT[c] - 1); }

  // 3. radial redistribution along every grid line i. The adapted mesh gets up to extra * NC more
  // rings; the far field is never coarsened (coarsening the far wake made URANS unstable). Columns
  // that need fewer extra rings than the others put the rest beyond ~2 chords.
  const g = new Float64Array(NC);
  for (let j = 0; j < NC; j++) g[j] = smooth01(1.5, 4, 0.5 * (d[j] + d[j + 1]));
  let Emax = 0;
  for (let i = 0; i < NI; i++) {
    let E = 0;
    for (let j = 0; j < NC; j++) E += 0.5 * (wR[j * NI + i] + wR[j * NI + (i - 1 + NI) % NI]) - 1;
    if (E > Emax) Emax = E;
  }
  const add = Math.min(Math.round(o.extra * NC), Math.ceil(Emax)), NCn = NC + add, NJn = NCn + 1;
  const s = Emax > 0 ? Math.min(1, add / Emax) : 0;
  const col = new Float64Array(NC), ETn = new Float64Array(NI * NJn);
  for (let i = 0; i < NI; i++) {
    let Sw = 0, Swg = 0;
    for (let j = 0; j < NC; j++) { col[j] = 1 + s * (0.5 * (wR[j * NI + i] + wR[j * NI + (i - 1 + NI) % NI]) - 1); Sw += col[j]; Swg += col[j] * g[j]; }
    const lam = 1 + (NCn - Sw) / Swg;   // >= 1
    for (let j = 0; j < NC; j++) col[j] *= 1 + (lam - 1) * g[j];
    for (let pass = 0; pass < o.rsm; pass++) {   // even out the ring-to-ring size changes (keeps the sum)
      let prev = col[0];
      for (let j = 0; j < NC; j++) { const c0 = col[j], nx = j < NC - 1 ? col[j + 1] : c0; col[j] = 0.25 * (prev + 2 * c0 + nx); prev = c0; }
    }
    // invert the cumulative measure: node k sits where the measure reaches k
    let j = 0, acc = 0;
    ETn[i] = 0;
    for (let k = 1; k < NCn; k++) {
      while (j < NC - 1 && acc + col[j] < k) { acc += col[j]; j++; }
      ETn[k * NI + i] = j + Math.min(1, Math.max(0, (k - acc) / col[j]));
    }
    ETn[NCn * NI + i] = NC;
  }
  // Neighbouring grid lines can be very close (wake band, trailing edge), so independent radial
  // distributions would shear the cells between them. Smooth each ring across the grid lines with a
  // stiffness (hR / hT)^2: implicit, (I + kappa L) eta = eta0, which keeps the rings in order.
  {
    const la = new Float64Array(NI), lb = new Float64Array(NI), lc = new Float64Array(NI), lr = new Float64Array(NI), lx = new Float64Array(NI), cs = new Float64Array(NI);
    for (let k = 1; k < NCn; k++) {
      let em = 0; for (let i = 0; i < NI; i++) em += ETn[k * NI + i];
      const j0 = Math.min(NC - 1, Math.max(0, Math.floor(em / NI)));
      for (let i = 0; i < NI; i++) {   // stiffness of the link i -> i+1
        const ip = (i + 1) % NI, a = j0 * NI + i, b = j0 * NI + ip, e = (j0 + 1) * NI + i, f = (j0 + 1) * NI + ip;
        const hT = Math.hypot(0.5 * (X[b] + X[f] - X[a] - X[e]), 0.5 * (Y[b] + Y[f] - Y[a] - Y[e]));
        const hR = Math.hypot(0.5 * (X[e] + X[f] - X[a] - X[b]), 0.5 * (Y[e] + Y[f] - Y[a] - Y[b]));
        cs[i] = o.kappa * Math.min(1e4, (hR / hT) ** 2);
      }
      for (let i = 0; i < NI; i++) {
        const cm = cs[(i - 1 + NI) % NI], cp = cs[i];
        la[i] = -cm; lc[i] = -cp; lb[i] = 1 + cm + cp; lr[i] = ETn[k * NI + i];
      }
      cyclic(NI, la, lb, lc, -cs[NI - 1], -cs[NI - 1], lr, lx);
      for (let i = 0; i < NI; i++) ETn[k * NI + i] = lx[i];
    }
  }
  // 4. tangential redistribution along every new ring k, sampling wT along its path
  // With supersonic flow, shocks can stand anywhere around the section (on the upper surface, ahead
  // of the nose), so the whole ring is redistributed, anchored at the leading-edge node.
  const wn = shocks ? NI / 2 : Math.round(o.win * NI / 2), XIn = new Float64Array(NI * NJn), seg = new Float64Array(NI), xe = new Float64Array(NI + 1);
  const wTat = (x, y) => {   // bilinear on the base cell-centre lattice
    let xx = x - 0.5, yy = y - 0.5; if (yy < 0) yy = 0; if (yy > NC - 1) yy = NC - 1;
    xx = ((xx % NI) + NI) % NI;
    const i0 = Math.floor(xx), j0 = Math.min(NC - 2, Math.floor(yy)), fx = xx - i0, fy = yy - j0, i1 = (i0 + 1) % NI;
    return (1 - fx) * (1 - fy) * wT[j0 * NI + i0] + fx * (1 - fy) * wT[j0 * NI + i1] + (1 - fx) * fy * wT[(j0 + 1) * NI + i0] + fx * fy * wT[(j0 + 1) * NI + i1];
  };
  const dAt = e => { const j = Math.min(NC - 1, Math.max(0, Math.floor(e))), f = e - j; return d[j] + f * (d[j + 1] - d[j]); };
  for (let k = 0; k < NJn; k++) {
    let em = 0;
    for (let i = 0; i < NI; i++) em += ETn[k * NI + i];
    // keep the near-wall lines orthogonal; a shock reaches down to the boundary layer, so with
    // supersonic flow the lines may move from just outside it
    const bl = k === 0 ? 0 : shocks ? smooth01(o.bt0S, o.bt1S, dAt(em / NI)) : smooth01(o.bt0, o.bt1, dAt(em / NI));
    if (bl === 0) { for (let i = 0; i < NI; i++) XIn[k * NI + i] = i; continue; }
    // only the sector behind the trailing edge (nodes -wn..wn) is redistributed; its ends stay put,
    // so grid lines are never dragged round from the far side of the section
    for (let i = 0; i < NI; i++) XIn[k * NI + i] = i;
    let Sx = 0;
    for (let m = 0; m < 2 * wn; m++) {
      const i = m - wn, eA = ETn[k * NI + (i + NI) % NI], eB = ETn[k * NI + (i + 1 + NI) % NI];
      seg[m] = wTat(i + 0.5, 0.5 * (eA + eB)); Sx += seg[m] - 1;
    }
    const sk = Sx > (o.kmax - 1) * 2 * wn ? (o.kmax - 1) * 2 * wn / Sx : 1;
    let tot = 0;
    for (let m = 0; m < 2 * wn; m++) { seg[m] = 1 + sk * (seg[m] - 1); xe[m] = tot; tot += seg[m]; }
    xe[2 * wn] = tot;
    let q = 0;
    for (let m = 1; m < 2 * wn; m++) {
      const target = m * tot / (2 * wn);
      while (q < 2 * wn - 1 && xe[q + 1] < target) q++;
      const xq = q + (target - xe[q]) / seg[q];
      XIn[k * NI + (m - wn + NI) % NI] = (m - wn) + bl * (xq - m) + (m - wn < 0 ? NI : 0);
    }
  }
  // the same along each grid line: smooth the tangential positions from ring to ring with a stiffness
  // (hT / hR)^2, so that closely spaced rings cannot shear against each other. One stiffness per ring
  // (averaged over the sector), so all grid lines get the same smoothing and cannot cross.
  if (wn > 0) {
    const M = NJn - 1, la = new Float64Array(M), lb = new Float64Array(M), lc = new Float64Array(M), lr = new Float64Array(M), lx = new Float64Array(M), cs = new Float64Array(NJn);
    for (let k = 0; k < NJn - 1; k++) {   // link k -> k+1
      let sum = 0;
      for (let q = -wn; q < wn; q++) {
        const i = (q + NI) % NI, ip = (i + 1) % NI, j0 = Math.min(NC - 1, Math.max(0, Math.floor(ETn[k * NI + i])));
        const a = j0 * NI + i, b = j0 * NI + ip, e = (j0 + 1) * NI + i, f = (j0 + 1) * NI + ip;
        const hT = Math.hypot(0.5 * (X[b] + X[f] - X[a] - X[e]), 0.5 * (Y[b] + Y[f] - Y[a] - Y[e]));
        const hR = Math.hypot(0.5 * (X[e] + X[f] - X[a] - X[b]), 0.5 * (Y[e] + Y[f] - Y[a] - Y[b]));
        sum += Math.log(Math.max(1e-4, Math.min(1e4, (hT / hR) ** 2)));
      }
      cs[k] = (shocks ? o.kappaS : o.kappa) * Math.exp(sum / (2 * wn));   // geometric mean
    }
    for (let q = -wn; q <= wn; q++) {
      const i = (q + NI) % NI;
      for (let k = 1; k < NJn; k++) {
        const r = k - 1, cm = cs[k - 1], cp = k < NJn - 1 ? cs[k] : 0;
        la[r] = -cm; lc[r] = -cp; lb[r] = 1 + cm + cp; lr[r] = XIn[k * NI + i];
        if (k === 1) { lr[r] += cm * XIn[i]; la[r] = 0; }
      }
      thomas(M, la, lb, lc, lr, lx);
      for (let k = 1; k < NJn; k++) XIn[k * NI + i] = lx[k - 1];
    }
  }
  // 5. nodes: (xi, eta) with eta taken from the radial distribution at that xi; then physical positions
  const comp = { XI: new Float64Array(NI * NJn), ET: new Float64Array(NI * NJn), NJ: NJn }, Xn = new Float64Array(NI * NJn), Yn = new Float64Array(NI * NJn), pt = [0, 0];
  for (let k = 0; k < NJn; k++) for (let i = 0; i < NI; i++) {
    const n = k * NI + i, xi = XIn[n], xw = ((xi % NI) + NI) % NI, i0 = Math.floor(xw), f = xw - i0, i1 = (i0 + 1) % NI;
    const et = (1 - f) * ETn[k * NI + i0] + f * ETn[k * NI + i1];
    comp.XI[n] = xi; comp.ET[n] = et;
    evalBase(B, xi, et, pt); Xn[n] = pt[0]; Yn[n] = pt[1];
  }
  // 6. validity: every cell must keep the orientation of the base mesh
  let sgn = 0, bad = 0;
  const area = (X, Y, j, i) => { const ip = (i + 1) % NI, a = j * NI + i, b = j * NI + ip, e = (j + 1) * NI + ip, f = (j + 1) * NI + i; return (X[e] - X[a]) * (Y[f] - Y[b]) - (Y[e] - Y[a]) * (X[f] - X[b]); };
  for (let i = 0; i < NI; i++) sgn += area(X, Y, 0, i);
  sgn = Math.sign(sgn);
  for (let j = 0; j < NCn; j++) for (let i = 0; i < NI; i++) if (!(area(Xn, Yn, j, i) * sgn > 0)) bad++;
  let refR = 0, refT = 0;
  for (let c = 0; c < N; c++) { refR = Math.max(refR, 1 + s * (wR[c] - 1)); refT = Math.max(refT, wT[c]); }
  return { X: Xn, Y: Yn, NJ: NJn, comp, bad, stats: { rings: add, scale: s, refR, refT, shocks } };
}

/* ---------- the case ---------- */
function createCase(P) {
  const NI = P.NI, NJ = P.NJ, NC = NJ - 1, N = NI * NC;
  const nu = 1 / P.Re, R = P.R || 40;
  const cfEst = 0.026 / P.Re ** (1 / 7);
  const ds1 = P.ds1 || (P.yplus || 1) / (P.Re * Math.sqrt(0.5 * cfEst));
  // P.grid: an adapted mesh with the same topology (see adaptGrid)
  const G = P.grid || buildGrid(P.code, P.invert, NI, NJ, ds1, R, P.te || 0), X = G.X, Y = G.Y, HB = G.hb;
  const F64 = () => new Float64Array(N);

  // cells
  const xc = F64(), yc = F64(), vol = F64();
  let neg = 0, pos = 0;
  for (let j = 0; j < NC; j++) for (let i = 0; i < NI; i++) {
    const c = j * NI + i, ip = (i + 1) % NI;
    const a = j * NI + i, b = j * NI + ip, e = (j + 1) * NI + ip, f = (j + 1) * NI + i;
    xc[c] = (X[a] + X[b] + X[e] + X[f]) / 4; yc[c] = (Y[a] + Y[b] + Y[e] + Y[f]) / 4;
    const A = 0.5 * ((X[e] - X[a]) * (Y[f] - Y[b]) - (Y[e] - Y[a]) * (X[f] - X[b]));
    vol[c] = Math.abs(A); A < 0 ? neg++ : pos++;
  }
  const badCells = Math.min(neg, pos);

  // neighbours
  const nW = new Int32Array(N), nE = new Int32Array(N), nS = new Int32Array(N), nN = new Int32Array(N);
  for (let j = 0; j < NC; j++) for (let i = 0; i < NI; i++) {
    const c = j * NI + i;
    nW[c] = j * NI + (i - 1 + NI) % NI; nE[c] = j * NI + (i + 1) % NI;
    nS[c] = j > 0 ? c - NI : -1; nN[c] = j < NC - 1 ? c + NI : -1;
  }

  // i-faces (f = j*NI + i): between L = cell (i-1, j) and R = cell (i, j)
  const iSx = F64(), iSy = F64(), iFx = F64(), iFy = F64(), iW = F64(), iG = F64();
  for (let j = 0; j < NC; j++) for (let i = 0; i < NI; i++) {
    const f = j * NI + i, a = j * NI + i, b = (j + 1) * NI + i, L = nW[f], Rc = f;
    let sx = Y[b] - Y[a], sy = -(X[b] - X[a]);
    const dx = xc[Rc] - xc[L], dy = yc[Rc] - yc[L];
    if (sx * dx + sy * dy < 0) { sx = -sx; sy = -sy; }
    const fx = (X[a] + X[b]) / 2, fy = (Y[a] + Y[b]) / 2;
    const dL = Math.hypot(fx - xc[L], fy - yc[L]), dR = Math.hypot(fx - xc[Rc], fy - yc[Rc]);
    iSx[f] = sx; iSy[f] = sy; iFx[f] = fx; iFy[f] = fy; iW[f] = dR / (dL + dR);
    iG[f] = (sx * sx + sy * sy) / (sx * dx + sy * dy);
  }
  // j-faces (f = j*NI + i, j = 0..NC): j=0 wall, j=NC far field; normal points outward (+j)
  const NF = NI * NJ;
  const jSx = new Float64Array(NF), jSy = new Float64Array(NF), jFx = new Float64Array(NF), jFy = new Float64Array(NF), jW = new Float64Array(NF), jG = new Float64Array(NF);
  for (let j = 0; j < NJ; j++) for (let i = 0; i < NI; i++) {
    const f = j * NI + i, ip = (i + 1) % NI, a = j * NI + i, b = j * NI + ip;
    let sx = Y[b] - Y[a], sy = -(X[b] - X[a]);
    const fx = (X[a] + X[b]) / 2, fy = (Y[a] + Y[b]) / 2;
    let dx, dy, w = 0.5;
    if (j === 0) { dx = xc[i] - fx; dy = yc[i] - fy; }
    else if (j === NC) { const L = (NC - 1) * NI + i; dx = fx - xc[L]; dy = fy - yc[L]; }
    else {
      const L = (j - 1) * NI + i, Rc = j * NI + i;
      dx = xc[Rc] - xc[L]; dy = yc[Rc] - yc[L];
      const dL = Math.hypot(fx - xc[L], fy - yc[L]), dR = Math.hypot(fx - xc[Rc], fy - yc[Rc]);
      w = dR / (dL + dR);
    }
    if (sx * dx + sy * dy < 0) { sx = -sx; sy = -sy; }
    jSx[f] = sx; jSy[f] = sy; jFx[f] = fx; jFy[f] = fy; jW[f] = w;
    jG[f] = (sx * sx + sy * sy) / (sx * dx + sy * dy);
  }

  // wall distance
  const wd = F64();
  for (let c = 0; c < N; c++) {
    let best = Infinity;
    for (let i = 0; i < NI; i++) {
      const ip = (i + 1) % NI, ax = X[i], ay = Y[i], bx = X[ip] - ax, by = Y[ip] - ay;
      let t = ((xc[c] - ax) * bx + (yc[c] - ay) * by) / (bx * bx + by * by);
      t = t < 0 ? 0 : t > 1 ? 1 : t;
      const ex = xc[c] - ax - t * bx, ey = yc[c] - ay - t * by, d2 = ex * ex + ey * ey;
      if (d2 < best) best = d2;
    }
    wd[c] = Math.sqrt(best);
  }

  // state
  const u = F64(), v = F64(), p = F64(), nt = F64(), mt = F64(), pc = F64();
  const Fi = F64(), Fj = new Float64Array(NF);
  const ubF = new Float64Array(NI), vbF = new Float64Array(NI), inF = new Uint8Array(NI);
  const gux = F64(), guy = F64(), gvx = F64(), gvy = F64(), gpx = F64(), gpy = F64(), gtx = F64(), gty = F64();
  const aP = F64(), aW = F64(), aE = F64(), aS = F64(), aN = F64(), bu = F64(), bv = F64(), aPu = F64(), dP = F64(), mres = F64(), sNb = F64(), dR = F64();
  const ciF = F64(), cjF = new Float64Array(NF);   // pressure-correction face coefficients
  const ntInf = 3 * nu;
  /* Compressible (all-speed) mode, freestream Mach number P.mach > 0. Still non-dimensional with
     chord, U∞, ρ∞ and T∞; p is the pressure relative to p∞ in units of ρ∞U∞², so the absolute
     pressure is p∞(1 + γM²p) and ρ = (1 + γM²p) / T. Mass fluxes Mi, Mj = ρ_face × volume flux.
     hh: total enthalpy minus its freestream value, in U∞²; T = 1 + (γ-1)M²(hh - (|u|² - 1)/2).
     With M = 0 the arrays alias the incompressible ones and the code paths are unchanged. */
  const MACH = P.mach > 0 ? P.mach : 0, comp = MACH > 0, GAM = 1.4, M2 = MACH * MACH, SUTH = 110.4 / 288.15;
  const inviscid = !!P.inviscid;
  // The Mach number is ramped up over the first iterations of a run from uniform flow: the early
  // SIMPLE iterations overshoot the velocity several times, which at full Mach number would
  // drive the temperature to zero. M2r, Mr: the Mach number in use.
  let M2r = M2, Mr = MACH;
  const RAMP = P.machRamp || 300;
  function rampMach() {
    // Starts at 1/10 of the target Mach number. At the impulsive start the velocity overshoots (up to
    // ~30x freestream at a blunt base corner), so a supersonic target first holds Mach 0.05 for 30
    // iterations until that has decayed, then ramps as usual.
    const sup = M2 > 1, f = Math.min(1, Math.max(0, S.it + 1 - (sup ? 30 : 0)) / RAMP);
    const m2 = M2 * Math.max(sup ? 0.0025 / M2 : 0.01, f * f);
    if (m2 === M2r) return;
    M2r = m2; Mr = Math.sqrt(m2); farBC();
  }
  const rho = comp ? F64().fill(1) : null, Tc = comp ? F64().fill(1) : null, hh = comp ? F64() : null;
  const muL = comp ? F64().fill(nu) : null, me = comp ? F64() : null, psi = comp ? F64() : null;
  const grx = comp ? F64() : null, gry = comp ? F64() : null;
  const Mi = comp ? F64() : Fi, Mj = comp ? new Float64Array(NF) : Fj;
  const rfI = comp ? F64() : null, rfJ = comp ? new Float64Array(NF) : null;
  // far field: boundary density, pressure (relative to p∞) and whether p is imposed there
  const rhoB = comp ? new Float64Array(NI) : null, pB = comp ? new Float64Array(NI) : null, pDir = new Uint8Array(NI);
  const sutherland = T => T * Math.sqrt(T) * (1 + SUTH) / (T + SUTH);
  // k-omega SST (Menter, Kuntz & Langtry 2003) with ambient terms (SST-sust, Spalart & Rumsey 2007)
  const sst = (P.turb || "sa") === "sst" || P.turb === "sstg";
  // Menter, Smirnov, Liu & Avancha (2015): one-equation local correlation-based transition model (γ) on top of SST
  const trans = P.turb === "sstg";
  const KLp = P.kl ?? trans;   // Kato–Launder production (S·Ω): no spurious k at stagnation points
  const SST = { b1: 0.075, b2: 0.0828, bs: 0.09, sk1: 0.85, sk2: 1.0, sw1: 0.5, sw2: 0.856, g1: 5 / 9, g2: 0.44, a1: 0.31 };
  // with the transition model, ω∞ follows the usual setup for it: νt/ν = 1 in the free stream (NASA verification cases)
  const kInf = 1.5 * (P.tu ?? 1e-3) ** 2, wInf = P.wInf || (P.turb === "sstg" ? kInf * (P.Re || 1e6) : 5);   // Tu = 0.1 %, omega = 5 U/c: nu_t / nu below 2 for Re up to 6 million
  const tk = sst ? F64() : null, tw = sst ? F64() : null, F1a = sst ? F64() : null, S2a = sst ? F64() : null, gam = sst || comp ? F64() : null;
  const gkx = sst ? F64() : null, gky = sst ? F64() : null, gwx = sst ? F64() : null, gwy = sst ? F64() : null;
  const gmi = trans ? F64() : null, nwx = trans ? F64() : null, nwy = trans ? F64() : null, vnA = trans ? F64() : null, ReVa = trans ? F64() : null;
  let nwDone = false, warm = false, gReleased = false, holdFrac = 0, holdOn = false;   // warm: continuing from a loaded solution, no start-up freeze
  const DBGT = P.debugT ? { TuL: F64(), Rec: F64(), Fon: F64(), RT: F64(), ReV: F64(), lam: F64() } : null;
  const wWall = new Float64Array(NI);
  for (let i = 0; i < NI; i++) wWall[i] = 6 * nu / (SST.b1 * wd[i] * wd[i]);   // Wilcox near-wall solution, fixed in the wall cells
  const S = { cdw: undefined, alpha: 0, gust: 0, ca: 1, sa: 0, it: 0, cl: 0, cd: 0, cdp: 0, cdf: 0, cm: 0, res0: null, resC: 1, resU: 1, ypMax: 0, ypMean: 0, pinf: 0 };
  // SIMPLEC (default) drops the neighbour-correction term more consistently than SIMPLE,
  // which allows a much larger pressure under-relaxation factor.
  const simplec = (P.scheme || "simplec") === "simplec";
  // Multigrid converges the pressure correction well enough for velocity relaxation 0.95
  // with two V-cycles per iteration; the plain line sweeps need 0.85.
  const useMG = P.mg !== false;
  // SST is stiffer than SA: 0.95 often diverges with it, 0.9 does not
  const relax = { u: P.relaxU || (P.mach > 0.6 ? 0.7 : 0) || (simplec ? (useMG ? (sst ? 0.9 : 0.95) : 0.85) : 0.7), p: P.relaxP || (simplec ? 1.0 : 0.3), t: 0.8 };
  const nI = NI, nNC = NC;
  const mg = useMG ? makeMG(NI, NC, P.mgOpt || {}) : null, mgA = { aP, aW, aE, aS, aN };
  const lineA = new Float64Array(Math.max(NI, NC)), lineB = new Float64Array(Math.max(NI, NC)), lineC = new Float64Array(Math.max(NI, NC)), lineR = new Float64Array(Math.max(NI, NC)), lineX = new Float64Array(Math.max(NI, NC));

  function setAlpha(deg) {
    S.alpha = deg; const a = deg * Math.PI / 180; S.ca = Math.cos(a); S.sa = Math.sin(a);
    farBC();
  }
  function farBC() {
    const G0 = P.noVortex ? 0 : 0.5 * S.cl; // clockwise circulation from Kutta–Joukowski, Cl = 2 Gamma / (U c)
    const gc = Math.cos(S.gust), gs = Math.sin(S.gust);   // URANS start-up gust: inflow turned by S.gust (forces keep the nominal axes)
    if (comp) { farBCcomp(G0, gc, gs); return; }
    for (let i = 0; i < NI; i++) {
      const f = NC * NI + i, rx = jFx[f] - 0.25, ry = jFy[f], r2 = rx * rx + ry * ry;
      const k = G0 / (2 * Math.PI * r2);
      ubF[i] = S.ca * gc - S.sa * gs + k * ry; vbF[i] = S.sa * gc + S.ca * gs - k * rx;
      inF[i] = (ubF[i] * jSx[f] + vbF[i] * jSy[f]) < 0 ? 1 : 0;
    }
  }
  /* Compressible far field, following the characteristics normal to the boundary:
     subsonic inflow imposes velocity and temperature (total enthalpy) and takes the pressure
     from inside; subsonic outflow imposes the pressure and takes the velocity from inside;
     supersonic inflow imposes everything, supersonic outflow nothing. The imposed state is
     the freestream plus the point vortex of the current lift in its Prandtl–Glauert form
     (subsonic only; in supersonic flow the far field ahead of the section is undisturbed),
     with p and ρ isentropic from the local speed. */
  function farBCcomp(G0, gc, gs) {
    const sub = Mr < 1, be = sub ? Math.sqrt(1 - M2r) : 1, Gk = sub ? G0 : 0;
    const ex = S.ca * gc - S.sa * gs, ey = S.sa * gc + S.ca * gs;   // freestream direction
    for (let i = 0; i < NI; i++) {
      const f = NC * NI + i, rx = jFx[f] - 0.25, ry = jFy[f];
      // vortex in freestream-aligned coordinates: u' = Γβ/2π · y/(x² + β²y²), v' = -Γβ/2π · x/(x² + β²y²)
      const xs = rx * ex + ry * ey, ys = -rx * ey + ry * ex, den = xs * xs + be * be * ys * ys;
      const k = Gk * be / (2 * Math.PI * den), us = 1 + k * ys, vs = -k * xs;
      ubF[i] = us * ex - vs * ey; vbF[i] = us * ey + vs * ex;
      const q2 = ubF[i] * ubF[i] + vbF[i] * vbF[i], T = 1 - 0.5 * (GAM - 1) * M2r * (q2 - 1);
      rhoB[i] = Math.pow(T, 1 / (GAM - 1)); pB[i] = (Math.pow(T, GAM / (GAM - 1)) - 1) / (GAM * M2r);
      const un = (ubF[i] * jSx[f] + vbF[i] * jSy[f]) / Math.hypot(jSx[f], jSy[f]), Mn = Math.abs(un) * Mr / Math.sqrt(T);
      inF[i] = un < 0 ? 1 : 0;
      pDir[i] = (un < 0) === (Mn >= 1) ? 1 : 0;   // supersonic inflow or subsonic outflow
    }
  }
  function initFlow() {
    S.cl = 0; S.cd = 0; farBC();   // also clears a far field left over from a diverged run
    for (let c = 0; c < N; c++) { u[c] = S.ca; v[c] = S.sa; p[c] = 0; nt[c] = ntInf; mt[c] = ntInf * fv1(ntInf / nu); }
    if (sst) for (let c = 0; c < N; c++) { tk[c] = kInf; tw[c] = wInf; mt[c] = kInf / wInf; }
    // Start laminar near the body: with γ = 1 there, the impulsive start trips one side of the
    // leading edge within a few iterations and the turbulent boundary layer then sustains itself.
    if (trans) for (let c = 0; c < N; c++) gmi[c] = wd[c] < (P.gInitD ?? 0.05) ? (P.gInit ?? 0.02) : 1;
    gReleased = false;
    for (let f = 0; f < N; f++) Fi[f] = u[f] * iSx[f] + v[f] * iSy[f];
    for (let f = 0; f < NF; f++) Fj[f] = f < NI ? 0 : S.ca * jSx[f] + S.sa * jSy[f];
    if (comp) { S.it = 0; M2r = -1; rampMach(); hh.fill(0); thermo(true); massFromVolume(); }
    S.it = 0; S.res0 = null;
  }
  // temperature, density and laminar viscosity from p, u, v and hh (relax < 1 under-relaxes ρ)
  function thermo(exact) {
    const wr = exact ? 1 : (P.relaxRho || 0.7);
    for (let c = 0; c < N; c++) {
      const q2 = u[c] * u[c] + v[c] * v[c];
      let T = 1 + (GAM - 1) * M2r * (hh[c] - 0.5 * (q2 - 1));
      if (T < 0.05) T = 0.05;
      let pa = 1 + GAM * M2r * p[c];
      if (pa < 0.02) pa = 0.02;
      Tc[c] = T; rho[c] += wr * (pa / T - rho[c]);
      muL[c] = inviscid ? 0 : nu * sutherland(T);
      me[c] = muL[c] + rho[c] * (inviscid ? 0 : mt[c]);
      psi[c] = GAM * M2r / T;
    }
  }
  // mass fluxes from the volume fluxes with central face densities (start-up and loaded states only)
  function massFromVolume() {
    for (let f = 0; f < N; f++) { const L = nW[f], w = iW[f]; rfI[f] = w * rho[L] + (1 - w) * rho[f]; Mi[f] = rfI[f] * Fi[f]; }
    for (let f = 0; f < NF; f++) {
      if (f < NI) { Mj[f] = 0; rfJ[f] = rho[f]; continue; }
      const L = f - NI;
      if (f >= NC * NI) { const i = f - NC * NI; rfJ[f] = inF[i] ? rhoB[i] : rho[L]; Mj[f] = rfJ[f] * Fj[f]; continue; }
      const w = jW[f]; rfJ[f] = w * rho[L] + (1 - w) * rho[f]; Mj[f] = rfJ[f] * Fj[f];
    }
  }
  function fv1(chi) { const c3 = chi * chi * chi; return c3 / (c3 + CV1_3); }

  /* Green–Gauss gradient. wallMode: 0 -> value 0 at wall, 1 -> zero gradient. far: array of inflow values or null */
  function grad(phi, wallMode, far, gx, gy, mask) {
    const msk = mask || inF;
    gx.fill(0); gy.fill(0);
    for (let f = 0; f < N; f++) {
      const L = nW[f], w = iW[f], pf = w * phi[L] + (1 - w) * phi[f];
      const sx = pf * iSx[f], sy = pf * iSy[f];
      gx[L] += sx; gy[L] += sy; gx[f] -= sx; gy[f] -= sy;
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, w = jW[f], pf = w * phi[L] + (1 - w) * phi[f];
      const sx = pf * jSx[f], sy = pf * jSy[f];
      gx[L] += sx; gy[L] += sy; gx[f] -= sx; gy[f] -= sy;
    }
    for (let i = 0; i < NI; i++) {
      const pw = wallMode === 0 ? 0 : phi[i];
      gx[i] -= pw * jSx[i]; gy[i] -= pw * jSy[i];
      const f = NC * NI + i, L = (NC - 1) * NI + i;
      const pb = far && msk[i] ? far[i] : phi[L];
      gx[L] += pb * jSx[f]; gy[L] += pb * jSy[f];
    }
    for (let c = 0; c < N; c++) { gx[c] /= vol[c]; gy[c] /= vol[c]; }
  }

  function sweep(b, phi, nsweep) {
    for (let s = 0; s < nsweep; s++) {
      for (let i = 0; i < NI; i++) {               // wall-normal lines
        for (let j = 0; j < NC; j++) {
          const c = j * NI + i;
          lineA[j] = -aS[c]; lineB[j] = aP[c]; lineC[j] = -aN[c];
          lineR[j] = b[c] + aW[c] * phi[nW[c]] + aE[c] * phi[nE[c]];
        }
        thomas(NC, lineA, lineB, lineC, lineR, lineX);
        for (let j = 0; j < NC; j++) phi[j * NI + i] = lineX[j];
      }
      for (let j = 0; j < NC; j++) {               // circumferential (periodic) lines
        const o = j * NI;
        for (let i = 0; i < NI; i++) {
          const c = o + i;
          lineA[i] = -aW[c]; lineB[i] = aP[c]; lineC[i] = -aE[c];
          let r = b[c];
          if (j > 0) r += aS[c] * phi[c - NI];
          if (j < NC - 1) r += aN[c] * phi[c + NI];
          lineR[i] = r;
        }
        cyclic(NI, lineA, lineB, lineC, -aE[o + NI - 1], -aW[o], lineR, lineX);
        for (let i = 0; i < NI; i++) phi[o + i] = lineX[i];
      }
    }
  }

  const clampHO = (h, a, b) => a < b ? (h < a ? a : h > b ? b : h) : (h < b ? b : h > a ? a : h);

  function momentum() {
    aP.fill(0); aW.fill(0); aE.fill(0); aS.fill(0); aN.fill(0); bu.fill(0); bv.fill(0);
    const face = (L, Rc, F, D, fx, fy, isI) => {
      const aLR = D + (F > 0 ? F : 0), aRL = D + (F < 0 ? -F : 0);
      if (isI) { aE[L] = aRL; aW[Rc] = aLR; } else { aN[L] = aRL; aS[Rc] = aLR; }
      aP[L] += aLR; aP[Rc] += aRL;
      if (P.firstOrder) return;
      const U = F > 0 ? L : Rc, dx = fx - xc[U], dy = fy - yc[U];
      const uh = clampHO(u[U] + gux[U] * dx + guy[U] * dy, u[L], u[Rc]);
      const vh = clampHO(v[U] + gvx[U] * dx + gvy[U] * dy, v[L], v[Rc]);
      const cu = F * (uh - u[U]), cv = F * (vh - v[U]);
      bu[L] -= cu; bu[Rc] += cu; bv[L] -= cv; bv[Rc] += cv;
    };
    if (comp) {   // mass fluxes, effective dynamic viscosity
      for (let f = 0; f < N; f++) { const L = nW[f], w = iW[f]; face(L, f, Mi[f], (w * me[L] + (1 - w) * me[f]) * iG[f], iFx[f], iFy[f], true); }
      for (let f = NI; f < NC * NI; f++) { const L = f - NI, w = jW[f]; face(L, f, Mj[f], (w * me[L] + (1 - w) * me[f]) * jG[f], jFx[f], jFy[f], false); }
    } else if (inviscid) {
      for (let f = 0; f < N; f++) face(nW[f], f, Fi[f], 0, iFx[f], iFy[f], true);
      for (let f = NI; f < NC * NI; f++) face(f - NI, f, Fj[f], 0, jFx[f], jFy[f], false);
    } else {
    for (let f = 0; f < N; f++) {
      const L = nW[f], w = iW[f];
      face(L, f, Fi[f], (nu + w * mt[L] + (1 - w) * mt[f]) * iG[f], iFx[f], iFy[f], true);
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, w = jW[f];
      face(L, f, Fj[f], (nu + w * mt[L] + (1 - w) * mt[f]) * jG[f], jFx[f], jFy[f], false);
    }
    }
    for (let i = 0; i < NI; i++) {
      if (!inviscid) aP[i] += (comp ? muL[i] : nu) * jG[i];   // no-slip wall (slip wall when inviscid)
      const f = NC * NI + i, L = (NC - 1) * NI + i, F = Mj[f];
      if (inF[i]) {
        const D = (inviscid ? 0 : comp ? me[L] : nu + mt[L]) * jG[f], a = D + (F < 0 ? -F : 0);
        aP[L] += D + (F > 0 ? F : 0); bu[L] += a * ubF[i]; bv[L] += a * vbF[i];
      } else {
        if (F > 0) aP[L] += F; else { bu[L] -= F * u[L]; bv[L] -= F * v[L]; }
      }
    }
    let ru = 0;
    for (let c = 0; c < N; c++) {
      bu[c] -= gpx[c] * vol[c]; bv[c] -= gpy[c] * vol[c];
      if (tr && comp) {   // d(ρu)/dt, BDF2
        const k = vol[c] / tr.dt, r1 = tr.r1[c], r2 = tr.r2[c];
        if (tr.n > 0) { aP[c] += 1.5 * k * rho[c]; bu[c] += k * (2 * r1 * tr.u1[c] - 0.5 * r2 * tr.u2[c]); bv[c] += k * (2 * r1 * tr.v1[c] - 0.5 * r2 * tr.v2[c]); }
        else { aP[c] += k * rho[c]; bu[c] += k * r1 * tr.u1[c]; bv[c] += k * r1 * tr.v1[c]; }
      } else if (tr) {   // physical time derivative: BDF2 (backward Euler on the first step)
        const k = vol[c] / tr.dt;
        if (tr.n > 0) { aP[c] += 1.5 * k; bu[c] += k * (2 * tr.u1[c] - 0.5 * tr.u2[c]); bv[c] += k * (2 * tr.v1[c] - 0.5 * tr.v2[c]); }
        else { aP[c] += k; bu[c] += k * tr.u1[c]; bv[c] += k * tr.v1[c]; }
      }
      let r = bu[c] + aW[c] * u[nW[c]] + aE[c] * u[nE[c]] - aP[c] * u[c];
      if (nS[c] >= 0) r += aS[c] * u[nS[c]];
      if (nN[c] >= 0) r += aN[c] * u[nN[c]];
      ru += Math.abs(r);
      sNb[c] = aW[c] + aE[c] + aS[c] + aN[c];
      const apr = aP[c] / relax.u;
      bu[c] += (apr - aP[c]) * u[c]; bv[c] += (apr - aP[c]) * v[c];
      aP[c] = apr; aPu[c] = apr;
    }
    sweep(bu, u, 1);
    sweep(bv, v, 1);
    return ru;
  }

  function pressure() {
    // dR: Rhie–Chow smoothing (same for both schemes, so the converged answer does not depend on the scheme)
    // dP: pressure-correction coefficient (SIMPLEC uses aP minus the neighbour sum)
    for (let c = 0; c < N; c++) {
      dR[c] = 0.7 * vol[c] / (relax.u * aPu[c]);   // = 0.7 V / aP(unrelaxed): independent of the relaxation factor
      if (simplec) { const den = aPu[c] - sNb[c], lo = 0.1 * aPu[c]; dP[c] = vol[c] / (den > lo ? den : lo); }
      else dP[c] = dR[c];
    }
    mres.fill(0);
    for (let f = 0; f < N; f++) {
      const L = nW[f], w = iW[f], w1 = 1 - w;
      const cr = (w * dR[L] + w1 * dR[f]) * iG[f], co = (w * dP[L] + w1 * dP[f]) * iG[f];
      const dx = xc[f] - xc[L], dy = yc[f] - yc[L];
      const F = (w * u[L] + w1 * u[f]) * iSx[f] + (w * v[L] + w1 * v[f]) * iSy[f]
        - cr * ((p[f] - p[L]) - ((w * gpx[L] + w1 * gpx[f]) * dx + (w * gpy[L] + w1 * gpy[f]) * dy));
      Fi[f] = F; ciF[f] = co; mres[L] += F; mres[f] -= F;
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, w = jW[f], w1 = 1 - w;
      const cr = (w * dR[L] + w1 * dR[f]) * jG[f], co = (w * dP[L] + w1 * dP[f]) * jG[f];
      const dx = xc[f] - xc[L], dy = yc[f] - yc[L];
      const F = (w * u[L] + w1 * u[f]) * jSx[f] + (w * v[L] + w1 * v[f]) * jSy[f]
        - cr * ((p[f] - p[L]) - ((w * gpx[L] + w1 * gpx[f]) * dx + (w * gpy[L] + w1 * gpy[f]) * dy));
      Fj[f] = F; cjF[f] = co; mres[L] += F; mres[f] -= F;
    }
    let fin = 0, fout = 0;
    for (let i = 0; i < NI; i++) {
      const f = NC * NI + i, L = (NC - 1) * NI + i;
      const F = inF[i] ? ubF[i] * jSx[f] + vbF[i] * jSy[f] : u[L] * jSx[f] + v[L] * jSy[f];
      Fj[f] = F;
      if (inF[i]) fin += F; else fout += F;
    }
    const sc = fout > 1e-12 ? -fin / fout : 1;
    for (let i = 0; i < NI; i++) {
      const f = NC * NI + i, L = (NC - 1) * NI + i;
      if (!inF[i]) Fj[f] *= sc;
      mres[L] += Fj[f];
    }
    // pressure-correction equation
    aP.fill(0); aW.fill(0); aE.fill(0); aS.fill(0); aN.fill(0);
    for (let f = 0; f < N; f++) { const L = nW[f], co = ciF[f]; aE[L] = co; aW[f] = co; aP[L] += co; aP[f] += co; }
    for (let f = NI; f < NC * NI; f++) { const L = f - NI, co = cjF[f]; aN[L] = co; aS[f] = co; aP[L] += co; aP[f] += co; }
    let rc = 0;
    for (let c = 0; c < N; c++) { rc += Math.abs(mres[c]); bu[c] = -mres[c]; pc[c] = 0; }
    if (useMG) mg.solve(mgA, bu, pc, P.mgCycles || 2);
    else sweep(bu, pc, P.pSweeps || (NI >= 384 ? 8 : NI >= 256 ? 6 : 4));   // finer meshes need more sweeps to stay stable with SIMPLEC
    const ref = pc[(NC - 1) * NI];
    for (let c = 0; c < N; c++) pc[c] -= ref;
    for (let f = 0; f < N; f++) Fi[f] -= ciF[f] * (pc[f] - pc[nW[f]]);
    for (let f = NI; f < NC * NI; f++) Fj[f] -= cjF[f] * (pc[f] - pc[f - NI]);
    grad(pc, 1, null, gtx, gty);
    for (let c = 0; c < N; c++) {
      u[c] -= dP[c] * gtx[c]; v[c] -= dP[c] * gty[c]; p[c] += relax.p * pc[c];
    }
    return rc;
  }

  /* Compressible pressure correction (all-speed SIMPLEC). The mass flux through a face is
     ρ_f F with the Rhie–Chow volume flux F and a bounded upwind-biased face density ρ_f. A
     pressure correction p' changes both: m' = -ρ_f d Δp' + F ψ_up p'_up with ψ = ∂ρ/∂p = γM²/T.
     The second term makes the equation convective where the flow is fast: at low Mach number
     it vanishes and the equation is the incompressible one, in supersonic regions it carries
     the correction downstream, which is what lets shocks form. */
  function pressureComp() {
    for (let c = 0; c < N; c++) {
      dR[c] = 0.7 * vol[c] / (relax.u * aPu[c]);
      if (simplec) { const den = aPu[c] - sNb[c], lo = 0.1 * aPu[c]; dP[c] = vol[c] / (den > lo ? den : lo); }
      else dP[c] = dR[c];
    }
    mres.fill(0);
    const faceRho = (L, Rc, F, fx, fy) => {
      const U = F > 0 ? L : Rc;
      return clampHO(rho[U] + grx[U] * (fx - xc[U]) + gry[U] * (fy - yc[U]), rho[L], rho[Rc]);
    };
    for (let f = 0; f < N; f++) {
      const L = nW[f], w = iW[f], w1 = 1 - w;
      const cr = (w * dR[L] + w1 * dR[f]) * iG[f], co = (w * dP[L] + w1 * dP[f]) * iG[f];
      const dx = xc[f] - xc[L], dy = yc[f] - yc[L];
      const F = (w * u[L] + w1 * u[f]) * iSx[f] + (w * v[L] + w1 * v[f]) * iSy[f]
        - cr * ((p[f] - p[L]) - ((w * gpx[L] + w1 * gpx[f]) * dx + (w * gpy[L] + w1 * gpy[f]) * dy));
      const r = faceRho(L, f, F, iFx[f], iFy[f]);
      Fi[f] = F; ciF[f] = co; rfI[f] = r; Mi[f] = r * F; mres[L] += r * F; mres[f] -= r * F;
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, w = jW[f], w1 = 1 - w;
      const cr = (w * dR[L] + w1 * dR[f]) * jG[f], co = (w * dP[L] + w1 * dP[f]) * jG[f];
      const dx = xc[f] - xc[L], dy = yc[f] - yc[L];
      const F = (w * u[L] + w1 * u[f]) * jSx[f] + (w * v[L] + w1 * v[f]) * jSy[f]
        - cr * ((p[f] - p[L]) - ((w * gpx[L] + w1 * gpx[f]) * dx + (w * gpy[L] + w1 * gpy[f]) * dy));
      const r = faceRho(L, f, F, jFx[f], jFy[f]);
      Fj[f] = F; cjF[f] = co; rfJ[f] = r; Mj[f] = r * F; mres[L] += r * F; mres[f] -= r * F;
    }
    // far field: imposed velocity on inflow, imposed pressure (Rhie–Chow towards p_B) on subsonic outflow
    const fco = new Float64Array(NI);
    for (let i = 0; i < NI; i++) {
      const f = NC * NI + i, L = (NC - 1) * NI + i;
      let F, r;
      if (inF[i]) { F = ubF[i] * jSx[f] + vbF[i] * jSy[f]; r = rhoB[i]; if (pDir[i]) fco[i] = 0; }
      else {
        F = u[L] * jSx[f] + v[L] * jSy[f]; r = rho[L];
        if (pDir[i]) {
          const ex = jFx[f] - xc[L], ey = jFy[f] - yc[L];
          F -= dR[L] * jG[f] * ((pB[i] - p[L]) - (gpx[L] * ex + gpy[L] * ey));
          fco[i] = dP[L] * jG[f];
        }
      }
      Fj[f] = F; rfJ[f] = r; Mj[f] = r * F; mres[L] += r * F;
    }
    // time derivative of the density
    if (tr) for (let c = 0; c < N; c++) {
      const k = vol[c] / tr.dt;
      mres[c] += tr.n > 0 ? k * (1.5 * rho[c] - 2 * tr.r1[c] + 0.5 * tr.r2[c]) : k * (rho[c] - tr.r1[c]);
    }
    // pressure-correction equation
    aP.fill(0); aW.fill(0); aE.fill(0); aS.fill(0); aN.fill(0);
    for (let f = 0; f < N; f++) {
      const L = nW[f], F = Fi[f], d = rfI[f] * ciF[f], fp = F > 0 ? F : 0, fm = F < 0 ? -F : 0;
      aE[L] = d + psi[f] * fm; aW[f] = d + psi[L] * fp; aP[L] += d + psi[L] * fp; aP[f] += d + psi[f] * fm;
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, F = Fj[f], d = rfJ[f] * cjF[f], fp = F > 0 ? F : 0, fm = F < 0 ? -F : 0;
      aN[L] = d + psi[f] * fm; aS[f] = d + psi[L] * fp; aP[L] += d + psi[L] * fp; aP[f] += d + psi[f] * fm;
    }
    let anyDir = false;
    for (let i = 0; i < NI; i++) {
      const f = NC * NI + i, L = (NC - 1) * NI + i, F = Fj[f];
      if (!inF[i]) aP[L] += rfJ[f] * fco[i] + psi[L] * F;   // outflow: density and (subsonic) pressure-driven flux respond
      if (pDir[i]) anyDir = true;
    }
    if (tr) for (let c = 0; c < N; c++) aP[c] += (tr.n > 0 ? 1.5 : 1) * vol[c] / tr.dt * psi[c];
    let rc = 0;
    for (let c = 0; c < N; c++) { rc += Math.abs(mres[c]); bu[c] = -mres[c]; pc[c] = 0; }
    if (useMG) mg.solve(mgA, bu, pc, P.mgCycles || 2);
    else sweep(bu, pc, P.pSweeps || (NI >= 384 ? 8 : NI >= 256 ? 6 : 4));
    if (!anyDir && !tr) { const ref = pc[(NC - 1) * NI]; for (let c = 0; c < N; c++) pc[c] -= ref; }
    for (let f = 0; f < N; f++) {
      const L = nW[f], F = Fi[f], dp = pc[f] - pc[nW[f]];
      Fi[f] = F - ciF[f] * dp;
      Mi[f] = rfI[f] * Fi[f] + F * (F > 0 ? psi[L] * pc[L] : psi[f] * pc[f]);
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, F = Fj[f], dp = pc[f] - pc[L];
      Fj[f] = F - cjF[f] * dp;
      Mj[f] = rfJ[f] * Fj[f] + F * (F > 0 ? psi[L] * pc[L] : psi[f] * pc[f]);
    }
    for (let i = 0; i < NI; i++) {
      const f = NC * NI + i, L = (NC - 1) * NI + i, F = Fj[f];
      if (inF[i]) continue;
      Fj[f] = F + fco[i] * pc[L];
      Mj[f] = rfJ[f] * Fj[f] + F * psi[L] * pc[L];
    }
    grad(pc, 1, zeroNI, gtx, gty, pDir);
    for (let c = 0; c < N; c++) {
      u[c] -= dP[c] * gtx[c]; v[c] -= dP[c] * gty[c]; p[c] += relax.p * pc[c];
    }
    return rc;
  }
  const zeroNI = new Float64Array(NI);

  function turbulence() {
    grad(u, 0, ubF, gux, guy); grad(v, 0, vbF, gvx, gvy);
    const ntb = new Float64Array(NI).fill(ntInf);
    grad(nt, 0, ntb, gtx, gty);
    aP.fill(0); aW.fill(0); aE.fill(0); aS.fill(0); aN.fill(0); bu.fill(0);
    // compressible: d(ρν̃)/dt + div(ρu ν̃) = div((μ + ρν̃)/σ grad ν̃) + ρ (sources); gam holds μ + ρν̃
    if (comp) for (let c = 0; c < N; c++) gam[c] = muL[c] + rho[c] * (nt[c] > 0 ? nt[c] : 0);
    const Dc = (L, Rc, w) => comp ? (w * gam[L] + (1 - w) * gam[Rc]) : (nu + w * nt[L] + (1 - w) * nt[Rc]);
    for (let f = 0; f < N; f++) {
      const L = nW[f], w = iW[f], F = Mi[f];
      const D = Dc(L, f, w) / SIG * iG[f];
      const aLR = D + (F > 0 ? F : 0), aRL = D + (F < 0 ? -F : 0);
      aE[L] = aRL; aW[f] = aLR; aP[L] += aLR; aP[f] += aRL;
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, w = jW[f], F = Mj[f];
      const D = Dc(L, f, w) / SIG * jG[f];
      const aLR = D + (F > 0 ? F : 0), aRL = D + (F < 0 ? -F : 0);
      aN[L] = aRL; aS[f] = aLR; aP[L] += aLR; aP[f] += aRL;
    }
    for (let i = 0; i < NI; i++) {
      aP[i] += (comp ? muL[i] : nu) / SIG * jG[i];
      const f = NC * NI + i, L = (NC - 1) * NI + i, F = Mj[f];
      if (inF[i]) {
        const D = (comp ? muL[L] + rhoB[i] * ntInf : nu + ntInf) / SIG * jG[f];
        aP[L] += D + (F > 0 ? F : 0); bu[L] += (D + (F < 0 ? -F : 0)) * ntInf;
      } else if (F > 0) aP[L] += F; else bu[L] -= F * nt[L];
    }
    for (let c = 0; c < N; c++) {
      const rc = comp ? rho[c] : 1, nuc = comp ? muL[c] / rho[c] : nu;
      const n0 = nt[c] > 0 ? nt[c] : 0, chi = n0 / nuc, f1 = fv1(chi), fv2 = 1 - chi / (1 + chi * f1);
      const d = wd[c], d2 = d * d, kd = 1 / (KAP * KAP * d2);
      const Om = Math.abs(gvx[c] - guy[c]);
      let St = Om + n0 * fv2 * kd;
      if (St < 0.3 * Om) St = 0.3 * Om;
      let r = St > 1e-14 ? n0 * kd / St : 10;
      if (r > 10) r = 10;
      const g = r + CW2 * (r ** 6 - r), fw = g * Math.pow((1 + CW3_6) / (g ** 6 + CW3_6), 1 / 6);
      if (comp) {
        bu[c] += (CB1 * St * n0 + CB2 / SIG * (gtx[c] * gtx[c] + gty[c] * gty[c])) * vol[c] * rc;
        aP[c] += CW1 * fw * n0 / d2 * vol[c] * rc;
      } else {
        bu[c] += (CB1 * St * n0 + CB2 / SIG * (gtx[c] * gtx[c] + gty[c] * gty[c])) * vol[c];
        aP[c] += CW1 * fw * n0 / d2 * vol[c];
      }
      if (tr) {
        const k = vol[c] / tr.dt * rc;
        if (tr.n > 0) { aP[c] += 1.5 * k; bu[c] += k * (2 * tr.t1[c] - 0.5 * tr.t2[c]); }
        else { aP[c] += k; bu[c] += k * tr.t1[c]; }
      }
      const apr = aP[c] / relax.t;
      bu[c] += (apr - aP[c]) * nt[c]; aP[c] = apr;
    }
    sweep(bu, nt, 1);
    if (comp) for (let c = 0; c < N; c++) { if (!(nt[c] > 0)) nt[c] = 0; mt[c] = nt[c] * fv1(nt[c] * rho[c] / muL[c]); }
    else for (let c = 0; c < N; c++) { if (!(nt[c] > 0)) nt[c] = 0; mt[c] = nt[c] * fv1(nt[c] / nu); }
  }

  /* Assemble convection (first-order upwind) and diffusion for a cell scalar with
     per-cell diffusivity gam. Wall: Dirichlet value wallVal[i] (array) with wall
     diffusivity nu, or zero flux when wallVal is null. Far field: inflow value inf. */
  function scalarCoeffs(phi, wallVal, inf) {
    aP.fill(0); aW.fill(0); aE.fill(0); aS.fill(0); aN.fill(0); bu.fill(0);
    for (let f = 0; f < N; f++) {
      const L = nW[f], w = iW[f], F = Mi[f];
      const D = (w * gam[L] + (1 - w) * gam[f]) * iG[f];
      const aLR = D + (F > 0 ? F : 0), aRL = D + (F < 0 ? -F : 0);
      aE[L] = aRL; aW[f] = aLR; aP[L] += aLR; aP[f] += aRL;
    }
    for (let f = NI; f < NC * NI; f++) {
      const L = f - NI, w = jW[f], F = Mj[f];
      const D = (w * gam[L] + (1 - w) * gam[f]) * jG[f];
      const aLR = D + (F > 0 ? F : 0), aRL = D + (F < 0 ? -F : 0);
      aN[L] = aRL; aS[f] = aLR; aP[L] += aLR; aP[f] += aRL;
    }
    for (let i = 0; i < NI; i++) {
      if (wallVal) { const D = (comp ? muL[i] : nu) * jG[i]; aP[i] += D; bu[i] += D * wallVal[i]; }
      const f = NC * NI + i, L = (NC - 1) * NI + i, F = Mj[f];
      if (inF[i]) {
        const D = gam[L] * jG[f];
        aP[L] += D + (F > 0 ? F : 0); bu[L] += (D + (F < 0 ? -F : 0)) * inf;
      } else if (F > 0) aP[L] += F; else bu[L] -= F * phi[L];
    }
  }
  function timeAndRelax(phi, h1, h2) {
    for (let c = 0; c < N; c++) {
      if (tr) {
        const k = comp ? vol[c] / tr.dt * rho[c] : vol[c] / tr.dt;
        if (tr.n > 0) { aP[c] += 1.5 * k; bu[c] += k * (2 * h1[c] - 0.5 * h2[c]); }
        else { aP[c] += k; bu[c] += k * h1[c]; }
      }
      const apr = aP[c] / relax.t;
      bu[c] += (apr - aP[c]) * phi[c]; aP[c] = apr;
    }
  }
  const kZero = new Float64Array(NI), kFar = new Float64Array(NI), wFar = new Float64Array(NI);
  function turbulenceSST() {
    const C = SST;
    grad(u, 0, ubF, gux, guy); grad(v, 0, vbF, gvx, gvy);
    kFar.fill(kInf); wFar.fill(wInf);
    grad(tk, 0, kFar, gkx, gky); grad(tw, 1, wFar, gwx, gwy);
    // blending functions and strain rate with the current k and omega
    for (let c = 0; c < N; c++) {
      const k = tk[c] > 0 ? tk[c] : 0, w = tw[c], d = wd[c], d2 = d * d, nuc = comp ? muL[c] / rho[c] : nu;
      S2a[c] = 2 * (gux[c] * gux[c] + gvy[c] * gvy[c]) + (guy[c] + gvx[c]) ** 2;
      const CD = Math.max(2 * C.sw2 / w * (gkx[c] * gwx[c] + gky[c] * gwy[c]), 1e-20);
      const sk = Math.sqrt(k), a = Math.max(sk / (C.bs * w * d), 500 * nuc / (d2 * w));
      const arg1 = Math.min(a, 4 * C.sw2 * k / (CD * d2));
      F1a[c] = Math.tanh(arg1 ** 4);
      if (trans) { const Ry = d * sk / nuc, F3 = Math.exp(-((Ry / 120) ** 8)); if (F3 > F1a[c]) F1a[c] = F3; }   // keeps k-ω in the laminar boundary layer
    }
    // Start-up: k and ω run, but turbulence production stays off (γ held at 0 everywhere) for the first
    // 300 iterations, until the impulsive start has settled and ω has built its boundary-layer profile.
    // Otherwise the start-up transient trips the boundary layer and the turbulent state then sustains
    // itself (the model is bistable). Then the γ equation is released from γ = 0.02.
    // hold: the floor ramps out over the first tRamp iterations, then the boundary layer relaxes to its
    // laminar thickness (no floor) until γ is released at t0h + th
    const t0h = P.tFreeze ?? 0, th = P.tHold ?? 500, tRamp = P.tRamp ?? 300, hold = trans && !warm && S.it < t0h + th;
    holdFrac = hold ? Math.max(0, 1 - Math.max(0, S.it - t0h) / tRamp) : 0; holdOn = hold;
    if (hold) { const dh = P.gHoldD ?? 1e9; for (let c = 0; c < N; c++) if (wd[c] < dh) gmi[c] = 0; }
    else if (trans) { if (!gReleased) { gReleased = true; for (let c = 0; c < N; c++) if (gmi[c] < 0.02) gmi[c] = 0.02; } gammaEq(); }
    // k equation
    if (comp) for (let c = 0; c < N; c++) { const f1 = F1a[c]; gam[c] = muL[c] + (f1 * C.sk1 + (1 - f1) * C.sk2) * rho[c] * mt[c]; }
    else for (let c = 0; c < N; c++) { const f1 = F1a[c]; gam[c] = nu + (f1 * C.sk1 + (1 - f1) * C.sk2) * mt[c]; }
    scalarCoeffs(tk, kZero, kInf);
    for (let c = 0; c < N; c++) {
      const k = tk[c] > 0 ? tk[c] : 0, w = tw[c];
      const Pk0 = Math.min((KLp ? Math.sqrt(S2a[c]) * Math.abs(gvx[c] - guy[c]) : S2a[c]) * mt[c], 10 * C.bs * k * w), rv = comp ? rho[c] * vol[c] : vol[c];
      const Pk = Pk0;
      if (trans) {
        // intermittency scales production; destruction keeps at least 10 %; Pk_lim lets separation-induced transition happen quickly
        const g = gmi[c], nuc = comp ? muL[c] / rho[c] : nu, Fl = Math.min(Math.max(ReVa[c] / (2.2 * 1100) - 1, 0), 3);
        const Plim = 5 * Math.max(g - 0.2, 0) * (1 - g) * Fl * Math.max(3 * nuc - mt[c], 0) * Math.sqrt(S2a[c]) * Math.abs(gvx[c] - guy[c]);
        bu[c] += (g * Pk + Plim + C.bs * kInf * wInf) * rv;
        aP[c] += Math.max(g, 0.1) * C.bs * w * rv;
      } else {
        bu[c] += (Pk + C.bs * kInf * wInf) * rv;
        aP[c] += C.bs * w * rv;
      }
    }
    timeAndRelax(tk, tr && tr.k1, tr && tr.k2);
    sweep(bu, tk, 1);
    for (let c = 0; c < N; c++) if (!(tk[c] > 1e-16)) tk[c] = 1e-16;
    // omega equation
    if (comp) for (let c = 0; c < N; c++) { const f1 = F1a[c]; gam[c] = muL[c] + (f1 * C.sw1 + (1 - f1) * C.sw2) * rho[c] * mt[c]; }
    else for (let c = 0; c < N; c++) { const f1 = F1a[c]; gam[c] = nu + (f1 * C.sw1 + (1 - f1) * C.sw2) * mt[c]; }
    scalarCoeffs(tw, null, wInf);
    for (let c = 0; c < N; c++) {
      const f1 = F1a[c], k = tk[c], w = tw[c];
      const g = f1 * C.g1 + (1 - f1) * C.g2, b = f1 * C.b1 + (1 - f1) * C.b2;
      const Pw = g * (mt[c] > 1e-30 ? Math.min(S2a[c], 10 * C.bs * k * w / mt[c]) : S2a[c]), rv = comp ? rho[c] * vol[c] : vol[c];
      bu[c] += (Pw + b * w * w + b * wInf * wInf) * rv;   // destruction linearised: 2 b w (implicit) - b w^2
      aP[c] += 2 * b * w * rv;
      const cd = 2 * (1 - f1) * C.sw2 / w * (gkx[c] * gwx[c] + gky[c] * gwy[c]);
      if (cd > 0) bu[c] += cd * rv; else aP[c] -= cd / w * rv;
    }
    timeAndRelax(tw, tr && tr.w1, tr && tr.w2);
    if (comp) for (let i = 0; i < NI; i++) wWall[i] = 6 * muL[i] / rho[i] / (SST.b1 * wd[i] * wd[i]);
    for (let i = 0; i < NI; i++) { aP[i] = 1e30; bu[i] = 1e30 * wWall[i]; }
    sweep(bu, tw, 1);
    for (let c = 0; c < N; c++) if (!(tw[c] > 1e-8)) tw[c] = 1e-8;
    // eddy viscosity with the SST limiter
    for (let c = 0; c < N; c++) {
      const k = tk[c], w = tw[c], d = wd[c], d2 = d * d;
      const arg2 = Math.max(2 * Math.sqrt(k) / (C.bs * w * d), 500 * (comp ? muL[c] / rho[c] : nu) / (d2 * w)), F2 = Math.tanh(arg2 * arg2);
      mt[c] = C.a1 * k / Math.max(C.a1 * w, Math.sqrt(S2a[c]) * F2);
    }
    // During the start-up hold the flow would be laminar at full Reynolds number, which separates and
    // goes unsteady at incidence. A viscosity floor, ramped from 50 ν to zero over the hold, keeps it steady.
    // Near the wall (d < 0.01 c) the floor ramps out first, so the boundary layer reaches its laminar thickness
    // before γ is released; the wake and outer flow keep it for the whole hold.
    if (trans && holdOn) {
      const fw = (P.holdVisc ?? 50) * holdFrac * nu, fo = (P.holdVisc ?? 50) * nu, dF = P.holdD ?? 0.01;
      for (let c = 0; c < N; c++) { const fl = wd[c] < dF ? fw : fo; if (mt[c] < fl) mt[c] = fl; }
    }
  }

  /* Intermittency γ (Menter et al. 2015): d(ργ)/dt + div(ρuγ) = Pγ − Eγ + div((μ + μt) grad γ),
     Pγ = 100 ρ S γ(1 − γ) Fonset, Eγ = 0.06 ρ Ω γ Fturb (50γ − 1). Onset where the vorticity Reynolds
     number d²S/ν exceeds 2.2 Reθc, with Reθc from the local turbulence intensity and a pressure-gradient
     parameter built from the wall-normal derivative of the wall-normal velocity. Zero flux at the wall,
     γ = 1 in the inflow. */
  function gammaEq() {
    if (!nwDone) {   // wall-normal direction = gradient of the wall distance
      grad(wd, 0, null, nwx, nwy);
      for (let c = 0; c < N; c++) { const l = Math.hypot(nwx[c], nwy[c]) || 1; nwx[c] /= l; nwy[c] /= l; }
      nwDone = true;
    }
    for (let c = 0; c < N; c++) vnA[c] = u[c] * nwx[c] + v[c] * nwy[c];
    grad(vnA, 0, null, vnGx, vnGy);
    for (let c = 0; c < N; c++) gam[c] = comp ? muL[c] + rho[c] * mt[c] : nu + mt[c];
    scalarCoeffs(gmi, null, 1);
    for (let c = 0; c < N; c++) {
      const k = tk[c] > 0 ? tk[c] : 0, w = tw[c], d = wd[c], nuc = comp ? muL[c] / rho[c] : nu, rv = comp ? rho[c] * vol[c] : vol[c];
      const Sm = Math.sqrt(S2a[c]), Om = Math.abs(gvx[c] - guy[c]);
      const ReV = d * d * Sm / nuc, RT = k / (nuc * w);
      ReVa[c] = ReV;
      const TuL = Math.min(100 * Math.sqrt(2 * k / 3) / (w * d), 100);
      const dVdy = vnGx[c] * nwx[c] + vnGy[c] * nwy[c];
      const lam = Math.min(Math.max(-7.57e-3 * dVdy * d * d / nuc + 0.0128, -1), 1);
      const FPG = Math.max(lam >= 0 ? Math.min(1 + 14.68 * lam, 1.5) : Math.min(1 - 7.34 * lam, 3), 0);
      const Rec = 100 + 1000 * Math.exp(-TuL * FPG);
      const Fon = Math.max(Math.min(ReV / (2.2 * Rec), 2) - Math.max(1 - (RT / 3.5) ** 3, 0), 0);
      const Ftu = Math.exp(-((RT / 2) ** 4)), g = gmi[c];
      if (DBGT) { DBGT.TuL[c] = TuL; DBGT.Rec[c] = Rec; DBGT.Fon[c] = Fon; DBGT.RT[c] = RT; DBGT.ReV[c] = ReV; DBGT.lam[c] = lam; }
      const P0 = 100 * Sm * Fon * g;          // Pγ = P0 (1 − γ): explicit P0, implicit P0 γ
      const E0 = 0.06 * Om * Ftu * g;         // Eγ = E0 (50γ − 1): implicit 50 E0 γ, explicit + E0
      bu[c] += (P0 + E0) * rv;
      aP[c] += (P0 + 50 * E0) * rv;
    }
    timeAndRelax(gmi, tr && tr.g1, tr && tr.g2);
    sweep(bu, gmi, 1);
    for (let c = 0; c < N; c++) { const g = gmi[c]; gmi[c] = g > 1 ? 1 : g > 0 ? g : 0; }
  }
  const vnGx = trans ? F64() : null, vnGy = trans ? F64() : null;

  const cpW = new Float64Array(NI), cfW = new Float64Array(NI), ypW = new Float64Array(NI);
  function forces() {
    let pinf = 0;
    if (comp) pinf = 0;   // the far field fixes the pressure level
    else { for (let i = 0; i < NI; i++) { const c = (NC - 1) * NI + i; pinf += p[c] + 0.5 * (u[c] * u[c] + v[c] * v[c] - 1); } pinf /= NI; }
    S.pinf = pinf;
    let fpx = 0, fpy = 0, fvx = 0, fvy = 0, mz = 0, ypm = 0, yps = 0;
    for (let i = 0; i < NI; i++) {
      const sx = jSx[i], sy = jSy[i], A = Math.hypot(sx, sy), nx = sx / A, ny = sy / A;
      const pw = p[i] - pinf;
      const dn = (xc[i] - jFx[i]) * nx + (yc[i] - jFy[i]) * ny;
      const un = u[i] * nx + v[i] * ny, utx = u[i] - un * nx, uty = v[i] - un * ny;
      const k = inviscid ? 0 : (comp ? muL[i] : nu) / dn;
      const dfx = -pw * sx + k * utx * A, dfy = -pw * sy + k * uty * A;
      fpx -= pw * sx; fpy -= pw * sy; fvx += k * utx * A; fvy += k * uty * A;
      mz += (jFx[i] - 0.25) * dfy - jFy[i] * dfx;
      cpW[i] = pw / 0.5;
      // skin friction along the surface, positive when the wall flow runs towards the trailing edge
      const ip = (i + 1) % NI, tl = Math.hypot(X[ip] - X[i], Y[ip] - Y[i]), ts = i < NI / 2 ? -1 : 1;
      cfW[i] = k * (utx * (X[ip] - X[i]) + uty * (Y[ip] - Y[i])) * ts / tl / 0.5;
      const yp = inviscid ? 0 : comp ? dn * Math.sqrt(k * Math.hypot(utx, uty) * rho[i]) / muL[i] : dn * Math.sqrt(k * Math.hypot(utx, uty)) / nu;
      ypW[i] = yp; if (yp > ypm) ypm = yp; yps += yp;
    }
    const ca = S.ca, sa = S.sa;
    S.cdp = 2 * (fpx * ca + fpy * sa); S.cdf = 2 * (fvx * ca + fvy * sa);
    S.cd = S.cdp + S.cdf;
    S.cl = 2 * ((fpx + fvx) * -sa + (fpy + fvy) * ca);
    S.cm = -2 * mz;
    S.ypMax = ypm; S.ypMean = yps / NI;
  }

  /* Wave drag from the entropy the shocks produce (Oswatitsch): D = T∞ Σ ṁ Δs over the
     boundary of the shock region, i.e. C_D,wave = 2/(γM²) Σ ṁ Δs/R. Shock cells are found with
     the normal Mach number u·∇p / (a|∇p|) > 1 (Lovely & Haimes) where the pressure rises
     steeply, outside the boundary layer and wake (νt/ν < 20), then widened by two cells. */
  const shockMark = comp ? new Uint8Array(N) : null, shockTmp = comp ? new Uint8Array(N) : null, sEnt = comp ? F64() : null;
  function waveDrag() {
    shockMark.fill(0);
    let n = 0;
    for (let c = 0; c < N; c++) {
      const gp = Math.hypot(gpx[c], gpy[c]); if (gp < 1e-12) continue;
      const a = Math.sqrt(Tc[c]) / Mr, Mn = (u[c] * gpx[c] + v[c] * gpy[c]) / (a * gp);
      const h = Math.sqrt(vol[c]);
      if (Mn > 1 && gp * h > 0.02 && (inviscid || mt[c] * rho[c] < 20 * muL[c])) { shockMark[c] = 1; n++; }
    }
    S.shockCells = n;
    if (!n) { S.cdw = 0; return; }
    const outer = c => inviscid || mt[c] * rho[c] < 20 * muL[c];
    for (let pass = 0; pass < (P.shockGrow || 8); pass++) {
      shockTmp.set(shockMark);
      for (let c = 0; c < N; c++) if (shockTmp[c]) for (const q of [nW[c], nE[c], nS[c], nN[c]]) if (q >= 0 && outer(q)) shockMark[q] = 1;
    }
    for (let c = 0; c < N; c++) sEnt[c] = GAM / (GAM - 1) * Math.log(Tc[c]) - Math.log(Math.max(1e-6, 1 + GAM * M2r * p[c]));
    let out = 0;
    const fl = (L, Rc, m) => { const a = shockMark[L], b = shockMark[Rc]; if (a === b) return; const sf = m > 0 ? sEnt[L] : sEnt[Rc]; out += a ? m * sf : -m * sf; };
    for (let f = 0; f < N; f++) fl(nW[f], f, Mi[f]);
    for (let f = NI; f < NC * NI; f++) fl(f - NI, f, Mj[f]);
    // Wave drag is part of the pressure drag. Where the shock region borders the wake (trailing-edge
    // shocks), viscous entropy leaks into the sum, so the estimate is capped at the pressure drag.
    S.cdwRaw = 2 / (GAM * M2r) * out; S.cdw = Math.min(Math.max(0, S.cdwRaw), Math.max(0, S.cdp));
  }

  // total enthalpy (URANS only: in a steady flow with an adiabatic wall and Prandtl number 1 it is constant)
  function energy() {
    for (let c = 0; c < N; c++) gam[c] = me[c];
    scalarCoeffs(hh, null, 0);
    for (let c = 0; c < N; c++) {
      const k = vol[c] / tr.dt;
      if (tr.n > 0) { aP[c] += 1.5 * k * rho[c]; bu[c] += k * (2 * tr.r1[c] * tr.h1[c] - 0.5 * tr.r2[c] * tr.h2[c] + 1.5 * p[c] - 2 * tr.p1[c] + 0.5 * tr.p2[c]); }
      else { aP[c] += k * rho[c]; bu[c] += k * (tr.r1[c] * tr.h1[c] + p[c] - tr.p1[c]); }
      const apr = aP[c] / 0.9; bu[c] += (apr - aP[c]) * hh[c]; aP[c] = apr;
    }
    sweep(bu, hh, 1);
  }

  const lam = P.laminar || inviscid;
  function iterate() {
    grad(u, 0, ubF, gux, guy); grad(v, 0, vbF, gvx, gvy);
    if (comp) { rampMach(); grad(p, 1, pB, gpx, gpy, pDir); grad(rho, 1, rhoB, grx, gry); }
    else grad(p, 1, null, gpx, gpy);
    const ru = momentum();
    const rc = comp ? pressureComp() : pressure();
    // transition, optional first start-up phase: turbulence frozen (off by default; the hold below suffices)
    const frozen = trans && S.it < (P.tFreeze ?? 0) && !warm;
    if (!lam && !frozen) (sst ? turbulenceSST : turbulence)(); else { grad(u, 0, ubF, gux, guy); grad(v, 0, vbF, gvx, gvy); }
    if (comp) { if (tr) energy(); thermo(false); }
    forces();
    if (comp && S.it % 10 === 0) waveDrag();
    S.it++;
    if (S.it % 10 === 0) farBC();
    if (!S.res0) S.res0 = { c: rc || 1, u: ru || 1 };
    S.resC = rc / S.res0.c; S.resU = ru / S.res0.u; S.rcRaw = rc;
  }

  /* ---------- unsteady (URANS) ----------
     Dual time stepping: each physical step dt (in chord / U units) is converged with
     inner SIMPLEC iterations; the BDF2 time term is added to momentum and SA. */
  let tr = null;
  // Dual time stepping: the time term already stabilises the inner iterations, and 0.85
  // needs fewer of them per step than the steady-state 0.95.
  let relaxSteady = relax.u;
  function startUnsteady(dt) {
    relax.u = Math.min(relaxSteady, P.relaxUt || 0.85);
    tr = { dt, n: 0, time: 0, u1: Float64Array.from(u), v1: Float64Array.from(v), t1: Float64Array.from(nt),
           u2: Float64Array.from(u), v2: Float64Array.from(v), t2: Float64Array.from(nt) };
    if (sst) { tr.k1 = Float64Array.from(tk); tr.k2 = Float64Array.from(tk); tr.w1 = Float64Array.from(tw); tr.w2 = Float64Array.from(tw); }
    if (trans) { tr.g1 = Float64Array.from(gmi); tr.g2 = Float64Array.from(gmi); }
    if (comp) for (const [a, b] of [["r", rho], ["h", hh], ["p", p]]) { tr[a + "1"] = Float64Array.from(b); tr[a + "2"] = Float64Array.from(b); }
  }
  function timeStep(maxInner) {
    let rc0 = 0, k = 1;
    for (; k <= maxInner; k++) {
      iterate();
      if (k === 1) rc0 = S.rcRaw;
      else if (k >= 3 && S.rcRaw < 0.05 * rc0) break;   // mass imbalance down by ~1.3 orders: step converged
    }
    tr.u2.set(tr.u1); tr.v2.set(tr.v1); tr.t2.set(tr.t1);
    tr.u1.set(u); tr.v1.set(v); tr.t1.set(nt);
    if (sst) { tr.k2.set(tr.k1); tr.w2.set(tr.w1); tr.k1.set(tk); tr.w1.set(tw); }
    if (trans) { tr.g2.set(tr.g1); tr.g1.set(gmi); }
    if (comp) { tr.r2.set(tr.r1); tr.h2.set(tr.h1); tr.p2.set(tr.p1); tr.r1.set(rho); tr.h1.set(hh); tr.p1.set(p); }
    tr.n++; tr.time += tr.dt;
    farBC();
    return Math.min(k, maxInner);
  }

  function fields() {
    const sp = new Float32Array(N), ux = new Float32Array(N), cp = new Float32Array(N), nr = new Float32Array(N), vo = new Float32Array(N);
    grad(u, 0, ubF, gux, guy); grad(v, 0, vbF, gvx, gvy);
    for (let c = 0; c < N; c++) {
      sp[c] = Math.hypot(u[c], v[c]); ux[c] = u[c] * S.ca + v[c] * S.sa; cp[c] = (p[c] - S.pinf) / 0.5; nr[c] = mt[c] / nu; vo[c] = gvx[c] - guy[c];
    }
    const n = NI / 2, sx = new Float32Array(NI), upper = new Uint8Array(NI);
    for (let i = 0; i < NI; i++) { sx[i] = jFx[i]; upper[i] = i < HB || i >= NI - HB ? 2 : i < n ? 1 : 0; }   // 2 = TE base
    let mach = null, schl = null;
    if (comp) {
      mach = new Float32Array(N); schl = new Float32Array(N);
      grad(rho, 1, rhoB, grx, gry);
      for (let c = 0; c < N; c++) { mach[c] = sp[c] * Mr / Math.sqrt(Tc[c]); schl[c] = Math.hypot(grx[c], gry[c]); nr[c] = mt[c] * rho[c] / (muL[c] || nu); }
    }
    // intermittency near the wall: the lowest γ in each wall-normal column within 0.01 chord of the wall
    let gamma = null, gw = null;
    if (trans) {
      gamma = Float32Array.from(gmi); gw = new Float32Array(NI);
      // turbulent where the boundary layer carries eddy viscosity: largest νt/ν within 0.02 chord of the wall
      for (let i = 0; i < NI; i++) { let m = 0; for (let j = 0; j < NC; j++) { const c = j * NI + i; if (wd[c] > 0.02) break; const r = mt[c] / (comp ? muL[c] / rho[c] : nu); if (r > m) m = r; } gw[i] = m; }
    }
    return { mach, schl, gamma, speed: sp, ux, u: Float32Array.from(u), v: Float32Array.from(v), cp, nutRatio: nr, vort: vo, wall: { x: sx, upper, cp: Float32Array.from(cpW), cf: Float32Array.from(cfW), yp: Float32Array.from(ypW), gw } };
  }

  // continue from a solution interpolated from another mesh of the same case
  function loadState(st, from) {
    u.set(st.u); v.set(st.v); p.set(st.p); nt.set(st.nt); mt.set(st.mt);
    if (sst) { tk.set(st.tk); tw.set(st.tw); }
    if (trans) { if (st.gm) gmi.set(st.gm); else gmi.fill(1); warm = true; }
    if (comp && st.hh) hh.set(st.hh);
    S.cl = from.cl; S.cd = from.cd; S.it = from.it; S.res0 = from.res0; farBC();
    if (comp) { M2r = -1; rampMach(); thermo(true); }
    for (let f = 0; f < N; f++) {
      const L = nW[f], w = iW[f];
      Fi[f] = (w * u[L] + (1 - w) * u[f]) * iSx[f] + (w * v[L] + (1 - w) * v[f]) * iSy[f];
    }
    for (let f = 0; f < NF; f++) {
      if (f < NI) { Fj[f] = 0; continue; }
      const L = f - NI;
      if (f >= NC * NI) { Fj[f] = u[L] * jSx[f] + v[L] * jSy[f]; continue; }
      const w = jW[f];
      Fj[f] = (w * u[L] + (1 - w) * u[f]) * jSx[f] + (w * v[L] + (1 - w) * v[f]) * jSy[f];
    }
    if (comp) massFromVolume();
  }

  setAlpha(P.alpha || 0);
  initFlow();
  return {
    loadState, d: G.d, relaxSteady: () => relaxSteady,
    NI, NJ, N, X, Y, ds1, ratio: G.ratio, hb: HB, badCells, S, relax,
    iterate, setAlpha, setGust: deg => { S.gust = deg * Math.PI / 180; farBC(); }, fields, initFlow, xc, yc, vol, wd, mt, startUnsteady, timeStep,
    unsteady: () => tr, stopUnsteady: () => { tr = null; relax.u = relaxSteady; }, lowerSteadyRelax: v => { relaxSteady = v || 0.85; },
    state: () => ({ u, v, p, nt, tk, tw, hh, gm: gmi }), DBGT, turb: sst ? "sst" : "sa", comp, mach: MACH
  };
}

/* Move a solution to another mesh of the same case (both given by computational node coordinates) */
function transferState(Co, compO, Cn, compN) {
  const NI = Cn.NI, N = Cn.N;
  const co = cellComp(NI, Co.NJ, compO.XI, compO.ET), cn = cellComp(NI, Cn.NJ, compN.XI, compN.ET);
  const loc = makeLocator(NI, Co.NJ - 1, co.A, co.B), so = Co.state(), src = { u: so.u, v: so.v, p: so.p, nt: so.nt, mt: Co.mt };
  if (so.tk) { src.tk = so.tk; src.tw = so.tw; }
  if (so.gm) src.gm = so.gm;
  if (so.hh) src.hh = so.hh;
  const out = {}, keys = Object.keys(src);
  for (const k of keys) out[k] = new Float64Array(N);
  for (let c = 0; c < N; c++) { const q = loc(cn.A[c], cn.B[c]); for (const k of keys) out[k][c] = loc.val(src[k], q); }
  Cn.loadState(out, Co.S);
}

/* Adapt the mesh of case C to its current solution. AD: {P, base, cur}. Returns the new case or null. */
function adaptCase(C, AD) {
  const st = C.state(); st.comp = C.comp; st.mach = C.mach;
  let res = null;
  AD.base.round = (sharpSection(AD.P.code) || {}).round === true;
  for (const str of [1, 0.6, 0.35]) { res = adaptGrid(AD.base, AD.cur, st, str); if (!res.bad) break; }
  if (res.bad) return null;
  const B = AD.base;
  const Cn = createCase({ ...AD.P, NJ: res.NJ, alpha: C.S.alpha, grid: { X: res.X, Y: res.Y, d: B.d, ratio: B.ratio, hb: B.hb } });
  if (Cn.badCells) return null;
  transferState(C, AD.cur, Cn, res.comp);
  if (C.relax.u < Cn.relax.u) { Cn.relax.u = C.relax.u; Cn.lowerSteadyRelax(C.relax.u); }
  AD.cur = res.comp;
  Cn.adaptStats = res.stats;
  return Cn;
}

/* ---------- worker protocol ---------- */
if (typeof self !== "undefined" && typeof window === "undefined") {
  let C = null, running = false, hist = [], token = 0, maxIt = 8000, tol = 1;   // tol scales the force-change limits
  let mode = "steady", U = null;   // U: unsteady run settings
  // adaptive mesh: passes per solution, passes left, iteration where the current pass started
  let AD = null, pendU = null;
  const ADAPT_CAP = 1500;   // an intermediate solution only has to show where the shear layers are
  const gridMsg = (t0, extra) => ({ type: "grid", NI: C.NI, NJ: C.NJ, X: Float32Array.from(C.X), Y: Float32Array.from(C.Y), ds1: C.ds1, ratio: C.ratio, hb: C.hb, badCells: C.badCells, ms: performance.now() - t0, ...extra });
  function startU(m) {
    if (m.fromRest) C.initFlow();
    mode = "urans"; U = { gust: m.gust !== false, tEnd: m.tEnd, maxInner: m.maxInner || 15, fieldEvery: m.fieldEvery || 8, lastField: 0, recEvery: m.recEvery || 0, recFrom: m.recFrom || 0 };
    C.startUnsteady(m.dt);
    running = true; token++; loopU(token);
  }
  const post = (m, t) => self.postMessage(m, t || []);
  function snapshot() {
    const S = C.S;
    return { it: S.it, cl: S.cl, cd: S.cd, cdp: S.cdp, cdf: S.cdf, cm: S.cm, resC: S.resC, resU: S.resU, ypMax: S.ypMax, ypMean: S.ypMean, alpha: S.alpha, cdw: S.cdw };
  }
  function converged() {
    const n = hist.length;
    if (n < 30 || C.S.it < 400) return false;
    // compare the average of the last 100 iterations with the 100 before (sampled every 10):
    // averaging ignores small decaying wiggles but still catches a real drift
    let la = 0, lb = 0, da = 0, db = 0;
    for (let k = 0; k < 10; k++) { const A = hist[n - 1 - k], B = hist[n - 11 - k]; la += A.cl; lb += B.cl; da += A.cd; db += B.cd; }
    la /= 10; lb /= 10; da /= 10; db /= 10;
    return Math.abs(la - lb) < 2e-4 * tol * Math.max(1, Math.abs(la)) && Math.abs(da - db) < 2e-5 * tol && hist[n - 1].resC < 1e-2;
  }
  // forces averaged over the last 100 iterations, which removes small residual wiggles
  function averaged() {
    const s = snapshot(), n = hist.length;
    if (n < 10) return s;
    for (const k of ["cl", "cd", "cdp", "cdf", "cm"]) { let a = 0; for (let i = n - 10; i < n; i++) a += hist[i][k]; s[k] = a / 10; }
    return s;
  }
  function loop(tk) {
    if (!running || tk !== token) return;
    const t0 = performance.now();
    while (performance.now() - t0 < 40) {
      C.iterate();
      // Relaxation 0.95 is close to the stability limit in deep stall: if the run blows up,
      // restart once from uniform flow with the conservative 0.85.
      const fb = C.comp && C.mach > 0.6 ? 0.55 : 0.85;   // transonic cases already start at 0.7
      if ((!isFinite(C.S.cl) || C.S.resC > 1e8) && C.relax.u > fb + 0.001) {
        C.relax.u = fb; C.lowerSteadyRelax(fb); C.initFlow(); hist = [];
        post({ type: "fallback", relaxU: fb });
        continue;
      }
      if (C.S.it % 10 === 0) hist.push(snapshot());
      if (C.S.it >= maxIt) break;
    }
    post({ type: "progress", s: snapshot(), hist: hist.slice(-1) });
    if (AD && AD.left > 0 && isFinite(C.S.cl) && (converged() || C.S.it - AD.start >= ADAPT_CAP || C.S.it >= maxIt)) {
      const t0 = performance.now(), Cn = adaptCase(C, AD);
      if (Cn) {
        C = Cn; AD.left--; AD.done++; AD.start = C.S.it; hist = [];
        maxIt = Math.max(maxIt, C.S.it + 8000);
        post(gridMsg(t0, { adapted: AD.done, passes: AD.passes, stats: C.adaptStats }));
        post({ type: "fields", f: C.fields(), s: snapshot() });
      } else { AD.left = 0; post({ type: "adaptFailed" }); }
      setTimeout(() => loop(tk), 0);
      return;
    }
    const done = converged() || C.S.it >= maxIt || !isFinite(C.S.cl);
    if (C.S.it % 200 < 15 || done) { const f = C.fields(); post({ type: "fields", f, s: snapshot() }); }
    if (done) {
      running = false; const ok = converged(), div = !isFinite(C.S.cl);
      post({ type: "done", s: ok ? averaged() : snapshot(), converged: ok, diverged: div, next: pendU && !div ? "urans" : null });
      if (pendU && !div) { const m = pendU; pendU = null; if (m.fresh) m.fromRest = true; startU(m); } else pendU = null;
      return;
    }
    setTimeout(() => loop(tk), 0);
  }
  function loopU(tk) {
    if (!running || tk !== token) return;
    const T = C.unsteady(), t0 = performance.now(), recs = [];
    let inner = 0;
    while (performance.now() - t0 < 40 && T.time < U.tEnd - 1e-9) {
      // a short gust with zero mean (inflow angle +1°, then -1°, over 1.5 c/U) disturbs the flow, so a
      // steady state that is unstable to shedding cannot just persist for lack of a disturbance,
      // without pushing the section into stall the way a one-sided gust can
      if (U.gust) C.setGust(T.time < 1.5 ? Math.sin(2 * Math.PI * T.time / 1.5) : 0);
      inner = C.timeStep(U.maxInner);
      const S = C.S;
      recs.push([T.time, S.cl, S.cd, S.cdp, S.cdf, S.cm, inner]);
      if (!isFinite(S.cl)) break;
      // movie: a compact snapshot every recEvery steps inside the recording window
      if (U.recEvery && T.n % U.recEvery === 0 && T.time >= U.recFrom - 1e-9) {
        const f = C.fields(), fr = { type: "frame", t: T.time, cl: S.cl, cd: S.cd, u: f.u, v: f.v, cp: f.cp, nut: f.nutRatio, wcp: f.wall.cp, wcf: f.wall.cf };
        const tl = [fr.u.buffer, fr.v.buffer, fr.cp.buffer, fr.nut.buffer, fr.wcp.buffer, fr.wcf.buffer];
        if (f.mach) { fr.mach = f.mach; fr.schl = f.schl; tl.push(f.mach.buffer, f.schl.buffer); }
        post(fr, tl);
      }
    }
    const s = snapshot(); s.t = T.time; s.step = T.n; s.inner = inner;
    post({ type: "uprogress", recs, s });
    const done = T.time >= U.tEnd - 1e-9 || !isFinite(C.S.cl);
    if (T.n - U.lastField >= U.fieldEvery || done) { U.lastField = T.n; post({ type: "fields", f: C.fields(), s }); }
    if (done) { running = false; post({ type: "udone", s, diverged: !isFinite(C.S.cl) }); return; }
    setTimeout(() => loopU(tk), 0);
  }
  self.onmessage = e => {
    const m = e.data;
    if (m.type === "init") {
      token++; running = false;
      const t0 = performance.now();
      C = createCase(m.p); mode = "steady"; hist = []; maxIt = m.p.maxIt || 8000; tol = m.p.tol || 1; pendU = null;
      const np = m.p.adapt | 0;
      AD = np > 0 ? { P: m.p, passes: np, left: np, done: 0, start: 0, cur: identityComp(C.NI, C.NJ),
                      base: { NI: C.NI, NJ: C.NJ, X: C.X, Y: C.Y, d: C.d, ratio: C.ratio, hb: C.hb } } : null;
      post(gridMsg(t0, {}));
      post({ type: "fields", f: C.fields(), s: snapshot() });
    } else if (m.type === "urans") {
      if (!C) return;
      pendU = null; if (AD) AD.left = 0;
      startU(m);
    } else if (m.type === "run") {
      if (!C) return;
      if (mode === "urans" && C.unsteady() && C.unsteady().time < U.tEnd - 1e-9 && !m.steady) { running = true; token++; loopU(token); return; }
      mode = "steady"; C.stopUnsteady();
      if (typeof m.maxIt === "number") maxIt = C.S.it + m.maxIt;
      if (m.thenUrans) pendU = m.thenUrans;   // a steady (adapted) solution first, then the unsteady run
      else if (m.steady) pendU = null;
      running = true; token++; loop(token);
    } else if (m.type === "pause") { running = false; token++; }
    else if (m.type === "tol") { tol = m.tol || 1; }
    else if (m.type === "reset") {   // discard the flow: the next run starts from uniform flow (the mesh is kept)
      if (!C) return;
      running = false; token++; mode = "steady"; pendU = null; C.stopUnsteady();
      if (typeof m.alpha === "number") C.setAlpha(m.alpha);
      C.initFlow(); hist = [];
      if (AD) { AD.left = AD.passes; AD.done = 0; AD.start = 0; }
      post({ type: "fields", f: C.fields(), s: snapshot() });
    }
    else if (m.type === "alpha") { if (C) { const ch = m.alpha !== C.S.alpha; C.setAlpha(m.alpha); hist = []; if (AD && ch) { AD.left = AD.passes; AD.done = 0; AD.start = C.S.it; } } }
  };
}
