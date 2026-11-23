# Wing Section Solver

A 2D airfoil CFD solver that runs entirely in the browser. It solves the Reynolds-averaged Navier–Stokes equations on a body-fitted C/O grid in a Web Worker and shows the flow field, surface pressure and skin friction, and integrated forces.

## Features

- **Sections:** NACA 4-digit, transonic sections from published coordinates (RAE 2822, NASA SC(2), NLR 7301, NACA 64A), NASA NLF(1)-0416, Eppler 387, sharp supersonic sections, and a circle (bluff body).
- **Turbulence models:** Spalart–Allmaras, k-ω SST, and k-ω SST with Menter's one-equation γ transition model (2015), with a free-stream turbulence setting.
- **Flow:** incompressible and compressible (subsonic and transonic), steady RANS and URANS with movie playback.
- **Mesh:** Coarse to Ultrafine grids and solution-based adaptation.
- **Visualisation:** pressure, velocity, vorticity, eddy viscosity and intermittency fields; streamlines, smoke and pulsed streaklines; Cp and Cf plots with transition markers; video export.

## Running locally

The page loads `rans-worker.js` as a Web Worker, which imports `sections.js`, so it must be served over HTTP rather than opened as a file:

```bash
python3 -m http.server 8765
```

Then open <http://localhost:8765/wing-section-solver.html>.

## Files

| File | Contents |
|---|---|
| `wing-section-solver.html` | User interface, plots and flow visualisation |
| `rans-worker.js` | Grid generation, adaptation and the RANS/URANS solver |
| `sections.js` | Airfoil geometry: analytic families and tabulated coordinates |
| `pocket-wind-tunnel.html` | The earlier single-file prototype |

## Validation (SST + γ, steady)

| Case | Transition x/c (upper / lower) | C_D with transition | C_D fully turbulent |
|---|---|---|---|
| NACA 0012, Re 3 M, α 0° | 0.44 / 0.44 | 0.0057 | 0.0088 |
| NLF(1)-0416, Re 4 M, α 0° | 0.35–0.36 / 0.39–0.42 (NASA ref. 0.41 / 0.52) | 0.0062–0.0064 | 0.0093 |
| NLF(1)-0416, Re 4 M, α 5° | 0.15–0.20 / 0.56–0.58 | 0.0075–0.0082 | 0.0106 |

Transition lands about 0.05–0.1 chord early on these meshes because a turbulent boundary-layer cell sustains itself. The first 500 iterations of a transition run are a laminar start-up hold.

## References

- Menter, Smirnov, Liu & Avancha (2015), *A One-Equation Local Correlation-Based Transition Model*, Flow Turbulence and Combustion. <https://doi.org/10.1007/s10494-015-9622-4>
- Hildebrand et al. (2023), *Recent Progress on RANS-Based Transition Model Verification*, NASA. <https://ntrs.nasa.gov/api/citations/20230010094/downloads/Hildebrand_ETMM14_RANS_final.pdf>
- Somers (1981), NASA TP-1861 (NLF(1)-0416). <https://ntrs.nasa.gov/citations/19810015487>
- McGhee et al. (1988), NASA TM-4062 (Eppler 387). <https://ntrs.nasa.gov/citations/19890001471>
