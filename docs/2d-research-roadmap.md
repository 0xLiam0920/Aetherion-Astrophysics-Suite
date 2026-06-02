# 2D Research Data Roadmap (Dwarf Galaxy & Generic Body Simulation)

Originally lived as a 73-line block comment at the top of
`src/2D/2D-simulation/research_data.hpp`. Moved here so the header stops
looking like missing code and can be planned separately from shipping
work.

## Phase 1 — Core Physics Validation
- Implement custom FITS unit parsing for `rg`, `M` in astropy (export compatibility).
- Export simulation parameters to PRIMARY HDU: `M_BH`, `c_s`, `r_Bondi/rg`, `beta`.
- Add geodesic convergence tests (ISCO radius, frame-dragging angles).
- Compute/export `ENERGY_DRIFT` per timestep for all integrators.

## Phase 2 — Dynamical Diagnostics
- Track stream morphology: `WIDTH(rg)`, `LENGTH(rg)`, density profile.
- Loss-cone statistics: `THETA_LC`, `FEED_RATE(stars/s)`, `J_SCATTER`.
- Angular momentum histogram: `J_MIN`, `J_ISCO` per particle.
- Precession module: `PERIAPSIS_PHI` evolution (fix empty table).

## Phase 3 — Accretion & Observables
- Bondi-Hoyle accretion: `MDOT_BONDI(t)`, `MDOT_HOLE(t/M)`.
- Gas dynamics grid: `RHO_GAS(r,phi)`, `VRAD(r)`, `CS_PROFILE(r)`.
- Bolometric light curve: `L_BOL(t)`, peak time/flux.
- Spectral energy distribution mockup (SED).

## Phase 4 — Numerical Convergence
- Adaptive timestep convergence: `DT_MIN` history, `CFL_NUMBER`.
- Spatial resolution sweep: `N_PARTICLES = 1e3 → 1e5`.
- Integrator comparison: RK4 vs GEOKON vs Bulirsch-Stoer.
- Particle noise analysis: shot noise vs dynamical friction.

## Phase 5 — Publication Exports
- Full FITS suite with WCS headers (astropy-compliant).
- HDF5 particle dump for community codes (GADGET / AREPO).
- LaTeX table generator for paper (parameters, convergence).
- JWST / ELT mock images of stream (surface brightness).

## Simulation Bounds (Sgr A* baseline)
- `M_BH = 4e6 Msun` → `rg = 1.2e10 cm`
- `r_Bondi ≈ 1e4 rg` at `c_s = 10 km/s`
- `r_SOI   ≈ 1e6 rg` (~3 pc)
- `t_final = 10^5 M` (~1 day)

## Open targets
- `ENERGY_DRIFT < 1e-10`
- `J_CONSERVATION < 1e-12`

> Status: deferred to **beta 3+**. Not part of the 2D beta-2 release scope.
