# cfd_solvers v0.15.4 - characteristic WENO compressible transport

v0.15.4 closes the remaining Phase-3 high-order compressible reconstruction item while preserving the existing MUSCL-minmod/Rusanov path.

## Added

- `CompressibleReconstruction::characteristic_weno5` for the 1-D ideal-gas Euler solver.
- Roe-averaged conservative-variable eigenbasis used to project five-cell stencils into characteristic coordinates.
- Jiang-Shu fifth-order WENO candidate polynomials, smoothness indicators and nonlinear weights.
- Local first-order face-state fallback when a reconstructed density/pressure state is nonphysical.
- `CompressibleTimeIntegrator::ssprk3` method-of-lines integration while retaining forward Euler.
- Arbitrary cell-profile initialization for smooth and shock/entropy regression cases.
- Conserved total momentum diagnostic.
- Normalized cell pressure-jump shock sensor intended for later AMR coupling.
- Focused smooth-wave, uniform-state, conservation, positivity and shock-sensor regression coverage.

## Validation

The periodic entropy-wave regression initializes exact finite-volume cell averages. At 40/80/160 cells the L1 density errors are `3.93392e-6`, `1.25593e-7` and `4.26201e-9`, corresponding to observed rates about 4.97 and 4.88 over this resolution range. The test also compares against the existing MUSCL path; the acceptance threshold is deliberately lower than five because SSPRK3 is third-order in time and asymptotic behavior depends on CFL and problem scaling.

The complete default CPU matrix passes **124/124 CTest targets**.

The Sod regression exercises characteristic reconstruction at a discontinuity and requires positive density and pressure throughout the evolved solution.

## Tracker

- Characteristic high-order/WENO compressible schemes: complete.
- Phase 3: 98/106 capabilities complete.
- Overall: 612/719 capabilities complete.

See `UPSTREAM_DOCUMENTATION_REVIEW_0_15_4.md` for the clean-room numerical-method mapping and explicit scope limits.
