# cfd_solvers 0.9.4 - 2-D electrostatic PIC foundation

## Summary

This release adds the first multi-dimensional PIC foundation layer after the v0.9.3 self-consistent 1-D/3-V EM-PIC baseline. The new code focuses on analytically testable 2-D building blocks before extending the full Yee electromagnetic update to 2-D/3-D.

## Added

- `PicParticle2D`, `ElectrostaticPic2DConfig` and `ElectrostaticPic2D`.
- Periodic 2-D CIC charge deposition.
- Periodic 2-D spectral Poisson electric-field solve.
- Bilinear 2-D electric-field gather and particle stepping.
- 2-D spectral charge-conserving current reconstruction satisfying `drho/dt + div(J) = 0` for old/new CIC charge density.
- 2-D particle wall handling with absorption, specular reflection and secondary-emission macro-weight reporting.
- 2-D electric-wall and absorbing-sponge field-boundary primitives.
- CLI smoke case `particle-pic2d`.
- Focused regression target `cfd-v094-pic2d-tests`.

## Validation

The focused v0.9.4 target validates:

- manufactured sinusoidal 2-D periodic Poisson fields;
- 2-D charge-continuity residual from spectral current reconstruction;
- finite 2-D field gather and explicit PIC stepping;
- particle wall reflection/absorption plus secondary-emission reporting;
- electric-wall tangential-field clamping and absorbing-sponge damping.

## Tracker

- Overall tracker: 500/642 capabilities, 77.9%.
- Phase 12 CST-class expansion: 34/54 capabilities, 63.0%.

## Still open

The next step is a full 2-D/3-V electromagnetic PIC update: Yee `Ez/Bx/By` and/or full 2-D transverse-field storage, charge-conserving current deposition at Yee locations, metallic/absorbing field boundaries coupled directly to particles, and then domain decomposition/GPU kernels.
