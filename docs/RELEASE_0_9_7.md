# Release 0.9.7 - 3-D electrostatic PIC foundation

## Summary

v0.9.7 extends the PIC/plasma roadmap from 2-D Yee EM-PIC toward full 3-D by adding a validated 3-D electrostatic particle/field foundation.

The release intentionally stops before full 3-D electromagnetic Yee fields. The goal is to validate trilinear 3-D deposition, periodic 3-D Poisson field reconstruction, trilinear field gather, particle stepping and spectral current reconstruction first.

## Added

- `PicParticle3D` and `ElectrostaticPic3D` in `include/cfd/particle/electromagnetic.hpp`.
- Trilinear 3-D CIC charge deposition on periodic structured grids.
- Periodic 3-D spectral Poisson electric-field reconstruction.
- Trilinear 3-D electric-field gather.
- Spectral 3-D charge-conserving current reconstruction from old/new CIC charge density.
- CLI smoke case: `particle-pic3d`.
- Focused regression target: `cfd-v097-pic3d-tests`.

## Validation

- Full debug/CI-speed CTest matrix: 66/66 passed.
- Focused ASan+UBSan+leak-detection executable for the new 3-D PIC slice: passed.
- New focused checks cover:
  - 3-D CIC total charge conservation;
  - manufactured sinusoidal periodic 3-D Poisson Ex/Ey/Ez;
  - 3-D spectral current continuity residual;
  - electrostatic 3-D PIC field gather and time stepping.

## Tracker

- Overall: 508/650 capabilities complete (78.2%).
- Phase 12: 42/62 capabilities complete (67.7%).

## Still open

- Full 3-D Yee EM-PIC field staggering.
- Local charge-conserving current deposition on a 3-D staggered grid.
- Metallic/PML boundaries for 3-D EM-PIC.
- Domain decomposition, GPU kernels and load balancing for particle/field updates.
- Plasma chemistry and breakdown workflows.
