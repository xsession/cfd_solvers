# Release 0.10.0 - local finite-volume 3-D EM-PIC current deposition

## Added

- Finite-volume local 3-D current reconstruction from old/new trilinear CIC charge density.
- `CurrentDeposition3DMode` selector for `StaggeredElectromagneticPic3D`.
- `particle-local-current3d` CLI smoke case.
- `cfd-v0100-local-current-tests` focused regression target.
- Configurable polynomial order for the 3-D absorbing-sponge boundary profile.

## Validation

- Full debug/CI-speed CPU CTest matrix: 72/72 passed.
- Focused ASan + UBSan + leak detection for the v0.10.0 local-current slice: passed.

## Notes

This release does not claim production-grade Esirkepov/Villasenor-Buneman deposition. The new local current path is a validated finite-volume bridge that satisfies periodic discrete continuity without relying on spectral Fourier reconstruction.
