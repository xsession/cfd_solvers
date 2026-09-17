# Release 0.9.9 - staggered 3-D Yee EM-PIC foundation

This release continues the CST-class particle/PIC roadmap by moving the 3-D EM-PIC branch from a compact nodal periodic reference to a compact staggered Yee reference.

## Added

- `StaggeredElectromagneticPic3D` with component-specific Ex/Ey/Ez and Bx/By/Bz Yee locations.
- 3-D Yee CFL guard.
- Staggered Faraday and Ampere curl updates.
- Component-offset trilinear E/B gather for 3-D particles.
- Spectral charge-continuity Jx/Jy/Jz coupling reused from the validated 3-D current reconstruction.
- 3-D electric-wall and absorbing-sponge field-boundary primitives.
- 3-D absorbing/specular particle-wall handling with secondary-emission macro-weight reporting.
- CLI smoke case `particle-staggered-em-pic3d`.
- Regression target `cfd-v099-staggered-pic3d-tests`.

## Validation

- Full debug/CI-speed CPU CTest matrix: 70/70 passed.
- Focused ASan + UBSan + leak detection: passed for `tests/test_v099_staggered_pic3d.cpp` plus `src/particle/electromagnetic.cpp`.

## Tracker

- Overall: 515/657 capabilities, 78.4%.
- Phase 12: 49/69 capabilities, 71.0%.

## Known limitations

This is a compact CPU reference, not yet a production EM-PIC engine. It still uses spectral continuity-current reconstruction rather than local Esirkepov/Villasenor-Buneman deposition, and it does not yet include CPML, conductor-material boundaries, particle sorting, MPI guard-cell exchange or GPU kernels.
