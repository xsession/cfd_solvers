# Release 0.9.6 - staggered Yee EM-PIC and self-study docs

## Summary

v0.9.6 adds the first true staggered 2-D Yee electromagnetic PIC baseline and a self-study documentation layer for the implemented solver families.

The numerical change builds on the v0.9.5 periodic centered EM-PIC baseline by adding component-specific field locations, a Yee-style staggered curl update, an explicit 2-D CFL guard, electric-wall and absorbing-sponge field boundaries, and integrated non-periodic particle wall handling.

The documentation change adds a learning map for the implemented code, with source files, tests, CLI cases and curated external resources.

## Added

- `StaggeredElectromagneticPic2D` in `include/cfd/particle/electromagnetic.hpp` and `src/particle/electromagnetic.cpp`.
- Component-specific Ex/Ey/Ez/Bx/By/Bz sampling locations exposed through compact same-sized arrays.
- Staggered forward/backward curl operators for 2-D/3-V EM-PIC.
- CFL guard for the 2-D Yee update.
- Integrated electric-wall field boundary enforcing tangential electric-field zeros.
- Integrated absorbing-sponge damping on all field components.
- Non-periodic absorbing/specular particle-wall handling inside the staggered PIC time step.
- Secondary-emission macro-weight reporting for absorbed particle impacts.
- Reuse of Vay particle push, spectral charge-conserving in-plane current reconstruction, transverse `Jz` deposition and MCC collision coupling.
- CLI smoke case: `particle-staggered-em-pic2d`.
- Focused regression target: `cfd-v096-staggered-pic2d-tests`.
- `docs/IMPLEMENTED_SOLVER_STUDY_GUIDE.md`.
- `docs/LEARNING_RESOURCES.md`.

## Validation

- Full debug/CI-speed CTest matrix: 64/64 passed.
- Focused ASan+UBSan+leak-detection executable for the new staggered EM-PIC slice: passed.
- New focused checks cover:
  - bounded source-free vacuum TM energy;
  - electric-wall tangential-E enforcement;
  - absorbing-sponge field damping;
  - particle absorption and specular reflection;
  - transverse current deposition coupling into `Ez`;
  - MCC collision invocation.

## Tracker

- Overall: 505/647 capabilities complete (78.1%).
- Phase 12: 39/59 capabilities complete (66.1%).

## Still open

- Local Esirkepov-style current deposition on the true staggered Yee grid.
- Full 3-D Yee EM-PIC deposition/gather.
- PML-grade absorbing boundaries for EM-PIC.
- Domain decomposition and GPU kernels for particle/field updates.
- Plasma chemistry and breakdown workflows.
