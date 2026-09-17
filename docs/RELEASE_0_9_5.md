# Release 0.9.5 - periodic 2-D/3-V electromagnetic PIC baseline

This release extends the CST-class particle/plasma track from the v0.9.4
2-D electrostatic PIC foundation to a first self-consistent 2-D/3-velocity
electromagnetic PIC baseline.

## Added

- `ElectromagneticPic2D`, a periodic 2-D/3-V electromagnetic PIC solver with:
  - collocated periodic `Ex/Ey/Ez` and `Bx/By/Bz` arrays;
  - centered curl updates for the electromagnetic fields;
  - bilinear field gather;
  - Vay relativistic particle push in 3 velocity components;
  - charge-conserving spectral in-plane `Jx/Jy` current reconstruction from old/new CIC charge density;
  - CIC transverse `Jz` deposition that couples into the `Ez` update;
  - optional periodic Poisson correction for the longitudinal in-plane electric field;
  - Monte-Carlo neutral-collision coupling that maps secondary particles back into the 2-D PIC population.
- `particle-em-pic2d` CLI/CTest smoke case.
- Focused `cfd-v095-em-pic2d-tests` regression target covering:
  - bounded short-run vacuum TM-wave energy;
  - in-plane current continuity residuals;
  - transverse `Jz -> Ez` coupling;
  - field gather and MCC coupling.

## Validation

- Focused v0.9.5 regression target: passes.
- `particle-em-pic2d` CLI smoke case: passes.
- Full normal CPU CTest matrix: 62/62 passes in this environment.
- Focused ASan+UBSan+leak-detection executable for the new v0.9.5 slice: passes.

## Tracker

- Phase 12 advances from 34/54 to 36/56 validated capabilities.
- Overall tracker advances from 500/642 to 502/644 validated capabilities.

## Notes and limitations

This is a deliberately compact baseline. It is periodic and uses centered curl
operators so it can validate particle/current/field coupling without mixing in
open-boundary or conductor-wall complexity. It is not yet a full production
2-D Yee/PML/decomposition implementation. The next logical step is a true
staggered 2-D Yee update with absorbing/metallic boundaries and then 3-D
deposition/field gather plus distributed/GPU execution.
