# Release 0.9.8 - periodic 3-D/3-V electromagnetic PIC baseline

## Summary

v0.9.8 converts the validated 3-D electrostatic PIC geometry from v0.9.7 into a first self-consistent 3-D electromagnetic PIC reference. The new `ElectromagneticPic3D` class is intentionally compact and periodic: it stores all electric and magnetic components on a nodal 3-D grid, advances centered Maxwell curls, gathers trilinear E/B fields to macro-particles, pushes particles with the Vay pusher and drives field sources through the existing spectral 3-D charge-continuity current reconstruction.

This is a reference baseline for validating particle-field plumbing before the production staggered-Yee/PML/GPU path.

## Added

- `cfd::particle::ElectromagneticPic3DConfig`.
- `cfd::particle::ElectromagneticPic3DDiagnostics`.
- `cfd::particle::ElectromagneticPic3D`.
- Full periodic 3-D Ex/Ey/Ez and Bx/By/Bz field storage.
- Centered periodic curl updates for Faraday and Ampere equations.
- Compact 3-D CFL guard for the reference time update.
- Trilinear 3-D E/B field gather.
- Vay-pusher coupling for `PicParticle3D` macro-particles.
- Spectral charge-conserving Jx/Jy/Jz source coupling using the old/new charge-density reconstruction from v0.9.7.
- Optional periodic Poisson correction for longitudinal electric fields.
- Monte-Carlo neutral collision coupling inside the 3-D EM-PIC step.
- CLI smoke case: `particle-em-pic3d`.
- Focused regression target: `cfd-v098-em-pic3d-tests`.

## Validation

- Full debug/CI-speed CPU CTest matrix: 68/68 passed.
- Focused ASan + UBSan + leak detection for the new v0.9.8 3-D EM-PIC slice: passed.
- Focused tests cover bounded vacuum-wave energy, trilinear E/B gather, charge-continuity residual, Poisson correction, MCC coupling and CFL rejection.

## Tracker impact

- Overall tracker: 511/653 capabilities, 78.3%.
- Phase 12 CST-class electromagnetic/particle/bioelectromagnetic coverage: 45/65, 69.2%.

## Known limitations

- The new solver is a compact periodic nodal-grid reference, not yet a production staggered 3-D Yee solver.
- It does not yet include local Esirkepov/Villasenor-Buneman current deposition, conductor boundaries, PML, MPI/GPU execution or plasma chemistry.
- Existing unrelated legacy translation units still emit warnings in full builds; the new v0.9.8 focused slice passes sanitizer validation.

## Next recommended step

Implement a true staggered 3-D Yee EM-PIC layer with component-specific field locations, local charge-conserving current deposition, PEC/PML boundaries, then parallel/GPU particle-field execution.
