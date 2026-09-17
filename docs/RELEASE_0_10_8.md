# Release 0.10.8 - Geant4-inspired track/process transport baseline

This release adds the first clean-room particle-through-matter transport layer
inspired by Geant4's public architecture.

## Added

- `TransportTrack`, `TransportStep`, `TransportProcess`, `TransportRegion`,
  `TransportMaterial` and `TransportWorld`.
- Geometry/material location over axis-aligned slab/box regions.
- Step arbitration between user maximum step, geometry boundary and
  process-limited physical interaction length.
- Continuous stopping-power energy loss process.
- Discrete interaction process with parent energy loss and secondary production.
- Production-cut and energy-cut handling.
- Run-level scoring for deposited energy, track length, steps, boundary crossings
  and secondaries.
- CLI smoke case: `particle-geant4-transport`.
- Focused regression target: `cfd-v0108-geant4-transport-tests`.
- Geant4 clean-room learning notes in `docs/GEANT4_CLEANROOM_NOTES.md`.

## Validation

- Continuous loss, GPIL/discrete-process step selection, secondary production,
  geometry boundary limiting, run scoring and invalid-input guards are covered by
  `cfd-v0108-geant4-transport-tests`.

## Limitations

- This is not a Geant4 adapter and does not read Geant4 geometry or physics
  lists.
- Region navigation is an axis-aligned box/slab baseline, not a full CSG/BVH
  navigator.
- Process lengths are deterministic coefficients; stochastic sampled interaction
  lengths and real cross-section tables are future work.

## Release validation in this workspace

- Full debug/CI-speed CPU CTest matrix: 92/92 passed.
- Focused ASan+UBSan+leak-detection executable for the v0.10.8 transport slice: passed.
- Integration tracker after this release: 538/679 = 79.2% overall; Phase 12: 72/91 = 79.1%.
