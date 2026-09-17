# cfd_solvers v0.10.9 release notes

## Theme

v0.10.9 adds a validated clean-room bridge from Geant4-style particle-through-matter transport into dose/SAR/bioheat scoring, then adds SU2/csauto-inspired campaign workflow primitives for DOE, template rendering, solver-adapter descriptors, registry summaries and finite-difference gradient studies.

## Implemented

- `TransportPhysicsList` and deterministic stochastic interaction-length sampling.
- Axis-aligned `TransportRegionBvh` region lookup.
- Sensitive detector declarations, hit collection and per-step detector scoring.
- 3-D dose-grid accumulation from deposited energy.
- Dose-rate to SAR conversion and 2-D SAR projection into the Pennes bioheat solver.
- Factorial and Latin-hypercube campaign generation.
- Template placeholder and IF/ENDIF conditional rendering.
- Solver-adapter descriptor/capability validation.
- Campaign registry updates/summaries.
- Finite-difference gradient and gradient-descent utility functions.
- CLI smoke cases `particle-transport-dose-bvh` and `particle-campaign-doe`.
- Focused regression target `cfd-v0109-transport-campaign-tests`.

## Validation

- Full debug/CI-speed CPU CTest matrix: 95/95 passed.
- Focused ASan + UBSan + leak detection for the v0.10.9 transport/campaign slice: passed.

## Clean-room basis

- Geant4 was used for track/step/process/scoring architecture only.
- SU2 was used for solver/iteration/workflow/optimization organization only.
- csauto was used for campaign automation and solver-adapter boundary ideas only.
- No Geant4, SU2 or csauto source code was copied or translated.

## Known limitations

- No production particle physics data tables.
- No hierarchical constructive geometry or repeated placements.
- No external solver launcher or dashboard.
- No real adjoint solver yet; only finite-difference utilities are present.
- No quantitative dose validation beyond regression-scale conservation/sanity checks.
