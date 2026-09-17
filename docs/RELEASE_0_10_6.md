# Release 0.10.6 - distributed PIC exchange round and EM field guards

## Scope

This release composes the previously separate PIC execution layers into a single validated rank-local exchange round and adds a six-component electromagnetic field guard wrapper for decomposed 3-D EM-PIC fields.

## Added

- `ElectromagneticFieldBlocks3D` and `ElectromagneticGuardedFieldBlocks3D`.
- `exchange_electromagnetic_field_guard_cells_3d(...)` for Ex/Ey/Ez/Bx/By/Bz guard exchange.
- `DistributedPicExchangeRound3D`.
- `run_serialized_distributed_pic_exchange_round_3d(...)`.
- `particle-distributed-round3d` CLI smoke case.
- `particle-field-guards3d` CLI smoke case.
- `cfd-v0106-distributed-round-tests` focused regression target.

## Validation

- Full debug/CI-speed CPU CTest matrix: 88/88 passed.
- Focused ASan + UBSan + leak detection for the v0.10.6 distributed-round slice: passed.

## Tracker

- Overall: 534/675 capabilities, 79.1%.
- Phase 12 CST-class expansion: 68/87 capabilities, 78.2%.

## Limitations

This is still a deterministic rank-local orchestration layer, not a real MPI multi-rank run and not GPU execution. The next implementation step should run these serialized envelopes through an MPI runtime, then connect the exchange round into rank-local staggered 3-D EM-PIC stepping.
