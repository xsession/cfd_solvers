# Release 0.10.7 - distributed staggered PIC timestep bridge

## Added

- `run_serialized_distributed_staggered_pic_step_3d`, the first rank-local distributed timestep bridge for staggered 3-D EM-PIC.
- Domain-local particle push from local Ex/Ey/Ez/Bx/By/Bz blocks using the existing Vay pusher.
- Serialized particle ownership migration after the local push.
- Six-component EM field guard reconstruction for the next local update.
- Diagnostics for particle counts, max displacement, serialized exchange and guarded fields.
- CLI smoke case: `particle-distributed-step3d`.
- Regression target: `cfd-v0107-distributed-step-tests`.

## Validation

- Full debug/CI-speed CPU CTest matrix: 90/90 passed.
- Focused ASan + UBSan + leak detection for the v0.10.7 distributed timestep bridge: passed.

## Limitations

- Still rank-local orchestration, not a real MPI multi-rank run.
- Field guard exchange is reconstructed in-process; serialized EM field-vector envelopes remain future work.
- The step bridge validates orchestration and migration ordering, not production Esirkepov current deposition across subdomain boundaries.
