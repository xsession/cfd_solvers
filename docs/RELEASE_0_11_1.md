# Release 0.11.1 - local campaign runner and runtime doctor checks

This release continues the SU2/csauto-inspired workflow track by moving from deterministic command plans to a small local execution substrate.

## Added

- `doctor_solver_runtime(...)` for validating adapter identity, case-directory availability and launch command availability.
- `discover_campaign_outputs(...)` for residual/history and performance/timing file discovery inside generated case folders.
- `run_local_solver_case(...)` for POSIX native solver launches through argv-based `fork`/`execvp`, with captured stdout/stderr logs.
- `run_local_campaign(...)` for executing pending cases and updating the persistent `registry.tsv` from parsed solver outputs.
- New CLI smoke case `particle-campaign-local-runner`.
- New focused regression target `cfd-v0111-local-campaign-runner-tests`.

## Validation

- Full CPU CTest matrix: 99/99 passed.
- Focused local-campaign runner regression: pass.
- CLI smoke case `particle-campaign-local-runner`: pass.
- Focused ASan + UBSan + leak detection for the v0.11.1 local-campaign runner slice: passed.

## Tracker movement

- Overall tracker: 552/693 = 79.7%.
- Phase 9 workflow/UX: 22/24 = 91.7%.

## Limitations

- The native runner is implemented for POSIX systems. Windows support needs a separate `CreateProcess` implementation.
- Docker, Apptainer/Singularity and Slurm are still command-planned and doctor-checked, but not executed by the local runner.
- There is still no web dashboard, live steering API or scheduler polling loop.
