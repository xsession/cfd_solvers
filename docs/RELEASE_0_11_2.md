# cfd_solvers v0.11.2

This checkpoint adds the next SU2/csauto-inspired campaign workflow layer: live-control directive files, scheduler-status parsing, and output-based registry refresh.

## Added

- `include/cfd/workflow/campaign_control.hpp` and `src/workflow/campaign_control.cpp`.
- Adapter-gated control requests for `stop`, `extend`, `checkpoint` and `flush`.
- Per-case `.cfd_control/*.directive` persistence and discovery.
- Slurm-style queue/status parser with normalized queued/running/done/failed/cancelled states.
- Registry status application from scheduler records.
- Registry refresh from discovered case logs, residual histories and performance files.
- CLI smoke case `particle-campaign-control`.
- Focused test target `cfd-v0112-campaign-control-tests`.

## Validation

- Full CPU CTest matrix: 101/101 passed.
- Focused ASan + UBSan + leak detection for the v0.11.2 campaign-control slice: passed.

## Current limitation

The new layer writes/reads durable live-control directives and parses scheduler status, but it does not yet submit Slurm jobs, poll a live scheduler command, or run Docker/Apptainer containers. Those runtimes can reuse the command-plan, doctor, directive and registry-refresh contracts added in v0.11.0-v0.11.2.
