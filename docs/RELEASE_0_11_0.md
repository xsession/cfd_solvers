# Release 0.11.0 - persistent campaign execution and external-solver workflow boundary

This release advances the SU2/csauto-inspired workflow layer from in-memory DOE primitives to a persistent, testable campaign execution substrate. It remains clean-room: no SU2 or csauto source is copied or translated.

## Added

- `CampaignTemplateFile`, `CampaignWriteOptions` and `CampaignPersistenceResult`.
- `write_campaign_case_folders(...)`, which creates one folder per campaign row, renders templates, writes `doe_row.csv`, records the template manifest and initializes a `registry.tsv`.
- `save_campaign_registry(...)` and `load_campaign_registry(...)` with deterministic TSV escaping for messages and status roundtrip.
- `SolverRuntime`, `SolverRunRequest` and `SolverCommandPlan`.
- `build_solver_command_plan(...)` for native, Docker, Apptainer/Singularity and Slurm-style wrappers.
- `detect_solver_outcome_from_log(...)`, `parse_residual_history(...)` and `parse_performance_metrics(...)`.
- `CampaignOptimizationConfig`, `CampaignOptimizationResult` and `run_gradient_descent_campaign(...)`.
- CLI smoke case `particle-campaign-execution`.
- Focused regression target `cfd-v0110-campaign-execution-tests`.

## Validation

- Full debug/CI-speed CPU CTest matrix: 97/97 passed.
- Focused ASan + UBSan + leak-detection run: `tests/test_v0110_campaign_execution.cpp` + `src/workflow/campaign.cpp` passed.

## Tracker impact

- Overall: 549/690 = 79.6%.
- Phase 9 workflow/UX: 19/21 = 90.5%.
- Phase 12 CST-class expansion remains 76/95 = 80.0%.

## Known limitations

- The command planner does not launch processes yet.
- Slurm/native/container commands are deterministic plans, not runtime-submitted jobs.
- Residual/performance parsers intentionally accept a compact common log grammar; solver-specific adapters can extend parsing later.
- The optimizer is finite-difference gradient descent, not a full SU2-style continuous/discrete adjoint implementation.
