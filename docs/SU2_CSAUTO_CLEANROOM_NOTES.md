# SU2 and csauto clean-room notes for v0.10.9

This note records what was learned from public repository inspection and what was implemented independently in `cfd_solvers`.

## Sources inspected

- `su2code/SU2` default branch `master`.
- `simvia-tech/csauto` default branch `main`.

No source code was copied or mechanically translated. SU2 is LGPL-2.1 licensed and csauto is GPL-3.0 licensed, so both are treated as clean-room architecture and workflow references for this MIT repository.

## SU2 lessons applied

SU2 is organized as a C++ PDE/CFD and optimization code with separate solver, iteration, numerics, output and Python workflow layers. The public tree shows dedicated C++ directories such as `drivers`, `iteration`, `numerics`, `solvers`, `variables` and Python-side workflows including continuous/discrete adjoint, direct differentiation, finite differences and uncertainty scripts.

Transferable clean-room ideas implemented here:

- keep primal solver kernels separate from campaign/optimization orchestration;
- represent a solver's externally visible capabilities with a small descriptor rather than hard-wiring workflow code to solver names;
- add simple finite-difference gradient and gradient-descent helpers as validation-facing building blocks before attempting full adjoints;
- keep gradient utilities regression-testable on tiny analytic objectives.

Implemented files:

- `include/cfd/workflow/campaign.hpp`
- `src/workflow/campaign.cpp`
- `tests/test_v0109_transport_campaign.cpp`

## csauto lessons applied

csauto demonstrates a campaign layer around CFD solvers: one DOE row generates one case, template files use placeholders and conditional blocks, runtime-specific logic is separated from solver-specific logic, and solver capabilities are provided through an adapter boundary.

Transferable clean-room ideas implemented here:

- factorial and Latin-hypercube DOE generation with deterministic seeds;
- `{placeholder}` rendering plus `<!-- IF key -->`, `<!-- IF key=value -->` and `<!-- IF key!=value -->` conditional blocks;
- `SolverAdapterDescriptor` and `SolverCapability` as a small boundary for future external solver adapters;
- in-memory campaign registry updates and summaries to support future dashboard/export work.

## Current limitations

- No web dashboard yet.
- No production Slurm/Docker/Apptainer execution yet.
- No live steering API yet.
- No web dashboard yet.
- No true continuous/discrete adjoint implementation yet.

These limitations are deliberate: the current release validates the reusable in-process primitives first.

## v0.11.0 follow-up - persistent campaign execution boundary

The next clean-room step applies the csauto idea that a campaign should have a stable on-disk shape before runtime orchestration is added. The implementation writes one directory per generated case, renders all template files into that directory, stores a deterministic `doe_row.csv`, records a template manifest, and initializes a small `registry.tsv` that can be roundtripped by tests.

The solver-adapter lesson is still kept at the boundary. `SolverAdapterDescriptor` declares identity, native/container binary names and capabilities; `SolverRunRequest` adds runtime, rank/thread counts and extra arguments; `build_solver_command_plan` returns argv/environment/display strings for native, Docker, Apptainer/Singularity and Slurm-style execution. It does not execute the solver yet, which keeps the core deterministic and testable without SU2, code_saturne, Slurm or containers installed.

The SU2 optimization lesson is represented by parser and optimization seams rather than a copied adjoint stack: residual histories and performance metrics can be extracted from common log shapes, objective values can update the campaign registry, and `run_gradient_descent_campaign` provides a small finite-difference design-loop baseline. Future solver-specific adapters should replace only the command/parser pieces, not the generic campaign folder, registry or optimizer contracts.


## v0.11.1 follow-up - safe local process execution

The next clean-room step adds a native local runner behind the existing solver-adapter and command-plan boundary. The implementation validates the runtime with `doctor_solver_runtime`, launches native processes through argv-based `fork`/`execvp` on POSIX systems, captures `solver.stdout.log` and `solver.stderr.log`, discovers residual/history and performance/timing files, and updates the persistent registry after each case.

This keeps csauto's useful separation of generic campaign orchestration from solver-specific conventions while avoiding any copied implementation. Container and scheduler execution remain planned/validated interfaces; the first real runner is native-only so the behavior can be tested without Docker, Apptainer, Slurm or external CFD software.


## v0.11.2 follow-up - live control and scheduler monitoring contracts

The next clean-room step extends the campaign layer toward daily operations without adding unvalidated external runtime dependencies. `write_campaign_control_directive` creates solver-adapter-approved stop/extend/checkpoint/flush requests as deterministic files under each case's `.cfd_control` directory, and `discover_campaign_control_directives` lets a runner or future dashboard list pending steering requests.

The scheduler side remains transport-neutral: `parse_slurm_queue_table` accepts compact Slurm-like status tables and maps queued/running/done/failed/cancelled states into the campaign registry through `apply_scheduler_records_to_registry`. `refresh_campaign_status_from_outputs` separately scans case logs and residual/performance files to update objective, iteration and status fields after an external runtime has produced output.

This follows csauto's useful runtime-monitoring idea while preserving the existing adapter boundary and keeping the implementation independent of csauto source code. Docker, Apptainer and real Slurm submission are still command-plan/doctor-level boundaries until a runtime-specific executor is validated.

## v0.11.3 deployment workflow extension

The csauto lesson extended in this release is that campaign orchestration should remain behind explicit boundaries: solver adapters describe how a solver runs, while deployment planning describes where a case should run. The new multi-server planner therefore emits auditable assignments and command plans instead of embedding remote execution into the solver adapter. The Docker scaffold follows the same split: Compose describes worker/manager services, while generated campaign folders and registry files remain external persistent data.


## v0.11.4 campaign supervision follow-through

The csauto-derived lesson of keeping solver/runtime details behind explicit boundaries is extended here from planning into supervised execution artifacts. The new layer does not hide remote execution behind implicit side effects: it writes health, launch, status and log-fetch scripts plus a job metadata table that can be inspected, tested, versioned and later replaced by a real remote executor.

Transferred clean-room ideas:

- keep scheduler/runtime status parsing separate from solver-specific outcome parsing;
- make dashboard/API-ready state transitions parseable from a small status table;
- keep live-control and log-fetch conventions outside solver numerics;
- keep Docker deployment healthchecks declarative and reproducible.


## v0.11.5 deployment-supervision follow-through

The csauto-inspired campaign layer now has another clean-room step beyond static Docker/remote command generation: a deterministic supervision contract. The design keeps generic workflow code separated from solver-specific launch details through the existing solver-adapter descriptor, then derives access checks, retry wrappers, log-tail commands and dashboard summaries from the multi-server job metadata.

The key design rule remains: generated files are inspectable artifacts, not a hidden remote daemon. Operators can review `multiserver_access_checks.sh`, `multiserver_retry_launch.sh`, `multiserver_tail_logs.sh`, `multiserver_dashboard.json` and `multiserver_supervision.tsv` before wiring them to SSH, Docker or a future web service.

The new redaction helper is intentionally conservative. It hides obvious inline values after markers such as `TOKEN=`, `PASSWORD=`, `SECRET=` and `KEY=`, but production deployments should still use a real secret manager or controlled environment files.

## v0.11.6 controller/API follow-through

The csauto-inspired workflow now has a generated controller scaffold over the existing multi-server supervision artifacts. This follows the same clean-room lesson as the solver-adapter boundary: the controller consumes generic status/log/control metadata rather than embedding solver-specific file conventions.

Transferred clean-room ideas:

- expose status and control through a stable HTTP/API seam after the registry and directive-file contracts are deterministic;
- keep mutating control actions capability-gated and token-gated;
- publish route metadata through OpenAPI and a TSV route table so a future UI can be generated or validated from the same source;
- keep the first controller dependency-free and operator-reviewable before adding a richer dashboard or daemon.

The generated controller is deliberately small. It serves status snapshots and writes `.cfd_control` directives; it does not replace the remote execution scripts, scheduler integration or future log-stream service.

## v0.11.7 dashboard follow-through

The v0.11.7 work continues the csauto-inspired operator-surface direction without copying its frontend. The clean-room transfer is architectural: a generated controller can expose status, logs, control actions and a lightweight dashboard while solver-specific behavior remains behind adapter/campaign abstractions.

The generated dashboard is intentionally simpler than csauto's web UI. It uses static HTML/CSS/JavaScript and the existing dependency-free controller routes so that the release tests can validate the contract offline. Future production work should replace polling with a supervised event/log stream and add real authentication/TLS at the deployment boundary.

## v0.11.8 live-event follow-through

The controller now adds a clean-room SSE status seam rather than adopting a framework-specific frontend or remote-execution implementation. Events are derived from the same generic campaign status contract, keeping solver adapters, placement logic and the operator surface separated. Polling remains as a fallback, and production TLS, identity and process supervision remain deployment concerns.

## v0.12.0 production-boundary follow-through

The production layer remains independent of solver internals. TLS, token roles, audit persistence, probe execution and daemon artifacts wrap the generic campaign-controller contract rather than entering SU2-style solver adapters or csauto-inspired campaign planning. Secrets remain operator-owned external inputs, and live probes reuse the generated supervision boundary.
