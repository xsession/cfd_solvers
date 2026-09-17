# cfd_solvers 0.11.4

## Focus

This release upgrades the multi-server campaign layer from static placement/deploy planning to a supervised execution contract that can be reviewed, tested and later connected to real SSH/Docker hosts.

## Added

- `RemoteSupervisedJobPlan` for per-case remote launch, status, cancel and log-fetch command plans.
- `RemoteHostHealthCheck` and `plan_multiserver_execution(...)` for remote-root, runtime and disk-space checks.
- `write_multiserver_execution_files(...)` for auditable shell artifacts and `multiserver_jobs.tsv`.
- `parse_remote_job_status_table(...)` and `apply_remote_job_status_to_registry(...)`.
- Docker Compose healthchecks plus `deploy/healthcheck.sh` and `deploy/worker-entrypoint.sh`.
- CLI smoke case `particle-multiserver-execution`.
- Focused regression target `cfd-v0114-multiserver-execution-tests`.

## Validation

- Focused multi-server execution regression: passed.
- Full CTest matrix: see packaged status JSON.

## Limitations

The new layer still generates deterministic command scripts and status parsers. It does not yet manage SSH keys/secrets, run persistent daemons, stream logs over a live API, retry failed remote jobs or supervise Docker daemons in production.
