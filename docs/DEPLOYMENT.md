# Deployment and multi-server campaign planning

The v0.11.3 workflow layer adds deterministic deployment planning for generated campaign folders. The goal is to keep the execution contract inspectable before adding privileged remote execution.

## Multi-server placement

`CampaignServerDescriptor` describes a worker host:

- `name`: stable identifier used in assignment files;
- `host` and `user`: SSH target metadata;
- `remote_root`: campaign root on the server;
- `runtime`: `native` or `docker`;
- `slots`: maximum concurrent case placements;
- `tags`: optional labels such as `cpu`, `gpu`, `linux`;
- `docker_image`: required for Docker workers;
- `online`: planner-side availability flag.

`plan_multi_server_campaign(...)` consumes a campaign registry and emits:

- case-to-server assignments;
- remote case paths;
- native or Docker solver command plans;
- mkdir, rsync and run command plans.

`write_multiserver_plan_files(...)` persists those plans as:

- `multiserver_assignments.tsv` for audit/review;
- `multiserver_commands.sh` for operator-controlled execution.

The generated script is intentionally simple. It is not a security boundary and should be reviewed before production use.

## Docker deployment scaffold

`write_docker_deploy_system(...)` writes a local Docker/Compose scaffold:

- `Dockerfile` builds a minimal runtime image layout;
- `compose.yaml` defines a campaign manager service and worker service;
- `deploy/deploy.sh` builds and starts the Compose project;
- `deploy/servers.example.tsv` documents the expected multi-server inventory shape;
- `deploy/env.example` captures common deployment variables.

Typical generated flow:

```bash
cd deploy_stack
./deploy/deploy.sh
```

For multi-server usage, copy or adapt `deploy/servers.example.tsv`, generate campaign case folders, then use the planner to produce `multiserver_assignments.tsv` and `multiserver_commands.sh`.

## Current limitations

The planner does not yet execute remote commands itself, track remote process IDs, stream remote logs, or manage credentials. Docker/Compose files are generated and structurally tested, but the release validation does not require a local Docker daemon.


## v0.11.4 supervised multi-server execution

After creating a multi-server placement plan, v0.11.4 can derive a supervised execution plan with:

- `multiserver_health.sh` for remote-root/runtime/disk-space checks;
- `multiserver_launch.sh` for per-case supervised `nohup` launch commands;
- `multiserver_status.sh` for parseable `case_id/server/state/exit_code` status rows;
- `multiserver_fetch_logs.sh` for collecting `.cfd_run` logs;
- `multiserver_jobs.tsv` with stdout/stderr/pid/exit-code paths and the generated command text.

The generated Docker/Compose scaffold now also includes container healthchecks, `deploy/healthcheck.sh`, and `deploy/worker-entrypoint.sh`. These are still intended as reviewable deployment scaffolds, not as a secrets manager or production-grade remote supervisor.


## v0.11.5 supervision hardening

The v0.11.5 layer keeps the remote execution model script-based, but adds the artifacts needed by a future dashboard or daemon:

- `multiserver_access_checks.sh` checks local/SSH access, Docker daemon availability for Docker workers, and campaign-root writability.
- `multiserver_retry_launch.sh` wraps each launch command with a bounded retry loop and configurable delay.
- `multiserver_tail_logs.sh` emits stdout/stderr tail commands for each remote job.
- `multiserver_dashboard.json` is a small dashboard-ready case/status/log-path summary.
- `multiserver_supervision.tsv` is a human-readable summary of case, server, state and log paths.

Command-display redaction helpers remove inline values following markers such as `TOKEN=`, `PASSWORD=`, `SECRET=` and `KEY=`. This is intended to avoid leaking obvious secrets into generated display/status files. Production deployments should still pass credentials through a real secret manager or a controlled environment file instead of inline command arguments.

The release tests validate the planning, redaction and file-generation contracts offline. They do not open SSH connections, contact a Docker daemon, or run remote jobs.

## v0.11.6 controller API scaffold

The v0.11.6 layer adds a generated controller directory on top of the v0.11.5 supervision files. The controller is intentionally small and dependency-free so it can be reviewed before it is exposed on a network.

Generated artifacts:

- `controller/cfd_controller.py` - Python standard-library HTTP controller using `ThreadingHTTPServer`.
- `multiserver_controller.sh` - root-level launcher that starts the controller from the generated controller directory.
- `controller/openapi.json` - machine-readable API route summary.
- `controller/status.json` - dashboard-ready case/status/log-path snapshot.
- `controller/routes.tsv` - compact route/security/mutation table.
- `controller/env.example` - environment variables for host, port, token and campaign root.
- `controller/README.md` - operator notes.

Routes:

- `GET /api/health` returns a health payload.
- `GET /api/status` returns the generated status JSON.
- `GET /api/cases` returns the case table.
- `GET /api/cases/{case_id}/logs/{stream}` tails stdout or stderr when the referenced path is readable from the controller host.
- `POST /api/cases/{case_id}/control/{action}` writes `.cfd_control/<action>.directive` and `latest.directive` when mutation is enabled and the token check passes.

Security notes:

- Set `CFD_CONTROLLER_TOKEN` before exposing mutating routes.
- Keep the controller bound to `127.0.0.1` unless it is behind an authenticated reverse proxy or private network.
- The controller writes directive files only; remote workers still need their own polling/control integration.
- This is not yet a production daemon, secrets manager, SSH orchestrator, WebSocket log streamer or live Docker probe service.

## v0.11.7 generated controller dashboard

The v0.11.7 layer adds a small browser UI to the generated controller directory. It is intentionally dependency-free: no Node.js, no Svelte/React build, and no Python packages beyond the standard library controller.

Generated files:

- `controller/index.html` - dashboard shell served at `/` and `/index.html`.
- `controller/dashboard.js` - polling dashboard client for `/api/health`, `/api/status`, log-tail routes and token-prompted control actions.
- `controller/dashboard.css` - compact dark-mode operator styling.
- `controller/events.ndjson` - newline-delimited event snapshot for downstream dashboards or log collectors.

Additional routes:

- `GET /` and `GET /index.html` serve the dashboard.
- `GET /dashboard.js` and `GET /dashboard.css` serve static assets.
- `GET /api/events` serves the persisted NDJSON event snapshot.

In v0.11.7 the UI remained a polling frontend over generated status files; v0.11.8 adds the SSE status seam described below.

## v0.11.8 server-sent status stream

The generated controller now exposes `GET /api/events/stream`. Browsers receive `status` events when `status.json` changes and heartbeat comments at the configured interval. `dashboard.js` uses `EventSource` for live updates while its existing polling loop remains active as a compatibility and recovery path.

This is still an operator-reviewable standard-library scaffold. Place it behind a TLS reverse proxy and production identity layer before exposing it outside a trusted network. SSE does not add secret storage, live SSH/Docker probes, durable event replay or daemon supervision.

## v0.12.0 production deployment baseline

The controller now defaults to TLS and authenticated API access. Copy `controller.tokens.example.json` to a root-owned location, replace placeholders with SHA-256 digests of independently generated random tokens, and point `CFD_CONTROLLER_TOKEN_FILE` to it. Certificate and key paths are supplied through `CFD_CONTROLLER_TLS_CERT` and `CFD_CONTROLLER_TLS_KEY`; missing values stop startup rather than silently falling back to plaintext.

Viewer tokens can read status, logs, history and SSE. Operator tokens can also issue control directives. `/api/probes` runs the generated access-check script with a timeout and bounded captured output. Audit records rotate when the configured byte limit is reached.

For persistent operation, install the generated systemd unit only after adapting paths, user/group ownership and `ReadWritePaths`. The Caddy example demonstrates TLS termination, SSE flushing and response security headers. Bind the direct controller listener to loopback when using a reverse proxy.
