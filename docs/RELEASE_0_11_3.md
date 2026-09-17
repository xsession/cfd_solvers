# Release 0.11.3 - multi-server campaign planning and Docker deployment scaffold

This checkpoint extends the SU2/csauto-inspired workflow layer from local/scheduler control contracts into deterministic deployment planning. It does not require real SSH or Docker during validation; instead it generates auditable server assignments, remote command scripts and a Docker/Compose deployment scaffold.

## Added

- `CampaignServerDescriptor` for multi-server campaign placement.
- `CampaignServerRuntime` for native vs Docker worker placement.
- `plan_multi_server_campaign(...)` with:
  - capacity-aware server slots;
  - online-server filtering;
  - required-tag filtering;
  - deterministic pending-case assignment;
  - native or Docker solver command planning per assigned case.
- Remote command plans for:
  - remote case-directory creation;
  - rsync-style case synchronization;
  - native or Docker worker invocation.
- `write_multiserver_plan_files(...)`, producing:
  - `multiserver_assignments.tsv`;
  - executable `multiserver_commands.sh`.
- Docker deployment generator producing:
  - `Dockerfile`;
  - `compose.yaml`;
  - `.dockerignore`;
  - `deploy/deploy.sh`;
  - `deploy/env.example`;
  - `deploy/servers.example.tsv`;
  - `deploy/README.md`.
- CLI smoke case `particle-multiserver-deploy`.
- Focused regression target `cfd-v0113-multiserver-deploy-tests`.

## Validation

- Full CPU CTest matrix: `103/103 passed`.
- Focused ASan + UBSan + leak detection for the v0.11.3 multi-server/deploy slice: passed.
- Integration tracker: `558/699 = 79.8%` overall.
- Phase 9 workflow/UX: `28/30 = 93.3%`.

## Limits

This is a deployment-planning and generated-script layer. It does not yet supervise real remote SSH sessions, validate Docker daemons on remote hosts, manage secrets, stream remote logs, or expose a web dashboard for multi-server operations. Those are the next workflow steps.
