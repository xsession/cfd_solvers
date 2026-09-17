# Release 0.11.6 - multi-server controller API scaffold

## Highlights

- Adds a generated dependency-free Python controller service over multi-server supervision artifacts.
- Provides health/status/case/log-tail read routes and an optional token-gated control route.
- Persists OpenAPI JSON, status JSON, routes TSV, env example, README and a root launch script.
- Supports read-only controller planning for dashboard-only deployments.

## Validation

- Full debug/CI-speed CPU CTest matrix: 109/109 passed.
- Focused ASan + UBSan + leak detection: `cfd-v0116-multiserver-controller-tests` passed.

## Limitations

This is an offline-generated controller scaffold, not a production service. It does not manage SSH keys, store secrets, provide TLS, stream logs over WebSockets, probe live Docker daemons at runtime or replace the generated remote execution scripts.
