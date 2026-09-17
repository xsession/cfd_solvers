# Release 0.12.0 - production controller deployment baseline

## Summary

This release closes the controller scaffold's previously listed deployment gaps with secure-by-default generated artifacts. It adds TLS, authenticated roles, bounded audit history, live access probes and service supervision without embedding credentials or certificates in generated source.

## Implemented

- TLS 1.2+ standard-library server wrapping; startup fails when enabled certificate/key paths are absent.
- External JSON token file containing only SHA-256 token digests and viewer/operator/admin roles.
- Constant-time digest comparison and authentication on protected API, SSE and control routes.
- Dashboard session-scoped token handling and authenticated EventSource connection.
- Bounded `events-history.ndjson` audit history with single-file rotation.
- Authenticated `/api/probes` execution of the generated SSH/Docker/filesystem access checks with a fixed timeout and bounded output.
- Hardened `cfd-controller.service` systemd unit.
- `Caddyfile.example` with TLS, streaming flush and security headers.
- Credential-free `env.example` and explicit token template.
- Focused `cfd-v0120-production-controller-tests` regression target.

## Security boundary

The project generates and validates deployment artifacts; it does not manufacture trust. Operators remain responsible for certificate issuance and rotation, high-entropy token generation, filesystem ownership, SSH host verification, Docker authorization, network policy, centralized logs and host patching. For enterprise SSO, place the controller behind an identity-aware proxy and disable direct network access to its listener.

## Validation

- Focused v0.11.6-v0.12.0 controller regressions compile and execute directly with C++20.
- Generated controller passes Python bytecode compilation.
- The full CTest matrix remains required in a CMake-equipped release environment.
