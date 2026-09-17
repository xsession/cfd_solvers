# Release 0.11.8 - server-sent controller status events

## Summary

This release adds a dependency-free server-sent event (SSE) seam to the generated multi-server controller. The browser dashboard receives status changes over a long-lived HTTP response and retains periodic polling as a compatibility and reconnect fallback.

## Implemented

- `RemoteControllerConfig::expose_event_stream` to include or omit the live route.
- `RemoteControllerConfig::event_stream_heartbeat_seconds` with non-zero validation.
- `GET /api/events/stream` using `text/event-stream`, per-connection event IDs and heartbeat comments.
- Change-aware `status` events generated from the same `status.json` contract used by polling clients.
- Dashboard `EventSource` client with automatic browser reconnect behavior and visible live/reconnecting state.
- OpenAPI and route-table metadata for the SSE endpoint.
- Focused regression target `cfd-v0118-controller-sse-tests`.

## Validation

- Direct C++ compilation and focused regression execution are supported without optional dependencies.
- Generated Python controller syntax and HTTP/SSE behavior are covered by release validation.
- The full CTest matrix should be run in a CMake-equipped build environment before tagging.

## Limitations

SSE is a one-way status stream, not a bidirectional WebSocket control plane. The generated controller still lacks TLS termination, user accounts, managed secrets, live Docker/SSH probes, durable event history and production daemon supervision.
