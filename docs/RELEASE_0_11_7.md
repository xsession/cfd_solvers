# Release 0.11.7 - generated multi-server controller dashboard

## Summary

This release adds a generated browser dashboard on top of the v0.11.6 dependency-free multi-server controller scaffold. The dashboard is still static/generated and offline-testable, but it gives the supervision artifacts a first operator-facing UI without requiring Node.js, FastAPI, WebSockets or a frontend build step.

## Implemented

- `RemoteControllerConfig::dashboard_title` and `dashboard_refresh_seconds` for operator-facing dashboard customization.
- `RemoteControllerConfig::expose_event_snapshot` to include or omit the generated event-snapshot route.
- Generated controller assets:
  - `controller/index.html`
  - `controller/dashboard.js`
  - `controller/dashboard.css`
  - `controller/events.ndjson`
- Controller static routes:
  - `GET /`
  - `GET /dashboard.js`
  - `GET /dashboard.css`
  - `GET /api/events`
- Dashboard JavaScript that polls `/api/health` and `/api/status`, renders summary cards and case rows, fetches stdout/stderr log tails, and posts token-authorized control actions.
- New CLI smoke case: `particle-multiserver-dashboard`.
- New focused regression target: `cfd-v0117-controller-dashboard-tests`.

## Validation

- Full CPU CTest matrix: 111/111 passed.
- Focused ASan + UBSan + leak-detection executable for the dashboard/controller slice passed.

## Limitations

The dashboard is a generated, static polling UI. It does not yet provide TLS, user accounts, WebSocket/SSE live streams, live Docker/SSH probing, stored secrets, or a persistent production controller daemon.
