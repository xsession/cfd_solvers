# Release 0.10.4 - in-memory PIC transport contract for future MPI

## Added

- Logical-rank topology mapping for 3-D PIC subdomains.
- Rank-addressed transport envelopes for particle migration and scalar guard-cell payloads.
- Deterministic `InMemoryPicTransport3D` inbox transport.
- Complete in-memory exchange helper that reproduces particle migration and scalar guard exchange through transport envelopes.
- `particle-transport3d` CLI smoke case.
- `cfd-v0104-transport-tests` focused regression target.

## Validation

- Full debug/CI-speed CPU CTest matrix: 80/80 passed.
- Focused ASan + UBSan + leak detection for the new v0.10.4 transport slice: passed.

## Tracker

- Overall: 527/669 validated capabilities, 78.8%.
- Phase 12: 61/81 validated capabilities, 75.3%.

## Limitation

This release still uses serial in-memory transport. Real MPI send/receive integration, rank-local distributed PIC stepping, GPU particle bins and guard-cell device exchange remain future work.
