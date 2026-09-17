# Release 0.10.3 - communicator-ready PIC exchange packets

## Added

- Communicator-ready 3-D particle migration messages.
- Retained-domain particle buckets plus source-domain/source-index tracking.
- Communicator-ready scalar guard-cell messages with destination padded coordinates.
- Apply helpers that reconstruct migrated particle buckets and ghost-padded scalar guard blocks from messages.
- `particle-comm-exchange3d` CLI smoke case.
- `cfd-v0103-comm-exchange-tests` focused regression target.

## Validation

- Full debug/CI-speed CPU CTest matrix: 78/78 passed.
- Focused ASan + UBSan + leak detection binary for the new communicator-message slice: passed.

## Scope note

This release does not perform MPI communication. It defines and validates the deterministic payload contract that the next MPI transport layer can send and receive.
