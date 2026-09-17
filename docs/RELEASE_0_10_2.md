# Release 0.10.2 - PIC domain decomposition and scalar guard exchange foundation

## Added

- Cartesian 3-D PIC domain descriptors with physical extents and uneven global-cell block coverage.
- Serial particle migration bucketing between subdomains.
- Periodic wrapping and nonperiodic outside-particle reporting for migration planning.
- Serial scalar guard-cell exchange that builds ghost-padded local blocks from neighbouring subdomains.
- `particle-domain-exchange3d` CLI smoke case.
- `cfd-v0102-domain-exchange-tests` focused regression target.

## Validation

- Full debug/CI-speed CPU CTest matrix: 76/76 passed.
- Focused ASan + UBSan + leak detection binary for the new domain-exchange slice: passed.

## Tracker

- Overall: 524/666 capabilities complete (78.7%).
- Phase 12 CST-class expansion: 58/78 complete (74.4%).

## Limitations

This release does not implement MPI particle exchange, communicator-aware field halo exchange, GPU particle bins or load balancing. It provides the serial decomposition, migration and guard-cell copy contracts that those later implementations should reuse.
