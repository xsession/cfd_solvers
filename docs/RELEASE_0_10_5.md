# Release 0.10.5 - Serialized PIC transport and MPI bridge contract

## Added

- ABI-local serialized 3-D PIC transport envelopes for particle migration and scalar guard-cell messages.
- Destination-rank serialized exchange planning with message and byte counts.
- Serialization/deserialization helpers with version/magic checks and malformed-envelope rejection.
- Guarded `CFD_HAS_MPI` allgather transport wrapper that can transport the same serialized envelopes when MPI is available.
- `particle-serialized-transport3d` CLI smoke case.
- `cfd-v0105-serialized-transport-tests` focused regression target.
- Multi-species plasma mass-action chemistry utilities with charge-density diagnostics.
- Paschen gas-breakdown and parallel-plate multipactor threshold estimators.
- `particle-plasma-chemistry` and `particle-breakdown-threshold` CLI smoke cases.
- `cfd-v0105-plasma-tests` focused regression target.

## Validation

- Full debug/CI-speed CPU CTest matrix: 85/85 passed in this environment.
- Focused ASan+UBSan+leak detection for the v0.10.5 serialized transport slice: passed.
- Focused ASan+UBSan+leak detection for the v0.10.5 plasma/breakdown slice: passed.

## Limitations

- The MPI wrapper is compile-time guarded and was not runtime-validated here because no MPI runtime/compiler is available in this environment.
- The transport is still an envelope/serialization layer; distributed PIC time stepping, optimized neighbour-only MPI, GPU particle bins/kernels and production local high-order current deposition remain open.
- Plasma chemistry and breakdown workflows are analytical/reference baselines, not a full multipactor/corona/PIC discharge solver yet.
