# Phase 3A results

## Implemented

- aspect-aware Cartesian process-grid selection;
- uneven 2-D/3-D brick decomposition;
- rank/coordinate and periodic/nonperiodic neighbour mapping;
- one-cell component-major halo storage;
- 26-neighbour face/edge/corner packing and unpacking;
- virtual-rank exchange backend for deterministic local validation;
- D3Q19/D3Q27 distributed two-grid pull LBM blocks;
- Guo body acceleration in distributed collision kernels;
- OpenMP execution inside each distributed rank;
- split strict-interior / boundary-shell kernels for communication overlap;
- MPI Cartesian topology and node-local rank discovery;
- nonblocking MPI point-to-point halo exchange;
- deterministic local-rank -> visible-SYCL-device assignment;
- rank-local checkpoint/restart;
- `cfd-distributed` MPI CLI;
- four-rank MPI CI job and serial virtual-rank regression coverage.

## Numerical validation performed locally

The development environment has no MPI implementation or AdaptiveCpp installation. Therefore MPI transport and SYCL device execution cannot be run directly here.

The transport-independent distributed algorithm is executed locally through the virtual-rank backend. On an uneven `17 x 15 x 13` periodic domain decomposed into `2 x 2 x 2` bricks:

- D3Q19 distributed fields match the monolithic two-grid reference after nine steps within `2e-6` maximum absolute field error;
- D3Q27 distributed fields match the monolithic two-grid reference after nine steps within the same gate;
- global mass agrees with the monolithic reference;
- the test exercises face, edge and corner ghost regions;
- checkpoint reload reproduces populations exactly and remains bitwise identical after the following step.

The new distributed code is also included in the ASan+UBSan smoke path via a small D3Q27 virtual decomposition.

## MPI validation strategy

When `CFD_ENABLE_MPI=ON`, `cfd-mpi-tests` runs under four processes. Each rank executes the overlapped halo/collision path, then rank 0 gathers global macroscopic fields and compares D3Q19 and D3Q27 against the existing monolithic `PrecisionPullSolver` reference.

The GitHub Actions MPI job installs OpenMPI, builds MPI + OpenMP, and executes the complete CTest suite. In this development environment the MPI translation units were additionally syntax-checked against a minimal MPI API shim to catch C++ compile errors; that is not a substitute for a real MPI run and is not reported as runtime validation.

## Current performance limitation

Phase 3A prioritizes correctness and debuggability. Each neighbour message contains every population for the affected boundary cells. This over-communicates compared with a production LBM transport. The next optimization pass will:

1. pack only populations crossing each neighbour boundary;
2. reuse persistent host buffers;
3. consider MPI persistent requests / neighborhood collectives where beneficial;
4. add pinned host staging for non-GPU-aware MPI;
5. add direct device-buffer exchange when the MPI implementation is GPU-aware.
