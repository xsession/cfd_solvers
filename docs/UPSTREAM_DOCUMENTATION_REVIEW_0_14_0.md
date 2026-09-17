# Upstream documentation review - v0.14.0 common HPC runtime

This is a clean-room documentation review used to close the remaining Phase-1 runtime items. No upstream implementation source was copied.

## OpenMP affinity and Linux CPU affinity

Reviewed:

- OpenMP 6.0/5.2 specifications and affinity examples: https://www.openmp.org/specifications/
- OpenMP affinity examples: https://www.openmp.org/wp-content/uploads/openmp-examples-6-0-1.pdf
- Linux `sched_setaffinity(2)`: https://man7.org/linux/man-pages/man2/sched_setaffinity.2.html

Transferable lessons applied:

1. OpenMP affinity is expressed in terms of places and binding policy, so the framework should expose a placement policy rather than assuming first-touch alone is enough.
2. Linux CPU affinity is per-thread; binding the calling OpenMP worker before first touch provides an explicit CPU-placement boundary without mutating global `OMP_*` environment variables at runtime.
3. The process may already be restricted by cpusets/cgroups or an MPI launcher, so topology discovery must intersect NUMA-node CPU lists with the process's current allowed affinity mask.
4. Fixed `cpu_set_t` is limited on very large systems. The Linux manual recommends dynamically allocated CPU sets when the kernel affinity mask can exceed 1024 CPUs; v0.14.0 uses `CPU_ALLOC`/`CPU_ALLOC_SIZE` accordingly.

Implemented:

- `NumaPlacementPolicy::{inherit,compact,spread,explicit_node}`;
- `/sys/devices/system/node/node*/cpulist` topology discovery on Linux with portable fallback;
- deterministic per-thread placement planning;
- explicit current-thread binding and OpenMP-team application;
- placement-aware first touch.

## SYCL 2020 USM and reductions

Reviewed:

- SYCL 2020 specification: https://registry.khronos.org/SYCL/specs/sycl-2020/pdf/sycl-2020.pdf
- Khronos SYCL USM reference: https://github.khronos.org/SYCL_Reference/iface/usm_allocations.html
- Khronos SYCL reduction reference: https://github.khronos.org/SYCL_Reference/iface/reduction-variables.html

Transferable lessons applied:

1. Device USM is appropriate for persistent CSR arrays and iterative vectors where explicit host/device transfers are desirable.
2. Queue `memcpy` is the explicit host/device boundary; the iterative loop should avoid copying whole vectors back each Krylov iteration.
3. Standard SYCL reductions provide the portable scalar dot-product primitive required by CG without introducing a vendor-specific reduction library.
4. Device capability checks should fail clearly when the selected device does not support the required USM allocation modes.

Implemented:

- persistent device CSR row offsets, column indices and values;
- SYCL CSR SpMV;
- SYCL dot-product reductions;
- device-resident CG vector updates and convergence loop;
- host-facing parity/solve regression when `CFD_ENABLE_SYCL=ON`.

## MPI 4.1 sparse communication

Reviewed:

- MPI 4.1 standard: https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report.pdf

Transferable lessons applied:

1. Sparse numerical operators should exchange only entries referenced by off-rank columns rather than gathering a full distributed vector.
2. Ownership metadata is setup-time state; global column requests can be exchanged once and cached for repeated Krylov iterations.
3. Nonblocking collectives provide a standard communication primitive while keeping solver code independent of MPI implementation details.
4. Krylov scalar quantities are global mathematical values and therefore require communicator-wide reductions.

Implemented:

- contiguous row-ownership discovery with validation;
- setup-time exchange of requested global sparse-vector columns;
- cached local response offsets and halo-slot mapping;
- `MPI_Ialltoallv` value refresh inside distributed CSR multiply;
- distributed CG with global dot products/residual norm through `MPI_Allreduce`;
- MPI regression solving a distributed SPD tridiagonal system.

## Remaining runtime work beyond Phase 1

Phase 1 is complete as tracked, but later phases still contain performance/scale work such as multi-GPU-per-rank execution, adaptive repartitioning, hardware CI on NVIDIA/AMD/Intel accelerators, and distributed FEM/FDTD/PIC features. Those remain phase-specific items rather than gaps in the common runtime baseline.
