# Release 0.14.0 - Phase 1 common HPC runtime completion

## Highlights

v0.14.0 closes the final three checklist items in the common HPC runtime phase:

- explicit NUMA node/thread placement;
- SYCL sparse matrix-vector and Krylov execution;
- MPI-owned sparse-vector exchange integrated into Krylov iterations.

The overall tracker moves from 580/719 to 583/719 capabilities, and Phase 1 moves from 31/34 to 34/34.

## NUMA/thread placement

`cfd/core/numa.hpp` now exposes topology discovery and deterministic placement plans. Linux discovery intersects `/sys/devices/system/node/node*/cpulist` with the calling process's current CPU affinity, which respects container/cpuset/MPI-launcher restrictions. Compact, spread and explicit-node plans map workers to concrete CPUs.

Linux affinity uses dynamically allocated CPU masks rather than fixed `cpu_set_t`, avoiding the 1024-CPU ceiling documented for fixed glibc masks. Placement-aware first touch binds each worker before it initializes its contiguous range.

## SYCL sparse linear algebra

`cfd::core::SyclCsrLinearAlgebra` keeps CSR structure/value arrays resident on the selected SYCL device and provides:

- CSR SpMV;
- device dot-product reductions;
- conjugate gradient for SPD systems;
- explicit USM host/device transfer boundaries;
- device capability checking and device-name diagnostics.

This is a common linear-algebra backend, not a solver-family-specific GPU implementation.

## MPI distributed Krylov

`cfd::distributed::MpiDistributedCsrOperator` builds on `DistributedCsrPartition`. During setup it gathers contiguous row ownership, determines the owner of each halo column, exchanges requested global indices, and caches response offsets. Each multiply then exchanges only needed vector entries and passes the refreshed halo to the existing local distributed CSR multiply.

`mpi_distributed_conjugate_gradient(...)` integrates that operator directly into each Krylov iteration and performs communicator-wide scalar reductions for dot products and convergence norms.

## Validation

Performed in the packaging environment:

- GCC/OpenMP compile and execution of `test_v0140_hpc_runtime.cpp`: passed.
- C++ syntax/interface check of `src/distributed/mpi_sparse.cpp`: passed against the MPI APIs used.
- C++ syntax/interface check of the extended `tests/test_mpi.cpp`: passed against the MPI APIs used.
- C++ syntax/interface check of `src/sycl/sparse_linalg_sycl.cpp`: passed against the SYCL API surface used.

Environment limitation: this host does not provide an MPI implementation or AdaptiveCpp/SYCL runtime, so real multi-rank execution and accelerator execution were not performed here. The CMake test matrix contains conditional runtime tests for those configurations.
