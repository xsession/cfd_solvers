# Deep-search: performance techniques from x3d2, AeroFVM, and cfd3d

This review compares three public repositories against the current `cfd_solvers`
architecture. The goal is to extract transferable performance techniques while
keeping the unified codebase readable, portable, and validation-backed.

Reviewed repositories:

- [xcompact3d/x3d2](https://github.com/xcompact3d/x3d2)
- [ankurjain123AEGit/AeroFVM](https://github.com/ankurjain123AEGit/AeroFVM)
- [chrismile/cfd3d](https://github.com/chrismile/cfd3d)

## Executive result

The strongest ideas are complementary rather than interchangeable:

| Repository | Best transferable idea | Main caution |
|---|---|---|
| x3d2 | Data-oriented grouped systems, cache blocking, fused RHS/solver passes, and communication-aware decomposition | The algorithm is specialized to compact finite differences and must not be copied into arbitrary unstructured operators |
| AeroFVM | Small, explicit CSR implementation with solver-level residual/error benchmarking | Primarily an educational baseline: the CSR SpMV and solver loops are serial and use wide `size_t` indices |
| cfd3d | One physical solver exposed through CPU, MPI, CUDA, and OpenCL execution paths; explicit face-packing and backend performance scripts | Several backend implementations duplicate physics and indexing logic, which increases maintenance and consistency risk |

## x3d2 findings

### 1. DistD2-TDS is a data-layout technique, not only a solver

The x3d2 theory documentation describes grouping many independent tridiagonal
systems so that the same position across systems is contiguous in memory. This
creates predictable access for CPU vectorization and GPU thread-level parallelism.
It also keeps intermediate data close enough for cache/shared-memory reuse.

The important transferable pattern is:

1. classify the repeated small systems;
2. pack them into a structure designed for the target vector width;
3. fuse right-hand-side construction with the first solver pass;
4. communicate only the reduced boundary state between neighboring ranks;
5. avoid moving the full field between solver phases.

References:

- [x3d2 theory: DistD2-TDS](https://github.com/xcompact3d/x3d2/blob/main/docs/source/user/theory.rst)
- [CPU distributed execution](https://github.com/xcompact3d/x3d2/blob/main/src/backend/omp/exec_dist.f90)
- [CUDA distributed execution](https://github.com/xcompact3d/x3d2/blob/main/src/backend/cuda/exec_dist.f90)
- [TDS operator state](https://github.com/xcompact3d/x3d2/blob/main/src/tdsops.f90)

### 2. 2D decomposition is valuable only when the transpose plan is explicit

x3d2 uses 2D pencil decomposition for FFT-based Poisson solves. The code keeps
separate physical/spectral layouts, persistent work arrays, and backend-specific
transpose paths. This improves scalability beyond slab decomposition but also
makes communication costs visible and testable.

For `cfd_solvers`, the practical lesson is to represent a distributed field as
an explicit layout state rather than hiding every redistribution behind a generic
copy. A future distributed FFT or structured pressure solver should expose:

- layout state and ownership ranges;
- planned send/receive counts and displacements;
- reusable transpose buffers;
- measured communication bytes and synchronization time.

References:

- [x3d2 OpenMP FFT path](https://github.com/xcompact3d/x3d2/blob/main/src/backend/omp/poisson_fft.f90)
- [x3d2 CUDA FFT path](https://github.com/xcompact3d/x3d2/blob/main/src/backend/cuda/poisson_fft.f90)

### 3. Persistent device allocation is preferable to per-step ownership changes

The CUDA backend allocates device coefficient/work arrays during setup and reuses
them across derivative, tridiagonal, transpose, and FFT operations. This is a
good match for the existing resident SYCL work in `cfd_solvers`.

The project should continue to move toward a lifecycle of:

`build plan -> allocate once -> execute many steps -> release plan`

rather than creating temporary accelerator buffers from inside time-step kernels.

## AeroFVM findings

AeroFVM is most useful as a compact correctness and benchmarking reference.
Its CSR class has a clear finalized representation, reserves storage before
assembly, and exposes solver-independent benchmark metrics: runtime, iterations,
residual, L2 error, L-infinity error, convergence rate, and status.

Those benchmark semantics are worth adopting broadly in `cfd_solvers`, especially
for comparing CPU, OpenMP, SYCL, and MPI variants. The implementation itself
should not be copied as a performance baseline: its CSR matrix uses
`std::size_t` indices and its SpMV loops are serial, so it is primarily a readable
reference implementation.

References:

- [AeroFVM CSR interface](https://github.com/ankurjain123AEGit/AeroFVM/blob/main/core/include/Sparse/SparseMatrix.h)
- [AeroFVM CSR implementation](https://github.com/ankurjain123AEGit/AeroFVM/blob/main/core/src/Sparse/SparseMatrix.cpp)
- [AeroFVM benchmark suite](https://github.com/ankurjain123AEGit/AeroFVM/blob/main/benchmark_solver/benchmark_solver.cpp)

## cfd3d findings

### 1. Separate execution backends are practical for radically different APIs

cfd3d keeps C++/OpenMP, MPI, CUDA, and OpenCL solver implementations and chooses
the execution path at runtime. This demonstrates a useful product boundary for
`cfd_solvers`: physics and discretization should be shared, while execution plans
may specialize memory layout, synchronization, and kernel launch details.

For the unified project, the safer version is a common operator contract with
backend-specific kernels, not four independent copies of the numerical model.

### 2. Face-wise MPI packing is simple and predictable

The MPI helper code explicitly packs the six neighbor faces into contiguous
buffers, performs `MPI_Sendrecv`, and unpacks into halo layers. This is easy to
profile and works for staggered fields with different extents. The next evolution
for `cfd_solvers` is to precompute those pack/unpack maps and use nonblocking or
persistent requests where overlap is beneficial.

Reference:

- [cfd3d MPI halo helpers](https://github.com/chrismile/cfd3d/blob/master/src/CfdSolver/Mpi/MpiHelpers.cpp)
- [cfd3d MPI index definitions](https://github.com/chrismile/cfd3d/blob/master/src/CfdSolver/Mpi/DefinesMpi.hpp)

### 3. GPU portability requires a deliberately limited kernel contract

cfd3d’s CUDA path supports Jacobi for the pressure solve, while the OpenMP and
OpenCL paths expose parallel stencil variants. This is a useful reminder that a
portable GPU backend should begin with a small set of regular, bandwidth-bound
kernels: pull/push updates, stencil application, reductions, and halo copies.
Advanced sparse or irregular kernels should have a correctness-first fallback.

References:

- [cfd3d CUDA SOR/Jacobi solver](https://github.com/chrismile/cfd3d/blob/master/src/CfdSolver/Cuda/SorSolverCuda.cu)
- [cfd3d OpenCL SOR/Jacobi kernel](https://github.com/chrismile/cfd3d/blob/master/src/CfdSolver/Opencl/SorSolverOpencl.cl)
- [cfd3d CPU solver](https://github.com/chrismile/cfd3d/blob/master/src/CfdSolver/Cpp/SorSolverCpp.cpp)

## Implemented in this pass

`DistributedCsrPartition::multiply()` previously performed a binary search over
the sorted halo-column list for every off-rank coefficient on every SpMV. The
constructor now builds a compact per-coefficient reference map once. The hot loop
uses direct local/halo indexing and runs through the existing OpenMP-aware
`parallel_for` utility.

This preserves the public API and numerical results while removing repeated
`O(log(H))` index lookups from each multiply. The MPI sparse operator now also
exposes an overlap-capable `begin_multiply()`/`finish_multiply()` lifecycle:
the packed value exchange starts first, independent local work can run while
MPI progresses, and the finish phase waits, unpacks the precomputed receive
map, and applies the local operator. The existing `multiply()` call remains a
synchronous wrapper for compatibility.

## Ranked next performance phases

1. Add benchmark counters for SpMV bytes, FLOPs, halo bytes, and synchronization
   time; report effective GB/s and scaling rather than wall time alone.
2. Extend the new nonblocking sparse lifecycle into production solver loops and
   measure communication/compute overlap on real MPI ranks.
3. Add a structured stencil operator for regular FVM/LBM pressure systems so
   regular grids do not pay general CSR indirection.
4. Fuse residual construction with stencil/SpMV application where the same field
   is scanned twice.
5. Add vector-width-aware blocked storage for batched 1D solves, inspired by
   DistD2-TDS, with scalar fallback and verification against the reference solver.
6. Extend the portable accelerator contract to CPU/OpenMP and SYCL first, then
   evaluate an OpenCL compatibility path without duplicating physics code.
7. Add performance-regression thresholds only after collecting stable baselines
   on representative CPU, NVIDIA, AMD, and Intel systems.

## Boundaries

The review does not claim that any repository's reported performance transfers
directly to this project. GPU, MPI, and cache improvements must be validated on
the target hardware and on problem sizes large enough to amortize setup costs.
