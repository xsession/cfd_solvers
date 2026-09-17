# Optimization plan and performance gates

Correctness precedes optimization. Each optimization must keep a reference path and a numerical-error threshold.

## Cross-cutting priorities

### Linear algebra

- [x] matrix-free CG.
- [x] CSR storage/builder.
- [x] Jacobi-PCG.
- [x] Jacobi-BiCGStab.
- [x] restarted GMRES for strongly nonsymmetric systems.
- [x] ILU(0).
- [ ] AMG adapter.
- [ ] block matrices for coupled vector/species systems.
- [ ] SYCL SpMV.
- [ ] distributed SpMV.
- [ ] mixed-precision refinement.

### Memory and allocation

- [x] aligned LBM SoA.
- [x] one-grid LBM streaming.
- [ ] reusable FVM face/cell scratch workspaces.
- [ ] reusable chemistry reaction scratch buffers.
- [ ] monotonic/pool allocator for per-step temporary fields.
- [ ] topology-dependent sparsity cache.
- [ ] field packing/reordering for cache locality.

### CPU

- [x] OpenMP cell/rank execution.
- [x] small-loop OpenMP serial cutoff and nested-team guard; current CPU baseline uses a conservative 2048-iteration threshold.
- [ ] explicit SIMD for fixed-size Vec3 and stencil loops.
- [ ] NUMA first-touch and rank/thread affinity helpers.
- [ ] graph/mesh reordering benchmark (RCM/METIS order).
- [ ] thread-scaling CI on a fixed CPU runner.

### GPU

- [x] LBM SYCL kernels.
- [x] distributed SYCL LBM path.
- [ ] FVM interpolation/gradient/divergence kernels.
- [ ] sparse algebra kernels.
- [ ] FEM assembly/matrix-free kernels.
- [ ] FDTD kernels.
- [ ] optics batched rays.
- [ ] electrochemistry transport kernels.
- [ ] device-memory pool and event dependency graph.

### Communication

- [x] compact LBM halos.
- [x] persistent MPI requests.
- [x] communication/interior overlap for distributed LBM.
- [ ] generic mesh ghost exchange.
- [ ] sparse matrix off-rank halo schedule.
- [ ] GPU-aware FVM/FEM communication.
- [ ] topology-aware rank placement.

## Numerical optimization gates

Every performance change must report:

- runtime/throughput;
- allocated bytes or peak memory;
- conservation error where applicable;
- physical/manufactured solution error;
- CPU/GPU parity if both paths exist;
- deterministic or bounded nondeterministic tolerance;
- before/after benchmark command and hardware metadata.

No fast-math or reduced-precision path becomes default until its error envelope is explicitly tested.
