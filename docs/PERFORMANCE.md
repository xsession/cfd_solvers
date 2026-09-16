# Performance plan

## Implemented in Phase 2

- 64-byte aligned host arrays.
- Structure-of-arrays LBM populations as `[population][cell]`.
- Static OpenMP work sharing.
- D2Q9, D3Q19 and D3Q27 compile-time descriptors.
- Single-grid Esoteric-Pull-style streaming implemented independently from the published algorithm.
- One-Step-Pull two-grid reference solver for numerical and performance comparison.
- Fused load/stream + macroscopic reconstruction + BGK collision + store kernel.
- Reuse of per-cell neighbor indices between the load and store halves of the in-place update.
- Optional device-resident SYCL population storage and single-grid kernels.
- Benchmark CLI with MLUPS, minimum DDF traffic estimate, memory footprint and mass drift.

The single-grid representation cuts population allocation from `2 * Q * cells * sizeof(T)` to `Q * cells * sizeof(T)`: exactly 50% for the DDF arrays. It does **not** guarantee a fixed speedup on every CPU/GPU; cache hierarchy, compiler, work-group mapping and device memory behavior matter.

## Benchmark methodology

Use a warmup interval before timing. Report:

- lattice and dimensions;
- compiler/backend/device;
- thread count or accelerator;
- update count and elapsed time;
- MLUPS;
- allocated DDF memory;
- absolute and relative mass drift;
- estimated minimum DDF GB/s.

The bandwidth estimate is `MLUPS * (2 * Q * sizeof(float)) / 1000`. It represents the idealized population load/store traffic only. It excludes cache-line effects, addressing, metadata, diagnostics and other traffic, so it must not be presented as measured DRAM bandwidth.

## CPU direction

- Profile before adding explicit SIMD intrinsics.
- Add NUMA-aware first touch and pinning for large dual-socket systems.
- Add blocked/brick traversal experiments for D3Q27 if cache/TLB behavior warrants it.
- Consider precomputed compact neighbor offsets only when their memory cost beats integer-address arithmetic.
- Keep OpenMP as the readable baseline; vendor-specific paths belong behind narrow hooks.

## GPU portability

AdaptiveCpp/SYCL is the preferred portable backend. The code keeps the DDF allocation on the device between time steps and performs a single kernel per BGK update. Future work:

- tune local/work-group sizes;
- subgroup-aware reduction/diagnostics;
- multi-device halos and compute/transfer overlap;
- measure occupancy/register pressure for D3Q27 local population arrays;
- vendor-specific fast paths only when benchmark evidence justifies them.

## Reduced storage

A future compact storage path may use 16-bit or custom encodings for populations while performing collision and global reductions in FP32/FP64. It must be gated by:

- mass/momentum conservation;
- velocity/density error against FP32 reference;
- benchmarked memory-bandwidth benefit;
- explicit case-level opt-in.

## Standard problems

Performance changes should ultimately be checked on Taylor-Green vortex, lid-driven cavity, Poiseuille/channel flow and multiphase cases, not only uniform periodic grids.

## Phase 3B distributed communication

The distributed LBM path packs only populations that actually cross each neighbour boundary. For a 64^3 local brick, exact message-layout accounting gives:

- D3Q19: 1,926,752 bytes/step for full-cell halos versus 494,592 bytes/step selective (`74.33%` reduction);
- D3Q27: 2,738,016 bytes/step versus 893,984 bytes/step selective (`67.35%` reduction).

These are exact payload counts, not network-throughput measurements. Persistent MPI requests and reusable staging buffers remove per-step request/buffer construction overhead. The scaling driver records halo messages and bytes alongside MLUPS so communication growth can be separated from kernel throughput.

## Phase 4A FVM performance direction

The polyhedral operators iterate through cell-owned adjacency rather than atomically accumulating by face. That preserves the owner/neighbour formulation while exposing independent cells to OpenMP. The pressure projection uses a shared matrix-free conjugate-gradient loop rather than a fixed-count Jacobi method; on the current 96x80 periodic projection smoke case the pressure correction converges in 7 iterations. This iteration count is a case-specific regression observation, not a general solver-performance guarantee.

Next performance work for FVM should prioritize sparse adjacency locality, reusable face interpolation data, matrix-free momentum/pressure operators, preconditioners, and MPI partition halos before considering architecture-specific intrinsics.

## Phase 4B transient-flow direction

The new face-scheme operators remain cell-owned and OpenMP-parallel, so conservative flux accumulation needs no atomics. The transient staggered solver keeps pressure application matrix-free and reuses the shared CG workspace across time steps. On the 64x64 Taylor-Green smoke case, 100 steps to `t=0.01` complete in roughly a tenth of a second on the development runner, but this is only a local regression observation. Future performance work must focus on collocated sparse adjacency locality, face-coefficient reuse, pressure preconditioning, and distributed unstructured halos rather than optimizing the periodic staggered oracle.
