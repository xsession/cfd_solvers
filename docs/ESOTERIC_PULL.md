# Esoteric-Pull implementation notes

This project implements the published Esoteric-Pull access pattern independently from the algorithm description in:

Moritz Lehmann, *Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs*, Computation 10(6), 92 (2022), DOI 10.3390/computation10060092.

No FluidX3D source code is copied or translated.

## Why this pattern

A conventional pull LBM uses two population fields:

```text
old populations -> pull + collide -> new populations
```

That is simple and parallel, but DDF storage is `2 * Q * N` values. LBM is usually memory-bandwidth bound, and the second full population grid also halves the maximum lattice size that fits in a fixed memory budget.

Esoteric Pull keeps only `Q * N` population values. Opposite directions are paired in the descriptor:

```text
0: rest
1,2: opposite pair
3,4: opposite pair
...
```

D2Q9 therefore has four moving pairs, D3Q19 nine, and D3Q27 thirteen.

## Parity-controlled ownership

For each opposite pair `(a,b)`, `b` is used to identify the periodic neighbor. The physical memory slot used for a logical direction alternates with time-step parity.

At one parity, the logical `a` value is pulled from the neighbor while `b` is local. After collision, the resulting values are written back with their pair orientation exchanged. At the other parity, the slot roles reverse.

The important invariant is:

**each thread stores to exactly the same pair of physical locations from which that logical pair was loaded.**

The neighbor map on a periodic Cartesian grid is bijective, so different lattice nodes do not write the same location. That is what makes the single-grid update thread-safe on both OpenMP and GPU kernels.

## Hot-loop sequence

For one lattice node:

1. reconstruct all logical populations using current parity;
2. accumulate density and momentum;
3. calculate velocity;
4. apply BGK/SRT collision in registers/local storage;
5. store each opposite pair using current parity;
6. toggle parity once the full kernel completes.

Streaming and collision are therefore one fused memory pass over the DDFs.

## Descriptors

`include/cfd/solvers/lbm/descriptors.hpp` currently defines:

- D2Q9, used primarily to compare the in-place path with the older two-grid implementation;
- D3Q19, intended as the default 3-D performance baseline;
- D3Q27, providing the full tensor-product velocity set.

The compile-time regression tests verify opposite-pair symmetry, weight normalization and second-order isotropy moments.

## Boundaries and forcing

The Phase-2 kernel is intentionally periodic. This keeps the first in-place implementation small enough to validate rigorously against the existing two-grid solver.

Next additions should be made one at a time with reference tests:

1. halfway bounce-back;
2. body-force formulation (Guo forcing first);
3. velocity/pressure inlet/outlet treatment;
4. interpolated curved-wall boundary;
5. TRT, then MRT/cumulant collision options.

Adding these after the parity kernel is validated is safer than debugging streaming, boundaries and forcing simultaneously.
