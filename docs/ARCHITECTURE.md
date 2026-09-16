# Architecture

## Goal

One readable C++20 codebase, one data/execution model, multiple numerical engines.

The framework deliberately does **not** force LBM, finite volume, finite element, FDTD and geometric optics into one discretisation. They share the expensive infrastructure: memory ownership, field layouts, parallel execution, decomposition, sparse/vector operations, geometry, materials, I/O, diagnostics and accelerator selection.

## Layers

### 1. Core

- aligned allocations and SoA fields
- Cartesian and, later, unstructured meshes
- parallel loop/reduction primitives
- execution/device discovery
- domain decomposition and halo metadata
- common vector/matrix/sparse storage
- matrix-free conjugate-gradient linear solver
- material/property database
- units and dimensional checks in configuration code, not hot loops

### 2. Numerical kernels

- LBM: dense structured grids and bandwidth-bound streaming/collision
- FVM: owner/neighbour polyhedral faces, Gauss operators, face-flux reconstruction and pressure/velocity coupling
- FEM: element integration, assembly, sparse linear/nonlinear solves
- FDTD: staggered-grid field updates and boundary/material kernels
- Optics: batched ray/surface intersection, refraction/reflection and analysis

### 3. Physics modules

Numerical engines are composed with physics models instead of duplicated into a new executable for every combination. Examples: incompressible Navier-Stokes, compressible flow, heat transfer, elasticity, Maxwell, multiphase interface transport, reacting flow and optical materials.

### 4. Runtime

- CPU: OpenMP static work sharing and SIMD-friendly data layouts.
- GPU: AdaptiveCpp/SYCL single-source C++ kernels.
- Distributed: MPI Cartesian brick decomposition, 26-neighbour nonblocking halo exchange, OpenMP within each rank, and deterministic node-local rank/device mapping.

## Data-layout rule

Hot fields are SoA by default. LBM populations are laid out `[population][cell]`, enabling contiguous traffic for a population across neighbouring work-items. The production LBM path stores one population lattice; a two-lattice pull solver is retained as a reference/oracle. Multi-component continuum fields will similarly expose component-major views. AoS is reserved for low-volume metadata where readability dominates.

## Solver API direction

Every solver will converge on:

```cpp
struct Solver {
    void initialize(const Case&);
    StepResult step(TimeStep);
    Diagnostics diagnostics() const;
    void checkpoint(Writer&) const;
};
```

Execution policy and memory space are selected without changing physics configuration.

## Distributed/multi-GPU direction

The CPU/MPI reference path implements brick partitioning, 26-neighbour halos, strict-interior/boundary-shell overlap, and persistent selective crossing-population exchange. It intentionally keeps MPI out of the numerical collision code. Phase 3B also provides staged MPI+SYCL transport through pinned host buffers and an explicit opt-in direct device-buffer MPI path.

Accelerator progression:

1. keep the current MPI + OpenMP path as the correctness oracle;
2. assign one visible SYCL accelerator per node-local rank;
3. keep local populations device-resident;
4. stage compact halo buffers through pinned host memory when MPI is not GPU-aware;
5. use direct device-buffer MPI when supported;
6. add multiple devices per rank only after one-device-per-rank scaling is stable.

The correctness oracle can still exchange complete boundary cells, while the production distributed path uses selective crossing-population packing. Host and device paths share the same constexpr halo linearization so face/edge/corner message ordering cannot silently diverge.

## Finite-volume foundation

Phase 4A adds a topology-first `PolyMesh` with cell volumes/centres, owner/neighbour face connectivity, oriented face-area vectors, boundary patches and prebuilt cell-to-face adjacency. The first reusable operators are cell-parallel Gauss gradient/divergence and an orthogonal two-point Laplacian. A periodic staggered pressure-projection kernel uses the shared matrix-free CG solver as the first incompressible pressure/velocity coupling baseline.

This layer is deliberately below OpenFOAM-class feature parity: non-orthogonal correction, convection schemes/limiters, collocated Rhie-Chow coupling and SIMPLE/PISO/PIMPLE algorithms remain subsequent work.

## Precision policy

- Reference: FP64 where discretisation/conditioning needs it.
- Production CFD/LBM/FDTD: FP32 by default where validated.
- Reduced storage: optional 16-bit distribution/field storage with FP32 compute after solver-specific error validation.
- Never conflate reduced storage with reduced accumulation precision for global conservation metrics.

## Transient incompressible baseline

Phase 4B adds reusable linear/upwind face interpolation and conservative scalar/vector convective-divergence operators on `PolyMesh`. The first time-dependent incompressible reference remains staggered and periodic: conservative flux-form momentum convection, centered viscosity, and a matrix-free CG pressure projection. It is intentionally a correctness oracle for later collocated arbitrary-polyhedron pressure/velocity coupling rather than a shortcut around Rhie-Chow/SIMPLE/PISO development.


## Collocated FVM pressure/velocity coupling

Phase 4C adds `pressure_velocity` geometry/BC/flux primitives and `CollocatedIncompressible`. Momentum prediction retains pressure-free `HbyA` and cell `V/aP`; Rhie-Chow-style face flux uses direct pressure differences, while the orthogonal pressure contribution forms a matrix-free SPD operator and non-orthogonal flux is iterated explicitly. SIMPLE, PISO and PIMPLE-style controls share this same discretization rather than duplicating solver kernels.
