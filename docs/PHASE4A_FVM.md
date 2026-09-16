# Phase 4A - polyhedral FVM foundation and pressure projection

Phase 4A begins the OpenFOAM-class finite-volume track while keeping the implementation clean-room C++20.

## Polyhedral mesh model

`cfd::fvm::PolyMesh` is topology-first rather than Cartesian-grid-first. It stores:

- cell center and volume;
- face owner and optional neighbour cell;
- face center;
- oriented face-area vector `Sf`;
- boundary-patch index for boundary faces;
- cell-to-face adjacency built and validated once.

For an internal face, `Sf` points from owner to neighbour. For a boundary face it points out of the owner cell. This gives every finite-volume operator one consistent sign convention.

`make_cartesian_hexa_mesh()` is currently the built-in mesh generator and is primarily a validation/bring-up path. The `PolyMesh` container itself accepts arbitrary cell/face geometry, but file import, non-convex validation and mesh-quality metrics are later milestones.

## Cell-parallel Gauss operators

Phase 4A adds:

- scalar Gauss gradient;
- vector Gauss divergence;
- orthogonal two-point scalar Laplacian.

The implementation loops over each cell's face adjacency rather than atomically accumulating through faces. This preserves readable owner/neighbour algebra and allows OpenMP parallelization without races.

On the generated orthogonal hexahedral mesh the regression suite verifies:

- exact linear scalar gradient;
- exact divergence of a linear vector field;
- exact interior Laplacian of a quadratic scalar field to floating-point tolerance.

Boundary values are supplied as face-indexed arrays. An empty boundary array selects zero-normal-gradient behavior for the current low-level operators.

### Current discretization limitation

The Laplacian is intentionally an orthogonal two-point flux approximation. Non-orthogonal correction, skewness correction, bounded interpolation and convection limiters are Phase 4B work. Phase 4A does not pretend a Cartesian verification result is sufficient for general skew polyhedra.

## Shared matrix-free conjugate gradient

`cfd/core/conjugate_gradient.hpp` adds a reusable matrix-free CG implementation with aligned work vectors and OpenMP reductions/updates. It is already used by the pressure projection and is intended to be shared with later FVM and FEM systems.

The direct regression solves a diagonal SPD system independently of the CFD path.

## Staggered incompressible pressure projection

`Projection2D` is a periodic MAC/staggered-grid finite-volume pressure-projection kernel:

- pressure is cell-centered;
- `u` is stored on x-normal faces;
- `v` is stored on y-normal faces;
- divergence is the conservative face-flux difference per cell;
- the pressure equation is solved matrix-free with CG;
- face velocities are corrected with pressure differences using the same staggered geometry.

A deliberately divergent periodic velocity field is projected to test the incompressibility kernel. On the 96x80 CLI case in this development environment:

- divergence RMS before projection: about `3.14091e-1`;
- divergence RMS after projection: about `2.85e-11`;
- pressure solve: 7 CG iterations with the default tolerance in that case.

Those values are regression observations on this case, not general performance guarantees.

Taylor-Green initialization is also projected. Its small discrete sampling divergence is driven to near machine precision while kinetic energy changes by less than the regression tolerance.

## What Phase 4A is not

This is not yet a full Navier-Stokes solver. It does not yet include:

- momentum convection;
- viscous momentum time integration on `PolyMesh`;
- collocated Rhie-Chow interpolation;
- SIMPLE, PISO or PIMPLE loops;
- non-orthogonal pressure correction;
- turbulence models;
- dynamic mesh or AMR;
- multiphase/compressible/reacting physics.

Those are built on top of the tested mesh, operator, flux and linear-solver foundation rather than being mixed into the first implementation step.
