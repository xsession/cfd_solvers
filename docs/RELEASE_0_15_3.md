# Release 0.15.3 - conservative FVM topology change and adaptive hex refinement

## Highlights

- Adds `AdaptiveHexMesh`, an axis-aligned leaf-cell representation that can be converted into the existing owner/neighbour `PolyMesh` format.
- Adds the polymorphic `MeshTopologyOperation` interface plus marked 2x2x2 hex refinement and complete-sibling coarsening operations.
- Adds optional 2:1 face-neighbour balancing after topology changes.
- Reconstructs coarse/fine interfaces as multiple finite-volume subfaces, preserving interface area and orientation.
- Adds `TopologyChangeMap` with exact source/target overlap volumes.
- Adds conservative overlap-volume remapping for scalar cell averages.
- Adds scalar jump-error indicators, Dörfler bulk marking and an automatic coarsen/refine AMR cycle.
- Closes all three remaining Phase-3 mesh-motion/adaptation tracker items.

## Topology-change boundary

`include/cfd/fvm/adaptive_mesh.hpp` introduces `MeshTopologyOperation`. An operation consumes an immutable `AdaptiveHexMesh` and returns `AdaptiveTopologyResult`, which contains both the replacement mesh and a `TopologyChangeMap`.

This keeps topology mechanics independent of any particular flow solver. Pressure, turbulence, VOF, species and thermal solvers can reuse the same operation/mapping boundary as their AMR integration is added.

## Local 2x2x2 refinement and coarsening

`MarkedHexRefinement` splits each marked leaf into eight equal-volume children. Every leaf retains its root ID, refinement level and compact octree lineage. `MarkedHexCoarsening` only merges complete groups of eight marked siblings, so unrelated cells cannot be combined accidentally.

When 2:1 balancing is enabled, face-adjacent leaves are checked after the requested topology operation and coarse cells are recursively split until neighbouring levels differ by at most one.

## PolyMesh coarse/fine connectivity

`AdaptiveHexMesh::poly_mesh()` directly reconstructs the cell-centered owner/neighbour representation. For a coarse cell adjacent to four fine cells on one face, four internal `Face` records tile the common plane. Their summed oriented area equals the original coarse face area, so downstream finite-volume flux assembly sees the correct control-surface measure without adding hanging-node concepts to `PolyMesh` itself.

The domain boundary keeps the existing six patches: `left`, `right`, `bottom`, `top`, `front`, `back`.

## Conservative field transfer

`build_topology_change_map(...)` computes all non-zero old/new cell intersection volumes. `conservative_adaptive_remap(...)` forms each new cell average from these overlap weights.

For a cell-average field `phi`, each target value is

`phi_new[j] = sum_i(phi_old[i] * overlap(i,j)) / V_new[j]`.

Because the adaptive cells partition the same domain, the integral `sum(phi*V)` is preserved to floating-point roundoff through both refinement and coarsening.

## Error indicator and automatic AMR

`scalar_jump_error_indicator(...)` accumulates squared internal-face jumps weighted by face area and inverse owner-neighbour distance. `mark_dorfler(...)` then selects the smallest deterministic high-error prefix whose squared indicators reach a requested fraction of total estimator energy.

`adapt_scalar_field(...)` combines:

1. low-error complete-sibling coarsening;
2. conservative remapping;
3. recomputed jump indicators;
4. Dörfler refinement up to `maximum_level`;
5. optional 2:1 balancing and a final conservative remap.

The marking policy is deliberately independent from topology change, leaving room for pressure, vorticity, interface, shock and adjoint indicators later.

## Validation

Focused target:

- `cfd-v0153-fvm-amr-tests`

The regression verifies:

- topology operations through the abstract interface;
- one-cell 2x2x2 splitting and total-volume conservation;
- four conservative subfaces at a level-0/level-1 interface;
- exact scalar-integral conservation during refinement;
- complete-sibling coarsening and refine/coarsen field round-trip;
- recursive 2:1 balancing;
- scalar-jump detection and Dörfler marking;
- automatic refinement at a discontinuity;
- low-error coarsening of a uniform field.

The complete default CPU regression matrix passes **123/123** tests.

See `UPSTREAM_DOCUMENTATION_REVIEW_0_15_3.md` for the public documentation review and clean-room design mapping.
