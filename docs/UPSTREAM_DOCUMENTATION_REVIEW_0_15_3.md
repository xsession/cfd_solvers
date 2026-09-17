# Upstream documentation review for v0.15.3 FVM topology change and AMR

This review records public mesh-adaptation concepts used to guide an independent C++20 implementation. No OpenFOAM source code was copied or mechanically translated. `cfd_solvers` keeps its own cell metadata, topology-change API, overlap mapping and `PolyMesh` reconstruction.

## Public references reviewed

### Dynamic finite-volume meshes and topological change

OpenFOAM's user-guide library reference identifies `dynamicFvMesh` as the finite-volume mesh layer that can move and undergo topological changes:

- https://www.openfoam.com/documentation/user-guide/a-reference/a.3-standard-libraries

Clean-room mapping: `MeshTopologyOperation` is a small project-native operation boundary. Each operation returns a new `AdaptiveHexMesh` and a geometric old/new overlap map instead of mutating OpenFOAM mesh primitives.

### Hex refinement

The OpenFOAM utility reference documents `refineHexMesh` as 2x2x2 splitting of selected hexahedral cells and separately lists `refineMesh` for directional cell refinement:

- https://www.openfoam.com/documentation/user-guide/a-reference/a.2-standard-utilities

OpenFOAM 13 also describes `hexRef8` as a refinement mode that supports further refinement:

- https://openfoam.org/release/13/

Clean-room mapping: `MarkedHexRefinement` independently splits an axis-aligned leaf cell into eight children. A compact root/level/lineage identity supports repeated refinement and complete-sibling coarsening without requiring point/edge storage in the current `PolyMesh` representation.

### Coarse/fine faces and field mapping

OpenFOAM's v1812 numerics notes describe mapping improvements for the new internal faces produced by 2x2x2 cell splitting, emphasizing that topology change requires deliberate field initialization/mapping rather than treating created faces as ordinary copies:

- https://www.openfoam.com/news/main-news/openfoam-v1812/numerics

The general `mapFields` documentation distinguishes source and target geometries and treats field transfer as an explicit operation:

- https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.6-mapping-fields-between-different-geometries

Clean-room mapping: v0.15.3 computes exact axis-aligned cell intersection volumes. New cell averages are volume-weighted from every overlapping old cell, making cell-integral conservation a direct property of the map. Coarse/fine finite-volume connectivity is rebuilt as multiple subfaces whose areas exactly tile the shared interface.

### Refinement-level transitions

OpenFOAM's v2406 preprocessing notes discuss `nCellsBetweenLevels` and the handling of slower-than-2:1 refinement transitions in `snappyHexMesh`:

- https://www.openfoam.com/news/main-news/openfoam-v2406/pre-processing

Clean-room mapping: the v0.15.3 baseline exposes a simpler optional 2:1 face-neighbour constraint. After marked refinement or coarsening, coarse neighbours are recursively refined when adjacent levels differ by more than one.

## Numerical/architectural conclusions carried into cfd_solvers

1. Mesh topology change is a separate operation boundary from mesh motion; an ALE coordinate update alone cannot represent cell creation/deletion.
2. Adaptation must return mapping metadata together with the new mesh so solver fields can be transferred deterministically.
3. For cell-average conservative quantities, overlap-volume mapping is preferable to nearest-centre copying plus global rescaling because it conserves locally over every target control volume.
4. Coarse/fine finite-volume interfaces need geometrically partitioned subfaces so flux area is not duplicated or lost.
5. Refinement and coarsening need lineage/sibling identity; geometric coincidence alone is not a robust coarsening criterion.
6. A residual/error indicator and marking policy should remain separate from topology mechanics so future flow-specific indicators can reuse the same mesh operations.

## Intentional scope boundaries

- The adaptive mesh is currently axis-aligned Cartesian hex/octree topology, not arbitrary-polyhedron cutting/remeshing.
- Refinement is isotropic 2x2x2; directional anisotropic split patterns are not yet implemented.
- The first estimator is a scalar internal-face jump indicator. Goal-oriented, Hessian/metric and flow-feature indicators remain future work.
- Parallel/distributed AMR and adaptive repartition are not claimed here.
- Face-centered/vector/tensor conservative remapping is not yet generalized; this release establishes the cell-average scalar path and topology map needed for those extensions.
