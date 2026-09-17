# cfd_solvers v0.16.3 - bonded/fracturing DEM and distributed foundations

v0.16.3 matures Phase 13 in two directions: a history-bearing bonded-particle constitutive law for fracture-capable DEM, and an ownership/migration/ghost-exchange layer for distributed-memory DEM.

## Bonded and fracturing DEM

`ParticleBond` connects two DEM bodies with a rest length and a per-bond `BondedParticleModel`. Each update tracks tangential displacement in the current bond tangent plane and evaluates normal and shear elastic forces plus viscous damping.

Failure is history-bearing rather than a one-shot distance cutoff:

1. normal extension and tangential spring displacement produce elastic force demand;
2. tensile and shear demand are normalized by independent failure-force thresholds;
3. damage begins at the configurable `damage_onset_ratio`;
4. damage grows monotonically and reduces elastic stiffness;
5. once either normalized demand reaches unity, the bond breaks irreversibly and stops transmitting force.

`ExplicitDemSystem::add_bond()` integrates the model into the ordinary DEM step. `last_bond_stats()` reports active, damaged and newly/previously broken bonds for diagnostics and fracture monitoring.

This is a force-threshold baseline; it intentionally does not assume a hidden bond area. A later cohesive-zone/beam-bond model can add stress-based calibration, bending/torsional bond moments and fracture-energy regularization.

## Distributed DEM foundations

`DemSlabDecomposition` partitions a global x extent into equal-width ownership slabs. `plan_dem_slab_exchange()` deterministically identifies:

- particles whose centers migrated to another owner rank;
- left/right adjacent ranks requiring a ghost copy because the expanded sphere AABB overlaps the slab boundary.

`DistributedDemParticle` preserves a stable 64-bit global ID together with rigid-body state, mass/inertia and radius.

When `CFD_ENABLE_MPI=ON`, `MpiDemDomainExchange` performs:

1. ownership reclassification,
2. packed particle migration with `MPI_Alltoallv`,
3. owner-rank normalization of incoming particles,
4. a second planning pass after migration,
5. adjacent-rank ghost exchange with `MPI_Alltoallv`.

The MPI test suite now contains a migration/ghost smoke path. The current execution environment does not provide a real MPI runtime, so this source was strict-compiled with a minimal MPI header harness but was not executed across ranks here. The tracker therefore keeps the full distributed-memory DEM checkbox open until real multi-rank contact/integration validation is available.

## Validation

The v0.16.3 regression covers:

- elastic bonded-particle action/reaction;
- progressive damage before failure;
- irreversible tensile fracture;
- persistent tangential/shear history;
- bond integration inside `ExplicitDemSystem`;
- deterministic slab owner selection;
- migration destination planning;
- left/right ghost planning;
- distributed-particle state conversion.

Modified multibody/distributed code and the MPI exchange source pass strict C++20 syntax compilation with `-Wall -Wextra -Wpedantic -Werror`. The complete default project matrix passes **135/135 CTest targets**, including Python ABI **0.16.3**.

## Tracker

The `bonded particles and fracture/bond failure` Phase-13 item is now complete. Phase 13 advances from **29/37 (78.4%)** to **30/37 (81.1%)**. Overall tracked completion advances from **642/756 (84.9%)** to **643/756 (85.1%)**.

The distributed-memory item remains open because multi-rank runtime/contact-loop validation is still pending despite the ownership, migration, ghost-planning and MPI exchange implementation being present.
