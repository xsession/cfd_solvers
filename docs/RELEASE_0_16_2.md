# cfd_solvers v0.16.2 - resident SYCL DEM

v0.16.2 applies the device-residency strategy used by the LBM/FDTD/FVM accelerator paths to the Phase-13 granular solver. The CPU `ExplicitDemSystem` remains unchanged as the portable reference implementation; the new SYCL backend is a separate resident spherical-DEM execution path.

## Resident DEM architecture

`ResidentDemSycl` stores particle position, translational/angular velocity, radius, inverse mass, spherical inverse inertia, force and torque as device-owned SoA arrays. The caller supplies a SYCL queue; as with the other resident backends, the queue is expected to be in-order so hot-loop kernels remain sequenced without host waits.

Each timestep performs:

1. device cell counting,
2. a device prefix/bucket-offset pass,
3. device particle bucket compaction,
4. 27-cell neighbor traversal,
5. one pair evaluation for `j > i`,
6. persistent Mindlin history lookup/update,
7. Hertz normal + Mindlin tangential + rolling-resistance force/torque accumulation,
8. stale contact-history cleanup,
9. semi-implicit particle integration.

Forces and torques use relaxed device atomics only where multiple contacts accumulate on one particle. Contact history avoids a global concurrent hash table: each lower-index particle owns a bounded neighbor-history table, so its tangential displacement is updated deterministically by one work-item. When the configured history table is exhausted, excess contacts retain normal/damping/friction response but do not retain a tangential spring across timesteps.

## Cell-linked neighbor search

The GPU grid uses a uniform cell size at least as large as the largest particle diameter. A device atomic count plus single-device-work-item prefix pass produces compact per-cell ranges in `sorted_indices`. No particle-position array is downloaded or host-sorted during a timestep.

A CPU `broad_phase_cell_linked()` reference implementation was added to validate the same 27-neighbor-cell topology against the existing sweep-and-prune broad phase.

## Residency contract

After `upload()`, `step()` performs no explicit host-field transfer. `download()` is an output/API boundary. `total_kinetic_energy()` performs a device reduction and synchronizes only the scalar result. Transfer accounting is exposed through `DeviceTransferStats`.

## Current limitations

- spherical particles only on the resident GPU path;
- particle diameter must not exceed cell size;
- fixed bounded per-particle contact-history capacity;
- no GPU plane/triangle/convex collider yet;
- no bonded/fracturing DEM yet;
- no distributed DEM halo exchange yet;
- single precision storage/contact kernels in this baseline;
- physical NVIDIA/AMD/Intel SYCL runtime/performance qualification remains open in this environment.

## Validation

The v0.16.2 regression checks the CPU cell-linked candidate set against sweep-and-prune. SYCL-enabled builds additionally exercise resident overlapping-particle contact and verify that the timestep hot loop records zero host field-transfer bytes. The source/test also passes the strict fake-SYCL C++20 syntax harness with `-Wall -Wextra -Wpedantic -Werror`. The complete default project matrix passes **134/134 CTest targets**, including Python ABI **0.16.2** and the earlier Phase-13 regressions.

Phase 13 advances from **28/37 (75.7%)** to **29/37 (78.4%)**. Overall tracked completion advances from **641/756 (84.8%)** to **642/756 (84.9%)**.
