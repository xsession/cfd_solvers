# cfd_solvers v0.16.6 - distributed resident DEM contact/force pipeline

v0.16.6 turns the v0.16.5 MPI/SYCL transport seam into a complete resident distributed DEM timestep. The implementation keeps stable-ID Mindlin/contact state on the accelerator, overlaps halo traffic with interior owned-owned contact work, evaluates local-ghost Hertz-Mindlin interactions after halo completion, reverse-exchanges compact equal-and-opposite forces/torques, and then integrates the owned device state.

## Split resident timestep

`ResidentDemSycl` now exposes a split-step API:

- `begin_distributed_contact_step()` builds the owned neighbor buckets, clears forces/gravity and advances the contact-history stamp once;
- `apply_owned_contacts_interior_device()` evaluates pairs whose two particles are outside the halo region;
- `apply_owned_contacts_boundary_device()` evaluates the complementary owned-owned boundary partition;
- `apply_ghost_contacts_device()` evaluates local-ghost Hertz/Mindlin contacts under deterministic lower-global-ID ownership and emits compact remote force records;
- `apply_remote_forces_device()` accumulates reverse-exchanged force/torque records by stable global ID;
- `finish_distributed_contact_step()` cleans stale Mindlin history and integrates the resident particles.

The interior and boundary phases use one contact-history stamp, so splitting communication/computation does not reset tangential spring memory.

## Nonblocking MPI overlap

`MpiResidentDemSyclExchange::step()` performs:

1. owned-particle/history migration;
2. bond-history migration with the lower-ID particle owner;
3. device halo classification and packing;
4. nonblocking `MPI_Ialltoallv` halo exchange;
5. interior owned-owned contact work while the halo request is active;
6. halo completion and device ghost compaction;
7. boundary owned-owned and local-ghost contact work;
8. compact reverse-force binning by destination rank;
9. nonblocking reverse `MPI_Ialltoallv`;
10. destination force accumulation and resident integration.

Both existing transport modes remain supported. `staged_host` stages only compact migration/halo/reverse-force payloads. `direct_device` passes SYCL USM pointers to MPI and still requires explicit `assume_gpu_aware_mpi=true` opt-in.

## Cross-rank bonds

`ResidentDemPackedBond` stores canonical stable particle IDs, rest length, normal/shear stiffness and damping, tensile/shear strengths, damage onset, tangential displacement, damage and broken state. The bond owner is the rank owning the lower global particle ID. When that particle migrates, the complete bond state migrates with it. Cross-rank partner forces use the same compact reverse-force exchange as Hertz-Mindlin contacts.

## Validation

The release is validated by:

- **138/138** default CTest targets;
- Python ABI **0.16.6**;
- strict `-std=c++20 -Wall -Wextra -Wpedantic -Werror` compilation of the resident SYCL and MPI/SYCL sources;
- CPU-executing SYCL emulation of a two-rank contact ownership scenario, including normal/tangential momentum conservation and persistent Mindlin history;
- cross-rank resident bond force/history regression;
- one-rank MPI/SYCL execution of the new nonblocking orchestration path.

A physical multi-rank GPU-aware MPI environment is not available here. Therefore runtime qualification of direct-device MPI and the machine-counted distributed-memory DEM checkbox remain open even though the end-to-end code path is implemented.

## Tracker

No capability checkbox is closed in this release. Phase 13 remains **30/37 (81.1%)** and overall completion remains **643/756 (85.1%)**. The remaining gap for the distributed-memory DEM item is real multi-rank physical-runtime qualification rather than another missing resident contact/force stage.
