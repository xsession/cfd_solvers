# cfd_solvers v0.16.4 - distributed DEM contact and force completion

v0.16.4 turns the v0.16.3 DEM migration/ghost layer into a complete CPU-side distributed contact timestep architecture while preserving deterministic pair/history ownership.

## Deterministic cross-rank pair ownership

Every distributed sphere has a stable 64-bit global ID. For a cross-rank pair, only the rank owning the lower global ID evaluates the contact. That rule gives each pair exactly one numerical owner and also defines where its persistent Mindlin state lives.

`evaluate_distributed_dem_contacts()` now handles local-local and local-ghost interactions with:

- Hertz nonlinear normal contact;
- normal damping and cohesion;
- history-dependent Mindlin tangential displacement;
- Coulomb return limiting;
- rolling-resistance torque;
- equal-and-opposite local/remote force and torque contributions.

The higher-ID owner deliberately skips the mirrored local-ghost pair, preventing duplicate interface forces.

## Reverse force exchange

`MpiDemDomainExchange::exchange_forces()` groups remote reaction forces by owner rank and performs an `MPI_Alltoallv` reverse exchange. Incoming contributions are accumulated into the locally owned particle before integration, so cross-rank contacts preserve Newton's third law and global linear momentum up to integration/roundoff error.

## Persistent state migration

Contact-history and bond records use the same lower-global-ID ownership rule. During particle migration, records owned by a leaving lower-ID particle migrate with it. This prevents loss of Mindlin tangential history or bond damage state merely because a particle crosses a slab boundary.

`DistributedDemBond` provides the same stable-ID representation for bonded-particle state. Cross-rank bond evaluation returns remote action/reaction contributions and retains progressive damage/fracture history on the deterministic owner rank. The remote bonded partner must be present in the current ghost halo; callers should size `ghost_width` to cover the active bonded neighborhood.

## Complete MPI CPU timestep driver

`MpiDemDomainExchange::step()` now performs the complete baseline sequence:

1. migrate ownership;
2. migrate lower-ID contact/bond state;
3. rebuild adjacent ghost particles;
4. evaluate local/local and local/ghost Hertz-Mindlin contacts;
5. evaluate distributed bonds;
6. reverse-exchange remote force/torque reactions;
7. accumulate local and incoming contributions;
8. integrate locally owned particle linear/angular state.

Particles that cross a slab during integration are migrated at the start of the next timestep.

## Validation

A deterministic two-rank emulation validates:

- exactly one owner for a cross-rank contact;
- equal-and-opposite interface force;
- global momentum conservation after applying reverse reactions;
- persistent Mindlin history age across timesteps;
- cross-rank bonded-particle action/reaction;
- migration detection after integration;
- ghost-width contact-halo planning.

The real MPI regression suite is extended to run the same ownership/momentum/history-migration scenario whenever an MPI build is available. In this environment, MPI is not installed, so `src/distributed/mpi_dem.cpp` is strict-compiled against the existing MPI compile harness rather than executed with multiple ranks.

The complete default project matrix passes **136/136 CTest targets**, including Python ABI **0.16.4**.

## Tracker

No new capability checkbox is claimed in v0.16.4. Phase 13 remains **30/37 (81.1%)** and overall completion remains **643/756 (85.1%)**. The distributed-memory DEM checkbox stays open until a real multi-rank runtime executes the end-to-end MPI test successfully. This is deliberate qualification discipline rather than an implementation gap in the CPU timestep path.
