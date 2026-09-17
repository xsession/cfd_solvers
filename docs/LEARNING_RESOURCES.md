# Learning resources for the implemented solver stack

This file collects learning resources for the domains implemented in this repository. Prefer official documentation, open-source project docs and primary method papers when extending the code. Do not copy incompatible source code into this MIT repository; use these resources for architecture, terminology, validation ideas and independent clean-room implementation.

## How to use these resources

- Start with `docs/IMPLEMENTED_SOLVER_STUDY_GUIDE.md` for the local implementation map.
- Use the links below to understand standard terminology and validation expectations.
- When implementing a feature inspired by an external project, record the source in `docs/upstreams.json` or in the relevant phase doc.
- Keep formulas and algorithms in your own words unless a license explicitly allows reuse.

## Electromagnetic PIC and plasma

- WarpX theory: explains the core PIC loop: particle push, deposition, field solve/update and field gather. Useful for mapping the implemented `ElectromagneticPic1D`, `ElectromagneticPic2D` and `StaggeredElectromagneticPic2D` loops to standard PIC terminology. <https://warpx.readthedocs.io/en/22.11/theory/picsar_theory.html>
- WarpX parameters: useful for seeing production options such as grid type, current deposition, field gather and implicit/semi-implicit controls. <https://warpx.readthedocs.io/en/24.08/usage/parameters.html>
- Smilei overview and algorithms: useful for understanding electromagnetic PIC code organization, macro-particle weighting, diagnostics and plasma use cases. <https://smileipic.github.io/Smilei/> and <https://smileipic.github.io/Smilei/Understand/algorithms.html>
- Smilei implementation notes: useful for high-level data-container and time-loop organization ideas. <https://smileipic.github.io/Smilei/implementation.html>
- PIConGPU PIC algorithm and current deposition docs: useful for charge-conserving deposition concepts such as Esirkepov/EZ and particle-shape options. <https://picongpu.readthedocs.io/en/latest/models/pic.html> and <https://picongpu.readthedocs.io/en/latest/usage/param/particles/current.html>
- PIConGPU PICMI intro: useful for seeing the standard Yee electromagnetic solver exposed through a high-level simulation interface. <https://picongpu.readthedocs.io/en/latest/usage/picmi/intro.html>

## FDTD and full-wave EM

- openEMS documentation: open-source EC-FDTD solver docs for simulation setup, ports, materials and FDTD workflow expectations. <https://docs.openems.de/en/latest/>
- openEMS introduction: concise overview of openEMS as a free/open electromagnetic field solver based on FDTD/EC-FDTD. <https://docs.openems.de/en/latest/intro.html>
- MFEM Maxwell notes: useful for understanding Nedelec H(curl) finite elements, tangential continuity and Maxwell discretization terminology. <https://mfem.org/maxwell-notes/>
- MFEM tutorials and examples: useful for H(curl) finite-element examples and Maxwell miniapp context. <https://mfem.org/tutorial/further/> and <https://mfem.org/examples/>

## SPICE, MNA and compact models

- ngspice documentation: primary open-source SPICE manual and control-flow documentation. Use it to check netlist syntax, analysis names and simulator workflow concepts. <https://ngspice.sourceforge.io/docs.html>
- ngspice beginner tutorial: useful for basic SPICE simulation terminology and passive/active device examples. <https://ngspice.sourceforge.io/ngspice-tutorial.html>
- Xyce documentation and tutorials: useful for large-scale SPICE-compatible solver architecture, parallel circuit simulation and reference-guide organization. <https://xyce.sandia.gov/documentation-tutorials/>
- Xyce home page: describes Xyce as an open-source, SPICE-compatible analog circuit simulator for large-scale parallel platforms. <https://xyce.sandia.gov/>

## CFD, FVM and LBM

- OpenFOAM user guide: use for case organization, discretization scheme vocabulary, solver controls and finite-volume workflow. <https://www.openfoam.com/documentation/user-guide>
- OpenFOAM standard solvers reference: useful for mapping solver-family coverage and naming conventions. <https://www.openfoam.com/documentation/user-guide/a-reference/a.1-standard-solvers>
- OpenFOAM schemes documentation: useful for conservative scalar-transport form and discretization-scheme terminology. <https://doc.openfoam.com/2606/tools/processing/numerics/schemes/>
- FluidX3D repository/documentation: useful as a performance and LBM architecture reference, especially OpenCL portability and memory-efficient DNS/LBM organization. <https://github.com/ProjectPhysX/FluidX3D>

## FEM and multiphysics

- Elmer FEM repository: useful as a broad open-source multiphysics FEM reference covering electromagnetics, structural mechanics, fluid mechanics and parallel FEM organization. <https://github.com/ElmerCSC/elmerfem>
- Elmer documentation archive and models manual: useful for solver-family organization and multiphysics model documentation patterns. <https://sourceforge.net/projects/elmerfem/files/ElmerDocumentation/> and <https://ftp.fi.netbsd.org/pub/csc/software/elmer/doc/ElmerModelsManual.pdf>
- MFEM getting started and feature pages: useful for modern finite-element software structure, examples and scalable solver terminology. <https://mfem.org/getting-started/> and <https://mfem.org/features/>

## RF, antennas, PEEC, EMC and SI/PI

- openEMS examples and docs are the closest open reference in this list for ports, S-parameters, FDTD antennas and RF workflows. <https://docs.openems.de/en/latest/>
- MFEM Maxwell notes and examples are useful for RF FEM and waveguide/cavity discretization checks. <https://mfem.org/maxwell-notes/>
- ngspice/Xyce resources are useful for circuit-level macromodel validation of RF network reductions and sampled transmission-line elements. <https://ngspice.sourceforge.io/docs.html> and <https://xyce.sandia.gov/documentation-tutorials/>

## Bioelectromagnetics and thermal coupling

- Start locally with `include/cfd/multiphysics/bioheat.hpp` and `src/multiphysics/bioheat.cpp`.
- For external context, use CST/SIMULIA-style domain descriptions captured in `docs/CST_COVERAGE_AUDIT.md`; they identify SAR, RF exposure and EM-to-thermal workflows as key coverage targets.
- Use general heat-transfer/FEM resources from Elmer and MFEM to validate thermal discretization patterns.

## Clean-room implementation checklist

Before adding a feature based on these resources:

1. Record the external resource and what it influenced.
2. Write a new public API in `include/cfd/...`.
3. Implement a minimal validated baseline.
4. Add a focused regression test.
5. Add a CLI smoke case if it is user-visible.
6. Update `docs/INTEGRATION_TRACKER.md` only after the test passes.
7. Update this learning-resource file if the feature introduces a new domain.


## 3-D PIC study checklist

Use `ElectrostaticPic3D` and `ElectromagneticPic3D` as stepping stones before production staggered 3-D EM-PIC:

1. Confirm trilinear CIC deposition conserves total macro-charge.
2. Verify periodic Poisson by using a manufactured sinusoidal potential.
3. Reconstruct spectral current from two charge snapshots and check the continuity residual.
4. Study `ElectromagneticPic3D` as the first full-field 3-D/3-V reference: full E/B storage, centered curls, Vay push and spectral Jx/Jy/Jz coupling.
5. Only then replace the nodal reference with a staggered 3-D Yee field layout and local charge-conserving deposition.

Local files to inspect:

- `include/cfd/particle/electromagnetic.hpp`
- `src/particle/electromagnetic.cpp`
- `tests/test_v097_pic3d.cpp`
- `tests/test_v098_em_pic3d.cpp`
- CLI: `./build/cfd-solve particle-pic3d`
- CLI: `./build/cfd-solve particle-em-pic3d`


## v0.9.9 learning checklist - 3-D Yee EM-PIC

Focus on the leapfrog structure before optimizing:
- Map electric components to edges and magnetic components to faces.
- Write the six curl equations on a Yee cell and compare them with `StaggeredElectromagneticPic3D::advance_magnetic` and `update_electric`.
- Re-run `cfd-v099-staggered-pic3d-tests`, then change the CFL factor to see where the guard rejects the step.
- Use the 3-D sponge/electric-wall primitives as a stepping stone toward CPML and conductor boundaries.
- Next reading target: charge-conserving current deposition such as Villasenor-Buneman/Esirkepov, because v0.9.9 still uses the validated spectral continuity reconstruction.

## v0.10.0 learning checklist - local current deposition

Focus for this release:

- Derive the discrete continuity equation on a finite-volume Yee grid.
- Compare global spectral current reconstruction against local cumulative finite-volume reconstruction.
- Study Villasenor-Buneman and Esirkepov charge-conserving current deposition as the next production-grade step.
- Review how absorbing sponge layers differ from CPML; v0.10.0 only adds a configurable sponge profile, not a full CPML.

Practical reading targets:

- WarpX theory notes on PIC and charge-conserving current deposition.
- Smilei documentation sections on interpolation/projection.
- PIConGPU documentation on particle shapes and current deposition.
- openEMS FDTD/absorbing-boundary documentation as a bridge from simple sponge damping to CPML/PML.


## v0.10.1 learning checklist - particle sorting and guard halos

Focus for this release:

- Study why production PIC codes periodically sort particles by cell or tile for cache/GPU locality.
- Understand the difference between particle sorting, field guard cells and particle halo exchange.
- Review how edge/corner particles require multiple halo destinations in Cartesian domain decomposition.
- Keep current deposition and sorting order separate: deposit current from matched old/new particle arrays, then reorder for the next step.

Practical reading targets:

- WarpX and PIConGPU documentation on particle containers/tiles and guard cells.
- Smilei documentation on patch decomposition and particle exchange.
- General MPI halo-exchange tutorials for Cartesian grids before mapping this to particle lists.


## v0.10.2 learning checklist - PIC decomposition and guard exchange

Study order:

1. Review why PIC codes split space into subdomains before splitting particles: field updates need neighbouring guard cells, while particles migrate between owning domains.
2. Study Cartesian block decomposition: global cell ranges, physical extents, uneven remainder handling and stable domain IDs.
3. Study particle migration separately from field halo exchange. The v0.10.2 code intentionally separates `plan_particle_domain_migration_3d` from `exchange_scalar_guard_cells_3d`.
4. Practice ghost-cell indexing with a scalar field before applying it to staggered EM fields.
5. Only after this serial pack/scatter contract is clear, move to MPI send/receive, guard-cell exchange and load balancing.

Local exercises:

- Run `./build/cfd-v0102-domain-exchange-tests`.
- Run `./build/cfd-solve particle-domain-exchange3d`.
- Inspect `ScalarGuardedBlock3D::padded_nx/padded_ny/padded_nz` and draw the interior plus one-cell ghost layer on paper.
- Compare the v0.10.1 guard-halo particle classification with the v0.10.2 full domain migration buckets.


## v0.10.3 learning checklist - communicator-ready PIC exchange

Focus for this release:

- Understand the difference between a **decomposition contract** and an actual communication backend. v0.10.3 creates explicit messages but does not call MPI.
- Study halo-exchange message anatomy: source rank/domain, destination rank/domain, payload coordinates and payload values.
- Study particle migration anatomy: retained particles, outbound messages, destination ownership and outside-particle diagnostics.
- Compare `exchange_scalar_guard_cells_3d` with `pack_scalar_guard_cell_messages_3d` + `apply_scalar_guard_cell_messages_3d`; the latter is what MPI should wrap.

Practical exercises:

- Run `./build/cfd-v0103-comm-exchange-tests`.
- Run `./build/cfd-solve particle-comm-exchange3d`.
- Print the guard messages for a 2x2x1 domain grid and classify face, edge and corner messages by `offset_x/y/z`.
- Build a small mock communicator that serializes each `ScalarGuardCellMessage3D` to bytes, deserializes it, and applies it.

Reading targets:

- MPI Cartesian topology and neighborhood collective examples.
- General structured-grid halo-exchange tutorials.
- WarpX/Smilei/PIConGPU documentation sections on guard cells, patch exchange and particle redistribution.


## v0.10.4 learning checklist - from serial PIC exchange to MPI

- Understand why a solver should first define deterministic message payloads before adding MPI calls. This keeps numerical ownership, communication, and transport backend separate.
- Trace the data path: particles/cell values -> migration/guard messages -> rank envelopes -> transport inboxes -> reconstructed domain buckets/guard blocks.
- Review MPI concepts needed for the next step: rank, communicator, nonblocking send/receive, message tags, derived datatypes vs byte payloads, and deadlock-free exchange ordering.
- Practical exercise: use `particle-transport3d` to compare same-rank and remote-rank message counts for different logical rank topologies.


## v0.10.5 addendum - learning distributed PIC transport

Study order for the new transport layer:

1. Domain decomposition: understand cells owned by each rank/subdomain.
2. Particle migration: particles leaving a subdomain become messages to the owning subdomain.
3. Guard exchange: field/charge ghost cells are copied from neighbours before stencil operations.
4. Serialization: MPI/GPU transports should move opaque byte records rather than depend on solver internals.
5. Validation: compare serialized transport with the serial apply path before enabling real MPI.

Practical exercises now included in this repository:

- `cfd-v0105-serialized-transport-tests` validates byte roundtrip and malformed input rejection.
- `cfd-solve particle-serialized-transport3d` prints message counts, serialized bytes, particle count and guard payload count.


## v0.10.5 addendum - plasma chemistry and RF breakdown study path

To understand the new plasma/breakdown helpers, study in this order:

1. Number-density reaction networks and stoichiometric source terms.
2. Charge conservation in ionization/recombination reactions.
3. Townsend avalanche and Paschen-law gas breakdown.
4. Secondary electron yield and multipactor resonance order.
5. How these threshold tools should later consume CST-like RF field maps from FEM/FDTD/MoM solvers.

Repository exercises:

- `cfd-v0105-plasma-tests` checks charge preservation for an electron-neutral ionization reaction.
- `cfd-solve particle-breakdown-threshold` reports Paschen and multipactor threshold estimates.


## v0.10.6 addendum - studying distributed PIC exchange rounds

Study the production distributed-PIC path in layers:

1. Particle ownership: map particles to the domain that owns their position.
2. Field ownership: split field arrays into domain-local blocks.
3. Guard cells: exchange neighbouring field values before stencil/curl operations.
4. Transport envelopes: add source/destination rank metadata without changing solver payloads.
5. Serialization: move opaque byte records through MPI/GPU/runtime layers.
6. Exchange round: pack -> serialize -> deliver -> deserialize -> apply -> compare with direct reference exchange.

Repository exercises:

- `cfd-v0106-distributed-round-tests` validates the exchange round without MPI.
- `cfd-solve particle-distributed-round3d` prints rank, message, byte and particle counts.
- `cfd-solve particle-field-guards3d` checks six-component EM field guard exchange.

Suggested next reading outside the repository: MPI neighborhood collectives, guard/halo exchange patterns in structured-grid solvers, and particle-in-cell domain-decomposition algorithms before attempting GPU kernels.

## v0.10.7 addendum - studying distributed PIC timestep orchestration

Focus topics for this step:

- Domain-local particle push before ownership migration.
- Guard-cell reconstruction before the next field update.
- Separation between numerical physics kernels and communicator/message contracts.
- Why deterministic serial transports are useful before MPI/GPU execution.

Study path:

1. Review the existing particle migration and guard-cell exchange sections.
2. Read `tests/test_v0107_distributed_step.cpp`.
3. Run `cfd-solve particle-distributed-step3d`.
4. Sketch the future MPI version: local push -> pack outbound particles/guards -> send/receive -> apply inbound payloads -> update next local fields.

## Geant4-style particle-through-matter transport

Study checklist:

1. Understand the track/step split: a track is the persistent particle state; a
   step is one proposed and accepted movement through geometry/material.
2. Understand GPIL-style step arbitration: user step, geometry boundary,
   continuous/discrete process limits and energy cuts compete for the next step.
3. Implement continuous processes as along-step energy/momentum changes and
   discrete processes as post-step interactions that may produce secondaries.
4. Keep scoring separate from physics: scorers consume completed steps and
   secondaries instead of owning the stepping logic.
5. Treat Geant4 as a reference architecture. Do not copy Geant4 source; use
   published physics formulas or independent derivations for process models.

Repository exercise in this tree:

```bash
./cfd-solve particle-geant4-transport
ctest -R cfd-v0108-geant4-transport-tests --output-on-failure
```


## v0.10.9 study queue - Geant4, SU2 and csauto clean-room lessons

### Geant4-inspired transport deepening
- Study physical-interaction-length sampling and how stochastic process arbitration differs from deterministic max-step limiting.
- Extend the current axis-aligned BVH to nested volumes, repeated placements and material overrides.
- Add hit collection exporters and dose validation cases before treating the bioheat coupling as quantitative dosimetry.

### SU2-inspired PDE and optimization workflow
- Keep solver kernels, iteration orchestration, numerics, output and workflow scripts separated.
- Add objective-function objects and finite-difference/adjoint-compatible gradient interfaces before tying optimization directly into any solver.
- Add regression cases for gradients, not just primal solution values.

### csauto-inspired campaign automation
- Keep solver-specific runtime knowledge behind a small adapter descriptor and test that generic workflow code does not depend on solver names.
- Treat DOE generation, template rendering, doctor checks, run registry, status, residual/performance extraction and cleanup as separate concerns.
- Add persistent campaign folders and JSON/CSV exports only after the in-memory campaign primitives are stable.
