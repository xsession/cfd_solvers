# Implemented solver study guide

This guide is a learning map for the solver families currently implemented in `cfd_solvers`. It is not a theory textbook. It tells you what to study, which local source files to inspect, which regression tests prove the baseline, and which CLI case gives a quick runnable example.

## How to study this repository

Use this loop for every solver family:

1. Read the short domain section below.
2. Open the listed public header first; treat it as the API contract.
3. Open the implementation file and identify the state variables, update equation and boundary handling.
4. Run the listed focused test and smoke case.
5. Change one physical parameter and predict how the test output should move.
6. Only then add a new feature.

Useful commands:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCFD_ENABLE_OPENMP=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/cfd-solve --list
./scripts/integration_status.py
```

## 1. LBM fluid core

**What is implemented**

- D2Q9, D3Q19 and D3Q27 lattice descriptors.
- BGK/TRT-style collision baselines.
- Esoteric-pull single-grid streaming.
- Two-grid reference solvers.
- Guo forcing, bounce-back walls, Zou-He inlet/outlet, Poiseuille/channel/cavity/Taylor-Green checks.
- MPI halo decomposition and optional SYCL paths.

**Study files**

- `include/cfd/solvers/lbm/descriptors.hpp`
- `include/cfd/solvers/lbm/one_step_pull.hpp`
- `src/lbm/esoteric_pull.cpp`
- `src/lbm/advanced_d2q9.cpp`
- `src/distributed/mpi_runtime.cpp`
- `docs/ESOTERIC_PULL.md`
- `docs/LBM_BOUNDARIES.md`

**Run**

```bash
./build/cfd-solve lbm-trt
./build/cfd-bench --lattice d3q19 --nx 64 --ny 64 --nz 64 --steps 100 --streaming both
```

**Key ideas to learn**

Population distribution functions, collision/streaming split, lattice units, relaxation time, Mach-number limits, bounce-back, halo exchange, memory bandwidth and single-grid streaming hazards.

## 2. Finite-volume CFD / OpenFOAM-like layer

**What is implemented**

- Polyhedral owner/neighbour mesh topology.
- Gauss gradients, divergence, orthogonal Laplacian.
- Scalar/vector face interpolation.
- Collocated pressure/velocity coupling with Rhie-Chow-style flux reconstruction.
- SIMPLE/PISO/PIMPLE-like loops.
- Scalar transport with workspace reuse and matrix/ILU caching.
- Compressible, reacting, radiation, VOF, turbulence and mesh-motion reference pieces.

**Study files**

- `include/cfd/fvm/poly_mesh.hpp`
- `include/cfd/fvm/operators.hpp`
- `include/cfd/fvm/schemes.hpp`
- `include/cfd/fvm/fv_matrix.hpp`
- `src/fvm/pressure_velocity.cpp`
- `src/fvm/collocated_incompressible.cpp`
- `src/fvm/scalar_transport.cpp`
- `docs/PHASE4A_FVM.md`
- `docs/PHASE4B_INCOMPRESSIBLE.md`
- `docs/PHASE4C_COLLOCATED.md`

**Run**

```bash
./build/cfd-solve fvm-operators3d
./build/cfd-solve fvm-taylor-green2d
./build/cfd-solve fvm-collocated-channel
./build/cfd-solve fvm-scalar-transport
```

**Key ideas to learn**

Conservative flux balance, cell/face topology, pressure projection, non-orthogonal correction, bounded reconstruction, under-relaxation, sparse linear solves and operator caching.

## 3. FEM / Elmer-like multiphysics core

**What is implemented**

- Reference elements for Line2, Tri3, Quad4, Tet4, Hex8, Prism6 and Pyramid5.
- 2-D/3-D Poisson, heat, elasticity, electrostatics and DC conduction.
- Axisymmetric elasticity.
- Nonlinear Poisson/reaction solved with Newton/ILU-GMRES.
- Darcy, Brinkman/Stokes, magnetostatics, eddy currents, modal bar dynamics, phase-change and nonlinear-geometry baselines.
- Adaptivity with error indicators and longest-edge refinement.

**Study files**

- `include/cfd/fem/reference_element.hpp`
- `include/cfd/fem/mesh2d.hpp`
- `include/cfd/fem/mesh3d.hpp`
- `src/fem/assembly.cpp`
- `src/fem/poisson3d.cpp`
- `src/fem/elasticity2d.cpp`
- `src/fem/nonlinear_poisson2d.cpp`
- `src/fem/eddy_current2d.cpp`
- `docs/PHASE4E_FEM_CORE.md`
- `docs/PHASE4F_FEM_ADVANCED.md`

**Run**

```bash
./build/cfd-solve fem-poisson3d
./build/cfd-solve fem-nonlinear-poisson
./build/cfd-solve fem-darcy
./build/cfd-solve fem-magnetostatic
./build/cfd-solve fem-modal-bar
```

**Key ideas to learn**

Weak forms, element quadrature, global assembly, essential/natural boundaries, nonlinear residual/Jacobian, error estimation, mesh refinement and modal/eigenvalue checks.

## 4. Time-domain EM / openEMS-like FDTD layer

**What is implemented**

- 1-D Maxwell FDTD with PEC, Mur and PML boundaries.
- Dispersive Debye/Drude/Lorentz ADE material updates.
- Hard/soft sources, probes, DFT monitors and wave-port helpers.
- 3-D Maxwell baseline and VTK output.
- Cylindrical TM baseline.
- CPML/TFSF and lumped RLC smoke cases.

**Study files**

- `include/cfd/solvers/fdtd/maxwell1d.hpp`
- `include/cfd/solvers/fdtd/maxwell3d.hpp`
- `include/cfd/solvers/fdtd/port.hpp`
- `src/fdtd/maxwell1d.cpp`
- `src/fdtd/maxwell3d.cpp`
- `src/fdtd/monitor.cpp`
- `src/fdtd/cylindrical_tm.cpp`
- `docs/PHASE5_FDTD_BASELINE.md`
- `docs/PHASE5B_ADVANCED_FDTD.md`

**Run**

```bash
./build/cfd-solve fdtd-mur1d
./build/cfd-solve fdtd-dispersive1d
./build/cfd-solve fdtd-port-vtk
./build/cfd-solve fdtd-cpml-tfsf
```

**Key ideas to learn**

Yee staggering, CFL condition, absorbing boundaries, material polarization state, source injection, DFT post-processing and port normalization.

## 5. Frequency-domain Maxwell and RF FEM

**What is implemented**

- Driven 1-D complex Maxwell/Helmholtz reference solver.
- Dielectric-loaded 1-D PEC cavity eigenmodes.
- 2-D Tri3 Nedelec edge FEM.
- 3-D Tet4 Whitney/Nedelec edge FEM.
- Complex sparse curl-curl solve through real-block CSR and ILU-GMRES.
- Rectangular wave-port modes and 1 W normalization.
- Resonator Q extraction with dielectric/conductor losses.
- RF adaptive mesh indicators and common-pole rational model reduction.

**Study files**

- `include/cfd/em/frequency_domain.hpp`
- `include/cfd/em/edge_fem2d.hpp`
- `include/cfd/em/edge_fem3d.hpp`
- `include/cfd/em/wave_port.hpp`
- `include/cfd/core/complex_sparse.hpp`
- `src/em/frequency_domain.cpp`
- `src/em/edge_fem2d.cpp`
- `src/em/edge_fem3d.cpp`
- `src/em/wave_port.cpp`
- `src/em/adaptivity3d.cpp`
- `docs/CST_COVERAGE_AUDIT.md`

**Run**

```bash
./build/cfd-solve em-frequency1d
./build/cfd-solve em-edge2d
./build/cfd-solve em-edge3d
./build/cfd-solve em-waveport
```

**Key ideas to learn**

Curl-conforming basis functions, edge orientation, PEC tangential constraints, complex sparse systems, wave-port cutoff, stored/lost energy and Q factor.

## 6. RF / antenna / PEEC / SI-PI-EMC tools

**What is implemented**

- Thin-wire MoM dipoles, arbitrary-orientation wires, endpoint junctions, line loads, PEC ground/image baseline and conductor/dielectric loss accounting.
- Antenna match, efficiency, coupling and optimization helpers.
- N-port S/Y/Z conversions, passivity and causality checks.
- PEEC capacitance and partial-inductance baselines.
- Multiconductor cable/harness RLCG propagation with shield transfer impedance and field coupling.
- Eye/BER signal-integrity metrics, PDN impedance/IR-drop/decoupling selection and EMC waveforms/probes.

**Study files**

- `include/cfd/rf/thin_wire_mom.hpp`
- `include/cfd/rf/antenna.hpp`
- `include/cfd/rf/nport.hpp`
- `include/cfd/rf/peec.hpp`
- `include/cfd/em/cable.hpp`
- `src/rf/thin_wire_mom.cpp`
- `src/rf/network.cpp`
- `src/rf/peec.cpp`
- `src/em/cable.cpp`
- `docs/RF_CIRCUIT_RESEARCH.md`

**Run**

```bash
./build/cfd-solve rf-dipole
./build/cfd-solve rf-microstrip
./build/cfd-solve rf-multiwire
```

**Key ideas to learn**

Green functions, impedance matrices, feed/load constraints, S-parameters, passivity, causality, PEEC equivalents, cable transfer functions, eye openings and PDN impedance.

## 7. SPICE / circuit / compact-model layer

**What is implemented**

- MNA circuit kernel with DC, AC, transient and adaptive transient.
- Sparse MNA stamping with dense fallback.
- Diodes, controlled sources, transformers, sampled TEM lines and behavioral sources.
- BDF2/Gear transient integration and LTE control.
- Periodic steady-state and pole-zero analysis.
- Noise, distortion, `.IC`, `.NODESET`, ASCII RAW I/O.
- Generic DAE device seam and electrothermal/self-heating wrapper.
- OSDI loader baseline.

**Study files**

- `include/cfd/circuit/analysis.hpp`
- `include/cfd/circuit/spice.hpp`
- `include/cfd/circuit/raw.hpp`
- `include/cfd/circuit/osdi_loader.hpp`
- `src/circuit/analysis.cpp`
- `src/circuit/spice.cpp`
- `src/circuit/raw.cpp`
- `tests/test_next_analysis.cpp`
- `tests/test_spice_io.cpp`
- `docs/PHASE11_ADAPTIVE_PSS_PZ.md`

**Run**

```bash
./build/cfd-solve spice-rc
./build/cfd-solve spice-diode
./build/cfd-solve spice-adaptive
./build/cfd-solve spice-pss-pz
```

**Key ideas to learn**

MNA unknown ordering, branch currents, Newton iteration, Jacobians, sparse stamping, charge-based DAE residuals, companion models, small-signal linearization and operating-point dependence.

## 8. Chemistry, electrochemistry and corrosion

**What is implemented**

- Arrhenius/mass-action chemistry networks.
- Reversible kinetics and equilibrium constants.
- Adaptive implicit isothermal reactors.
- Acid/base aqueous pH equilibrium.
- Nernst, Butler-Volmer, Faraday helpers.
- 1-D corrosion current/recession solver.
- 1-D and generic PolyMesh Nernst-Planck transport.
- Scharfetter-Gummel drift-diffusion flux and mixed-potential galvanic solver.

**Study files**

- `include/cfd/chemistry/kinetics.hpp`
- `include/cfd/chemistry/implicit_reactor.hpp`
- `include/cfd/chemistry/aqueous_equilibrium.hpp`
- `include/cfd/electrochemistry/electrochemistry.hpp`
- `include/cfd/solvers/electrochemistry/nernst_planck_poly.hpp`
- `src/electrochemistry/corrosion1d.cpp`
- `src/electrochemistry/nernst_planck_poly.cpp`
- `docs/CHEMISTRY_CORROSION_RESEARCH.md`
- `docs/PHASE4D_NUMERICS_CHEMISTRY.md`

**Run**

```bash
./build/cfd-solve chemistry-reactor
./build/cfd-solve chemistry-equilibrium
./build/cfd-solve electrochem-corrosion1d
./build/cfd-solve electrochem-pnp-poly
./build/cfd-solve electrochem-galvanic
```

**Key ideas to learn**

Stoichiometric invariants, stiffness, equilibrium residuals, electrochemical boundary fluxes, electroneutral vs Poisson formulations and corrosion feedback loops.

## 9. Particle, PIC, plasma and wakefields

**What is implemented**

- Nonrelativistic and relativistic Boris pushers.
- Vay relativistic pusher.
- 1-D electrostatic PIC and 1-D/3V EM-PIC.
- 2-D electrostatic PIC.
- 3-D electrostatic PIC with trilinear CIC deposition and spectral Poisson fields.
- 2-D/3-V centered electromagnetic PIC.
- 3-D/3-V centered electromagnetic PIC with full Ex/Ey/Ez and Bx/By/Bz fields.
- True staggered 2-D Yee EM-PIC baseline with electric-wall/sponge field boundaries and absorbing/specular particle walls.
- Spectral charge-conserving current reconstruction from old/new charge density in 1-D/2-D/3-D baselines.
- Monte-Carlo collision and ionization coupling.
- Secondary-yield reporting at particle walls.
- Wake potential convolution and wake impedance.

**Study files**

- `include/cfd/particle/electromagnetic.hpp`
- `include/cfd/particle/wakefield.hpp`
- `src/particle/electromagnetic.cpp`
- `src/particle/wakefield.cpp`
- `tests/test_v093_em_pic.cpp`
- `tests/test_v094_pic2d.cpp`
- `tests/test_v095_em_pic2d.cpp`
- `tests/test_v096_staggered_pic2d.cpp`
- `tests/test_v097_pic3d.cpp`
- `tests/test_v098_em_pic3d.cpp`

**Run**

```bash
./build/cfd-solve particle-pic1d
./build/cfd-solve particle-em-pic1d
./build/cfd-solve particle-pic2d
./build/cfd-solve particle-em-pic2d
./build/cfd-solve particle-staggered-em-pic2d
./build/cfd-solve particle-pic3d
./build/cfd-solve particle-em-pic3d
```

**Key ideas to learn**

Macro-particles, 1-D/2-D/3-D shape functions, charge/current deposition, field gather, centered and staggered Maxwell curl updates, leapfrog timing, Yee staggering, charge conservation, particle boundaries, MCC collisions and kinetic energy diagnostics.

## 10. Bioelectromagnetics and thermal coupling

**What is implemented**

- RMS electric field to SAR conversion.
- Heterogeneous voxel tissue grid and local SAR averaging.
- SAR to Pennes mesh projection.
- Implicit 2-D Pennes bioheat solver with perfusion, metabolic heat and spatial SAR.
- EM -> thermal -> structural examples through multiphysics couplers.

**Study files**

- `include/cfd/multiphysics/bioheat.hpp`
- `include/cfd/multiphysics/electro_thermal.hpp`
- `src/multiphysics/bioheat.cpp`
- `src/multiphysics/electro_thermal.cpp`
- `docs/PHASE5_8_FDTD_MULTIPHYSICS.md`

**Run**

```bash
./build/cfd-solve bioheat-sar
./build/cfd-solve multiphysics-electrothermal
./build/cfd-solve multiphysics-thermoelastic
```

**Key ideas to learn**

SAR units, tissue density, perfusion sink, implicit heat stepping, source transfer and coupled validation using analytical limiting cases.

## 11. Optics

**What is implemented**

- Sequential ray tracing.
- Paraxial ABCD and Gaussian beams.
- Sellmeier materials, Fresnel, Jones/Stokes polarization.
- Conic/even-asphere and freeform pieces.
- Diffraction, wavefront and optimization helpers.

**Study files**

- `include/cfd/solvers/optics/ray.hpp`
- `include/cfd/solvers/optics/gaussian_beam.hpp`
- `include/cfd/optics/freeform.hpp`
- `include/cfd/optics/wavefront.hpp`
- `src/optics/sequential.cpp`
- `src/optics/gaussian_beam.cpp`
- `docs/PHASE6_OPTICS_BASELINE.md`

**Run**

```bash
./build/cfd-solve optics-gaussian
./build/cfd-solve optics-lens
```

**Key ideas to learn**

Ray/surface intersections, optical path length, paraxial approximations, polarization matrices, spot diagrams and merit functions.

## 12. Runtime, validation and packaging

**What is implemented**

- Shared CSR, PCG, BiCGStab, ILU(0), restarted GMRES and complex sparse wrapper.
- CPU backend/profiler/task-graph pieces.
- Workflow schemas, provenance and validation catalog.
- Integration tracker and release docs.

**Study files**

- `include/cfd/core/csr_matrix.hpp`
- `include/cfd/core/iterative_solvers.hpp`
- `include/cfd/core/complex_sparse.hpp`
- `include/cfd/workflow/validation_catalog.hpp`
- `src/workflow/validation_catalog.cpp`
- `scripts/integration_status.py`
- `docs/VALIDATION.md`
- `docs/SOURCE_MAP.md`

**Run**

```bash
./build/cfd-validate
./scripts/integration_status.py
ctest --test-dir build --output-on-failure
```

## Recommended learning order

1. Sparse linear algebra and test style.
2. FVM scalar transport and pressure projection.
3. LBM descriptors and streaming.
4. FEM Poisson/heat/elasticity.
5. FDTD 1-D Maxwell.
6. Frequency-domain Maxwell and Nedelec edge FEM.
7. Circuit MNA and transient companion models.
8. Thin-wire MoM and S-parameters.
9. Electrochemistry and corrosion.
10. PIC particle push/deposition/gather.
11. Multiphysics couplers.
12. Performance, MPI and GPU backends.

## Practical exercises

- Change the `particle-staggered-em-pic2d` time step until the CFL guard rejects the run; explain why.
- Move the `particle-pic3d` particles across periodic boundaries and verify total deposited charge remains invariant.
- Add one new smoke case that reports a physically meaningful scalar and a pass/fail criterion.
- Add a manufactured solution test before adding any new PDE solver.
- For every new boundary condition, test a limiting case: zero source, constant field, symmetry or known analytical reflection.
- Keep the integration tracker honest: mark an item only after there is executable validation.


## v0.9.9 study addendum - staggered 3-D Yee EM-PIC

Study files:
- API: `include/cfd/particle/electromagnetic.hpp` (`StaggeredElectromagneticPic3D`).
- Implementation: `src/particle/electromagnetic.cpp`.
- Regression: `tests/test_v099_staggered_pic3d.cpp`.
- CLI: `cfd-solve particle-staggered-em-pic3d`.

Learning exercises:
1. Compare `ElectromagneticPic3D` with `StaggeredElectromagneticPic3D`: identify the nodal-grid curl locations versus the Yee Ex/Ey/Ez and Bx/By/Bz locations.
2. Derive the forward Faraday update for Bx, By and Bz from the component offsets.
3. Derive the backward Ampere update for Ex, Ey and Ez and verify the code uses the matching dual differences.
4. Toggle `GridBoundaryMode3D::electric_wall` and `absorbing_sponge` in the CLI-style case and inspect how boundary energy changes.
5. Replace the spectral continuity current with a local deposition method as the next production exercise.

## v0.10.0 study addendum - local 3-D current deposition

New source/API entry points:

- `include/cfd/particle/electromagnetic.hpp`
  - `deposit_charge_conserving_current_3d_local(...)`
  - `CurrentDeposition3DMode`
  - `StaggeredElectromagneticPic3DConfig::current_deposition`
  - `GridBoundary3DConfig::sponge_polynomial_order`
- `src/particle/electromagnetic.cpp`
  - finite-volume current reconstruction and selectable staggered-PIC coupling
- `tests/test_v0100_local_current.cpp`
  - continuity residual, solver integration and sponge-profile regressions
- CLI: `cfd-solve particle-local-current3d`

Exercises:

1. Move one macro-particle by a small known displacement and compare spectral vs local current distributions.
2. Recompute the finite-volume divergence `(Jx[i]-Jx[i-1])/dx + ...` by hand for a small 3x3x3 grid.
3. Change `sponge_polynomial_order` from 1 to 4 and inspect how quickly fields are damped across the sponge layer.
4. Switch `StaggeredElectromagneticPic3DConfig::current_deposition` between spectral and local modes, then compare continuity residuals and current smoothness.


## v0.10.1 study addendum - particle sorting and guard halos

New source/API entry points:

- `include/cfd/particle/electromagnetic.hpp`
  - `particle_cell_index_3d(...)`
  - `sort_particles_by_cell_3d(...)`
  - `ParticleCellSort3DResult`
  - `classify_particle_guard_halos_3d(...)`
  - `ParticleGuardHalo3DConfig` / `ParticleGuardHalo3DResult`
  - `StaggeredElectromagneticPic3DConfig::sort_particles_by_cell`
  - `StaggeredElectromagneticPic3DDiagnostics::particle_sort_passes`
- `src/particle/electromagnetic.cpp`
  - stable cell sorting, guard-halo classification and sorted staggered-PIC integration
- `tests/test_v0101_particle_sorting.cpp`
  - cell range/recovery, halo face/edge/outside and solver sorting regressions
- CLI: `cfd-solve particle-sort-halo3d`

Exercises:

1. Build a small 4x4x4 particle set and verify `cell_offsets[cell+1]-cell_offsets[cell]` equals `cell_counts[cell]`.
2. Put one particle near an x/y/z corner and confirm it appears in multiple halo lists.
3. Enable `sort_particles_by_cell` in `StaggeredElectromagneticPic3D`, step once, then verify particle cell indices are nondecreasing.
4. Explain why sorting must happen after current deposition for a step, not between the old/new particle snapshots.
5. Sketch how the guard-halo lists map to future MPI send buffers and GPU particle tiles.


## v0.10.2 study addendum - PIC domain decomposition and guard exchange

Source files:

- `include/cfd/particle/electromagnetic.hpp`
- `src/particle/electromagnetic.cpp`

Core APIs to study:

- `make_pic_domain_grid_3d` - builds Cartesian subdomain descriptors from a global grid and domain counts.
- `locate_pic_domain_3d` - maps physical macro-particle positions to destination subdomains.
- `plan_particle_domain_migration_3d` - buckets particles by destination domain and reports nonperiodic outside particles.
- `exchange_scalar_guard_cells_3d` - creates ghost-padded scalar blocks from neighbouring decomposed subdomains.

Validation entry points:

- `tests/test_v0102_domain_exchange.cpp`
- `./build/cfd-v0102-domain-exchange-tests`
- `./build/cfd-solve particle-domain-exchange3d`

Exercises:

1. Change the global grid to an uneven size such as `11x7x5` and verify every global cell is covered exactly once by the subdomains.
2. Add a particle just outside each side of the nonperiodic domain and confirm it appears in `outside_indices`.
3. Fill each local scalar block with a unique global-cell formula and verify every ghost cell copies the correct neighbour value.
4. Replace the scalar field with one component of `Ex/Ey/Ez` and sketch how the same pack/scatter contract becomes an MPI guard-cell exchange.


## v0.10.3 study addendum - communicator-ready PIC exchange messages

Source/API entry points:

- `include/cfd/particle/electromagnetic.hpp`
  - `ParticleMigrationMessage3D` / `PackedParticleMigration3D`
  - `pack_particle_migration_messages_3d(...)`
  - `apply_particle_migration_messages_3d(...)`
  - `ScalarGuardCellMessage3D` / `ScalarGuardCellValue3D`
  - `pack_scalar_guard_cell_messages_3d(...)`
  - `apply_scalar_guard_cell_messages_3d(...)`
- `src/particle/electromagnetic.cpp`
  - deterministic grouping of source-to-destination particle and guard-cell payloads
- `tests/test_v0103_comm_exchange.cpp`
  - message/apply parity with serial exchange, periodic wrapping and nonperiodic exterior fill regressions
- CLI: `./build/cfd-solve particle-comm-exchange3d`

What to learn:

1. Separate **ownership** from **transport**. `pack_particle_migration_messages_3d` decides which destination owns each particle, but it does not send anything.
2. Guard-cell messages carry destination padded coordinates, so a future MPI layer can remain a byte transport and does not need to understand ghost-index geometry.
3. Apply functions are deterministic reducers. This makes serial, MPI and GPU implementations comparable against the same contract.
4. Nonperiodic exterior cells are intentionally left as the caller-supplied exterior fill value when no message arrives.
5. The next implementation step is a real MPI wrapper that maps each message to rank-local send/receive buffers.

Exercises:

- Add a 2-cell guard width and verify the number of guard-message entries grows as expected.
- Change the domain grid to `3x2x2` and confirm that edge/corner guard messages still reproduce the serial exchange.
- Add one out-of-domain particle to a nonperiodic run and inspect the `outside_particles` source-domain/source-index record.
- Sketch the MPI calls needed to replace direct `apply_*_messages` with send/receive followed by apply.


## v0.10.4 study addendum - PIC transport layer before MPI

**Source/API:** `include/cfd/particle/electromagnetic.hpp`, `src/particle/electromagnetic.cpp`

Study path:
1. Start with `PicDomainGrid3DConfig` and `make_pic_domain_grid_3d` to understand how global cells are split into subdomains.
2. Read `make_pic_rank_topology_3d`; it groups one or more subdomains onto logical ranks before any MPI dependency is introduced.
3. Compare `ParticleMigrationMessage3D` and `ScalarGuardCellMessage3D` with `PicTransportEnvelope3D`; the envelope adds source/destination rank metadata but preserves solver payloads.
4. Inspect `InMemoryPicTransport3D`; each destination rank has a deterministic inbox, which models the contract that MPI will later implement with real send/receive calls.
5. Run `./build/cfd-solve particle-transport3d` and `./build/cfd-v0104-transport-tests` to verify that transported exchange matches the serial migration/guard exchange.

Exercise: change the rank topology from `3 x 1 x 1` to `1 x 2 x 1` in the v0.10.4 test and observe how same-rank and remote-rank message counts change while the final guard blocks remain identical.


## v0.10.5 study addendum - serialized PIC transport

**What to read in code**

- `include/cfd/particle/electromagnetic.hpp`: `PicSerializedEnvelope3D`, `PicSerializedExchangePlan3D`, `serialize_pic_transport_envelope_3d`, `deserialize_pic_transport_envelope_3d`, and the guarded MPI transport declaration.
- `src/particle/electromagnetic.cpp`: deterministic binary packing of particle migration and scalar guard-cell payloads.
- `tests/test_v0105_serialized_transport.cpp`: roundtrip, exchange-planning and malformed-envelope checks.

**Concept checklist**

1. Start from `ParticleMigrationMessage3D` and `ScalarGuardCellMessage3D`.
2. Wrap them as rank-addressed `PicTransportEnvelope3D` messages.
3. Serialize each envelope into an opaque byte record.
4. Bucket serialized records by destination rank.
5. Decode records and apply them with the existing migration/guard reconstruction helpers.

**Exercise**

Change the rank topology in `particle-serialized-transport3d` from `3x1x1` to `1x2x1`; confirm the same particles and guard values arrive, while same-rank vs remote-rank diagnostics change.


## v0.10.5 study addendum - plasma chemistry and breakdown thresholds

**Code path**

- `include/cfd/particle/plasma.hpp`
- `src/particle/plasma.cpp`
- `tests/test_v0105_plasma.cpp`

**Implemented ideas**

- Multi-species number-density state with charge/mass metadata.
- Mass-action plasma reactions with stoichiometric reactant/product coefficients.
- Charge-density diagnostics before and after each chemistry step.
- Paschen-law gas-breakdown voltage/field estimate.
- Parallel-plate resonant multipactor field and impact-energy estimate.

**Exercises**

1. Change the ionization rate coefficient in `particle-plasma-chemistry` and observe the electron/ion growth.
2. Sweep pressure*gap in `particle-breakdown-threshold` to recover the Paschen minimum qualitatively.
3. Compare first-, third- and fifth-order multipactor estimates and note how transit time changes resonant field.


## v0.10.6 study addendum - distributed PIC exchange round and EM field guards

**Code path**

- `include/cfd/particle/electromagnetic.hpp`
  - `ElectromagneticFieldBlocks3D`
  - `ElectromagneticGuardedFieldBlocks3D`
  - `exchange_electromagnetic_field_guard_cells_3d(...)`
  - `DistributedPicExchangeRound3D`
  - `run_serialized_distributed_pic_exchange_round_3d(...)`
- `src/particle/electromagnetic.cpp`
  - six-component EM field guard exchange built from the scalar guard contract
  - serialized rank-local distributed exchange orchestration
- `tests/test_v0106_distributed_round.cpp`
  - direct scalar-vs-vector guard parity
  - serialized distributed round vs direct guard/migration reconstruction

**Concept checklist**

1. A distributed EM-PIC step needs both particle migration and field guard exchange.
2. Scalar guard exchange is the smallest validated primitive; v0.10.6 lifts it to six EM components.
3. A rank-local exchange round should be validated without MPI first: pack, envelope, serialize, deserialize, apply and compare with direct reference exchange.
4. Serialized byte counts and destination-rank message counts are useful CI diagnostics before replacing the in-memory path with real MPI.
5. Keeping this orchestration separate from the Yee update makes it easier to test the communicator contract independently of field physics.

**Exercises**

- Run `./build/cfd-solve particle-distributed-round3d` and inspect serialized message/byte counts.
- Run `./build/cfd-solve particle-field-guards3d` and compare guarded-value counts with the scalar guard width.
- Change the logical rank topology in `tests/test_v0106_distributed_round.cpp` and confirm the final particles/guard blocks are unchanged while remote message counts change.
- Extend the field-guard wrapper to carry charge/current components and compare its output with six independent scalar guard exchanges.

## v0.10.7 study addendum - distributed PIC timestep bridge

### What was added

- API: `run_serialized_distributed_staggered_pic_step_3d`.
- Types: `DistributedStaggeredPicStep3DConfig` and `DistributedStaggeredPicStep3D`.
- Test: `tests/test_v0107_distributed_step.cpp`.
- CLI: `cfd-solve particle-distributed-step3d`.

### How to study it

1. Start from the v0.10.6 exchange round and identify its inputs: particle buckets, scalar field blocks, domains and topology.
2. Compare that to the v0.10.7 timestep bridge: particles are first pushed using domain-local Ex/Ey/Ez/Bx/By/Bz samples, then ownership migration and guard reconstruction are performed.
3. Run `cfd-solve particle-distributed-step3d` and check the printed particle counts, maximum displacement, message count and guarded-domain count.
4. Inspect `tests/test_v0107_distributed_step.cpp` to see why a ballistic zero-charge particle is useful: it validates timestep/migration ordering without adding Lorentz-force numerical ambiguity.

### Exercises

- Change the test to apply a uniform electric field and verify that particle velocity changes before migration.
- Add a non-periodic case and confirm outside-particle reporting remains deterministic.
- Extend the bridge to exchange all six field components through serialized envelopes rather than using direct in-process field guard exchange.

## v0.10.8 Geant4-inspired particle transport baseline

Source files:

- `include/cfd/particle/transport.hpp`
- `src/particle/transport.cpp`
- `tests/test_v0108_geant4_transport.cpp`

What to study:

- `TransportTrack` is the persistent particle state.
- `TransportStep` is one accepted movement with pre/post state, limiter,
  deposited energy and secondaries.
- `TransportProcess` is a clean-room process hook. The first baseline supports
  continuous stopping power, discrete secondary-producing interactions and a
  simple absorber threshold.
- `transport_one_step` mirrors the Geant4-style order: locate region, select the
  minimum allowed step, apply continuous losses, apply a post-step discrete
  process, update track state and scoring data.

Exercises:

1. Change `physical_interaction_length_m` and verify the discrete process no
   longer limits the first step.
2. Raise `production_cut_energy_ev` and verify secondaries are suppressed while
   energy deposition remains conservative.
3. Add a third slab and verify `geometry_boundary` limits track motion at the
   new interface.

## v0.10.9 - Geant4/SU2/csauto clean-room transport and campaign workflow

This checkpoint extends the Geant4-inspired particle transport seam from deterministic process arbitration to a small stochastic transport stack that remains independent of Geant4 source code. The implemented pieces are:

- `TransportPhysicsList`, which bundles process instances and can override the production cut in a copied world for a run.
- `TransportRandom` plus `sample_exponential_interaction_length_m`, which turns a mean free path into a sampled discrete process distance using a deterministic LCG for regression tests.
- `TransportRegionBvh`, built by median-splitting region AABBs, used to compare accelerated region lookup against the existing linear lookup.
- `TransportSensitiveDetector`, `TransportHitCollection` and `TransportDoseGrid3D`, which collect per-step detector hits, accumulate deposited energy as dose, derive dose-rate/SAR and project it to a 2-D Pennes bioheat SAR map.

SU2 was used as a clean-room workflow reference for PDE-constrained optimization organization: solver/iteration/numerics/output separation, Python-side study scripts, adjoint/finite-difference gradient patterns and regression-driven validation. In this project, the corresponding lightweight pieces are `finite_difference_gradient`, `gradient_descent_update`, and `SolverAdapterDescriptor` for workflow-facing solver declarations.

csauto was used as a clean-room campaign-automation reference: one table row becomes one case, templates use placeholders and conditional sections, the core talks through a solver-adapter boundary, and run state is summarized through a registry. The implemented C++ equivalents are `generate_factorial_campaign`, `generate_latin_hypercube_campaign`, `render_campaign_template`, `validate_solver_adapter_descriptor`, `update_case_status` and `summarize_campaign_registry`.

Try the new focused cases:

```bash
./build/cfd-solve particle-transport-dose-bvh
./build/cfd-solve particle-campaign-doe
./build/cfd-v0109-transport-campaign-tests
```

Limitations: the dose scorer uses axis-aligned regions and a regular grid, the stochastic transport only handles the compact process interfaces already implemented here, and the campaign layer is an in-process reference library rather than a persistent web dashboard or external solver launcher.

## v0.11.0 - persistent campaign execution workflow

Source files:

- `include/cfd/workflow/campaign.hpp`
- `src/workflow/campaign.cpp`
- `tests/test_v0110_campaign_execution.cpp`

Study path:

1. Start with `generate_factorial_campaign` or `generate_latin_hypercube_campaign` to create rows.
2. Pass those rows plus `CampaignTemplateFile` entries to `write_campaign_case_folders`.
3. Inspect the generated case folders: rendered input files, `doe_row.csv`, `registry.tsv` and `template_manifest.tsv`.
4. Build a `SolverRunRequest` and inspect the returned `SolverCommandPlan` before adding any real process launcher.
5. Feed a synthetic solver log into `parse_residual_history`, `parse_performance_metrics` and `detect_solver_outcome_from_log`.
6. Use parsed residuals/objectives to update the registry, then run `run_gradient_descent_campaign` to study a compact finite-difference optimization loop.

Exercises:

- Add a new runtime enum value for a local job queue and test only the command-plan contract first.
- Extend the residual parser with a solver-specific grammar while preserving the generic parser tests.
- Add JSON export next to `registry.tsv` and compare roundtrip robustness.


## v0.11.1 workflow study: local campaign execution

The campaign layer now has a minimal native execution path:

1. Generate or persist a campaign with `write_campaign_case_folders`.
2. Validate the selected solver/case combination with `doctor_solver_runtime`.
3. Launch a native executable with `run_local_solver_case`, which uses argv execution and captures stdout/stderr logs.
4. Discover residual/history and performance/timing files with `discover_campaign_outputs`.
5. Run all pending cases and update `registry.tsv` with `run_local_campaign`.

The regression target `cfd-v0111-local-campaign-runner-tests` builds a stub executable, runs a two-case campaign, verifies residual/performance parsing, and checks registry status/objective/iteration updates.


## v0.11.2 workflow study: control directives and scheduler refresh

The campaign workflow now has a control/monitoring seam on top of the local runner:

1. Use `validate_campaign_control_request` to ensure the selected solver adapter advertises both the `control` capability and the requested action.
2. Use `write_campaign_control_directive` to write deterministic `.cfd_control/<action>.directive` files for stop, extend, checkpoint or flush requests.
3. Use `discover_campaign_control_directives` to enumerate pending steering requests from a generated case folder.
4. Use `parse_slurm_queue_table` to convert compact scheduler output into typed `SchedulerJobRecord` values.
5. Use `apply_scheduler_records_to_registry` and `refresh_campaign_status_from_outputs` to update registry state from scheduler status or solver logs without relaunching cases.

The regression target `cfd-v0112-campaign-control-tests` checks control capability validation, directive persistence/discovery, scheduler-state mapping and output-based registry refresh.

## v0.11.3 workflow/deployment study notes

Study the new workflow deployment seam by starting with `include/cfd/workflow/deploy.hpp`. The core exercise is to generate a factorial campaign, persist it with `write_campaign_case_folders`, create several `CampaignServerDescriptor` entries, and call `plan_multi_server_campaign`. Inspect the resulting assignments, native/Docker solver plans, and generated `multiserver_commands.sh` before attempting any real remote execution.


## v0.11.4 supervised campaign execution layer

Study these functions together with the v0.11.3 placement planner:

- `plan_multiserver_execution(...)` converts case-to-server assignments into health checks and supervised jobs.
- Each `RemoteSupervisedJobPlan` owns launch, status, cancel and fetch-log command plans plus `.cfd_run` metadata paths.
- `parse_remote_job_status_table(...)` accepts tab-separated status rows.
- `apply_remote_job_status_to_registry(...)` maps distributed job state back into the campaign registry.

The design is intentionally deterministic so remote execution semantics are tested before adding real SSH/Docker daemon orchestration.


## v0.11.5 - Multi-server supervision hardening

New workflow APIs in `cfd/workflow/deploy.hpp` convert the v0.11.4 supervised job plan into operator/dashboard artifacts:

- `plan_multiserver_supervision(...)` creates access probes, retry launch plans, stdout/stderr tail plans and dashboard case summaries.
- `redact_sensitive_command_display(...)` removes obvious inline secret values from display strings.
- `build_multiserver_dashboard_json(...)` emits a small JSON summary for a future UI or service endpoint.
- `write_multiserver_supervision_files(...)` persists access-check, retry-launch and log-tail shell scripts plus JSON/TSV summaries.

The CLI smoke case `particle-multiserver-supervision` exercises the complete offline path. The focused regression target is `cfd-v0115-multiserver-supervision-tests`.

## v0.11.6 - Multi-server controller API scaffold

New workflow APIs in `cfd/workflow/deploy.hpp` build a controller layer above the supervision artifacts:

- `plan_multiserver_controller(...)` derives read and mutation routes from `MultiServerSupervisionPlan` metadata.
- `build_multiserver_controller_status_json(...)` converts dashboard cases into a stable status payload.
- `build_multiserver_controller_openapi_json(...)` documents the generated API surface.
- `build_multiserver_controller_script(...)` creates a Python standard-library controller script.
- `write_multiserver_controller_files(...)` persists the controller script, OpenAPI JSON, status JSON, routes TSV, env example, README and launch script.

Study path:

1. Generate a campaign and multi-server supervision plan as in v0.11.5.
2. Configure `RemoteControllerConfig` with host, port, token env var, tail length and allowed control actions.
3. Inspect `controller/routes.tsv` and `controller/openapi.json` before starting the service.
4. Run the generated controller in read-only mode first, then enable token-gated mutations.
5. Verify that POST control requests write `.cfd_control/*.directive` files rather than calling solver-specific control code directly.

The CLI smoke case `particle-multiserver-controller` exercises the offline generation path. The focused regression target is `cfd-v0116-multiserver-controller-tests`.

### v0.11.7 generated controller dashboard

The controller scaffold now includes a generated browser dashboard:

- `build_multiserver_controller_dashboard_html(...)` creates a dependency-free operator page.
- `build_multiserver_controller_dashboard_js(...)` polls health/status, renders case rows, fetches stdout/stderr tails and posts token-authorized control actions.
- `build_multiserver_controller_dashboard_css()` emits compact dark-mode styling for summary cards and case states.
- `build_multiserver_controller_events_ndjson(...)` emits a newline-delimited event snapshot suitable for downstream dashboard ingestion.
- `write_multiserver_controller_files(...)` persists the static assets alongside OpenAPI/status/routes/env artifacts.

Use `particle-multiserver-dashboard` as the smoke case and `cfd-v0117-controller-dashboard-tests` as the focused regression target.

### v0.11.8 server-sent controller status events

The controller's live-update path remains built on the same generic status document:

- `RemoteControllerConfig::expose_event_stream` controls route exposure.
- `event_stream_heartbeat_seconds` bounds idle time between SSE frames.
- `StreamingHandler` emits named `status` events only when the normalized JSON payload changes and sends heartbeat comments otherwise.
- The dashboard opens an `EventSource` connection and continues periodic polling as a compatibility fallback.

Inspect `/api/events/stream` in `routes.tsv` and `openapi.json`, then run `cfd-v0118-controller-sse-tests` to validate route gating, JavaScript wiring, heartbeat configuration and generated Python syntax.

### v0.12.0 production controller deployment

The generated controller now has an explicit production boundary:

- TLS fails closed when enabled certificate/key paths are absent.
- External token files store only SHA-256 digests with viewer/operator/admin roles.
- Protected routes use constant-time token verification.
- Audit and probe events are written to a bounded rotating NDJSON history.
- Live probes reuse the generated access-check script with timeout and output limits.
- systemd and reverse-proxy examples document persistent operation and network hardening.

Run `cfd-v0120-production-controller-tests` and inspect the generated `env.example`, token template, service unit and proxy configuration before adapting them to a real host.

### v0.13.0 XDMF/HDF5 output and Python bindings

Start with `write_xdmf_inline(...)`: it validates a `Mesh2D`, writes triangle connectivity, XY coordinates and optional nodal scalar data. Then compare `write_hdf5_xdmf(...)`, where the XDMF document references datasets in a companion HDF5 file.

The Python layer intentionally binds the stable C ABI rather than C++ classes. `cfd_solvers_capi` accepts flat coordinate/connectivity arrays, reconstructs the validated mesh and invokes the same XDMF writer. `python/cfd_solvers/bindings.py` performs Python-side shape conversion and error translation using only `ctypes`.

Run `cfd-v0130-xdmf-python-tests` and `cfd-v0130-python-binding-smoke`. On a host with HDF5 development files, configure with `CFD_ENABLE_HDF5=ON` and inspect the generated datasets with `h5dump` or an HDF5 viewer.

## v0.14.0 - common HPC runtime completion

Study the three Phase-1 closure paths together because they meet at sparse iterative solvers:

1. `cfd/core/numa.hpp` and `src/core/numa.cpp` separate topology discovery, deterministic placement planning, thread binding and first touch. Compare compact vs spread plans on a multi-socket host, and note that Linux discovery honors the process's pre-existing affinity mask.
2. `cfd/core/sycl_sparse.hpp` and `src/sycl/sparse_linalg_sycl.cpp` keep CSR data and Krylov vectors in device USM. Trace one CG iteration: SpMV -> dot reduction -> vector update -> residual reduction -> direction update.
3. `cfd/distributed/mpi_sparse.hpp` and `src/distributed/mpi_sparse.cpp` turn global CSR halo columns into an ownership-aware request plan. The expensive index-discovery exchange occurs once; repeated Krylov multiplies exchange values only.

A useful exercise is to instrument bytes moved per distributed SpMV and compare the cached sparse exchange against an `MPI_Allgatherv` full-vector baseline. For SYCL, compare the current generic row-per-work-item SpMV with a later device-tuned segmented/warp-aware kernel while preserving the same public API.

### v0.14.1 LBM advanced-model study path

Read `compressed_pull.hpp`, `particles.hpp`, `free_surface.hpp` and `advanced_d2q9.hpp` as independent numerical boundaries before combining them. The focused v0.14.1 regression intentionally checks analytic/simple invariants: half-way curved-wall equivalence, particle action/reaction balance, free-surface volume conservation, Q-criterion on known gradients, compressed-vs-float error gates and Taylor-Green mass/energy behavior for LES and cumulant collision.

The optional SYCL voxelizer is designed to reproduce the CPU mask contract; hardware CI is still required before treating accelerator parity as production-qualified.


### v0.14.2 portability/autotuning study path

Start with `cfd/core/autotune.hpp`: a tuning record combines a portable device key with a kernel-group name, launch/fusion parameters and an observed score. Follow `KernelTuningDatabase::autotune(...)` through save/load to see how benchmarking stays separate from kernel code.

Then inspect `cfd/distributed/device_assignment.hpp`. `assign_device_group(...)` distributes visible accelerators across local ranks, `split_rank_work(...)` maps a contiguous range onto the selected ordinals, and the SYCL helper materializes one queue per selected device. This is the boundary higher solver layers use for multi-GPU submission.

Finally, read `cfd/distributed/repartition.hpp`: weighted prefix load produces new contiguous ownership boundaries, migration segments describe old/new overlap, and the decision helper adds hysteresis. The focused regression is `cfd-v0142-portability-tests`.

### v0.15.0 FVM second-order time integration and pressure AMG seam

Start with `cfd/fvm/temporal.hpp`, then inspect `ScalarTransport::step()`. The reusable spatial matrix is combined with a scheme-specific transient diagonal and history RHS. BDF2 bootstraps with Euler; Crank-Nicolson combines the implicit new-state spatial operator with the old-state spatial residual.

Next follow `CollocatedIncompressible::momentum_predictor(...)`. PISO and PIMPLE pass a frozen physical time level into every coupling/outer iteration so PIMPLE does not accidentally advance history several times in one timestep. SIMPLE deliberately bypasses the physical second-order history and retains pseudo-time Euler behavior.

Finally inspect `correct_pressure(...)`: the pressure operator remains matrix-free, while CG can be replaced by PCG with Jacobi or with a supplied AMG/multigrid V-cycle. Run `cfd-v0150-fvm-temporal-pressure-tests` to exercise all three integration seams.

### v0.15.1 solid heat, coupled thermal sources and FVM restart I/O

Start with `SolidHeatConduction` and trace how physical `rho`, `cp`, `k` and heat generation are converted onto the existing implicit scalar-transport equation. The source-only energy-balance regression is useful because it has an exact integral result independent of spatial discretization.

Next inspect `ThermalTransport::set_reactive_radiative_source(...)`. Chemistry and radiation are deliberately composed at the energy-source boundary. The radiation model is held by `shared_ptr<const RadiationSourceModel>` so a temporary or externally destroyed model cannot leave the time-step callback dangling.

Finally inspect `cfd/io/fvm_checkpoint.hpp` and `src/io/fvm_checkpoint.cpp`. The public API exposes meshes, typed field arrays and metadata rather than HDF5 handles. Compare the fail-closed non-HDF5 regression with the HDF5-enabled round-trip regression, then inspect a generated file with `h5dump` to follow the `/Mesh`, `/Meta` and `/Fields` groups.

### v0.15.2 transported RANS, Reynolds stress and SST-DES study path

Start with `cfd/solvers/fvm/rans_transport.hpp` and the internal `advance_scalar(...)` path in `src/fvm/rans_transport.cpp`. It is the common numerical kernel for turbulence variables: implicit upwind convection, harmonic variable diffusion, semi-implicit sinks, the shared ILU0/GMRES stack and Euler/BDF2/Crank-Nicolson history.

Then compare the three eddy-viscosity transports. SA builds `fv1/fv2`, modified strain, nonlinear-gradient production and wall destruction around one transported working variable. k-epsilon uses production and dissipation coupling between two equations. SST adds F1/F2 blending, cross diffusion and a viscosity limiter; enable `des_enabled` and inspect `hybrid_dissipation_factor()` to see where the modeled length scale switches toward LES behavior.

Finally inspect `ReynoldsStressTransport`. Follow the velocity gradient into the exact tensor-production term, then the LRR-style pressure-strain source and the six segregated tensor equations. After the solve, trace `enforce_realizability()` to see how a numerically transported symmetric tensor is mapped back to a valid covariance tensor.

Run `cfd-v0152-rans-transport-tests`; the most useful invariants are exact simple-shear strain, positivity of scalar turbulence variables, F1/F2 bounds, stronger DES dissipation on a fine grid and positive-semidefinite Reynolds stresses.

### v0.15.3 conservative topology change and FVM AMR study path

Start with `AdaptiveHexMesh` in `cfd/fvm/adaptive_mesh.hpp`. Compare the leaf-cell bounds/root/level/lineage metadata with the much lighter `PolyMesh` representation. Then follow `AdaptiveHexMesh::poly_mesh()` and identify how a coarse/fine interface is expressed as several owner/neighbour faces without adding explicit hanging nodes to the solver-facing mesh.

Next inspect `MeshTopologyOperation`, `MarkedHexRefinement` and `MarkedHexCoarsening`. Trace one root cell through 2x2x2 splitting and back through complete-sibling merging, then enable 2:1 balancing and follow the recursive coarse-neighbour refinement.

Finally inspect `build_topology_change_map(...)`, `conservative_adaptive_remap(...)` and `adapt_scalar_field(...)`. The key invariant is local overlap-volume conservation, followed by estimator/marking/topology separation. Run `cfd-v0153-fvm-amr-tests` and check the field integral before and after every topology change.

### v0.15.4 characteristic WENO Euler study path

Start with the new enums and controls in `cfd/solvers/fvm/compressible1d.hpp`. Compare the legacy `muscl_minmod + forward_euler` combination with `characteristic_weno5 + ssprk3`; the old path is intentionally still the default so existing callers do not silently change numerics.

Then follow `roe_basis(...)` in `src/fvm/compressible1d.cpp`. Derive the three 1-D Euler right eigenvectors from Roe-averaged velocity, total enthalpy and sound speed. Inspect the small 3x3 inverse and trace a five-cell conservative stencil through left-eigenvector projection, scalar WENO5 reconstruction and right-eigenvector recovery.

Next inspect `weno5_left(...)`: identify the three third-order candidate polynomials, the Jiang-Shu smoothness indicators and ideal weights 0.1/0.6/0.3. Reverse the stencil to see how the right face state uses the same scalar routine.

Finally follow the SSPRK3 stages and `pressure_jump_sensor()`. Run `cfd-v0154-compressible-weno-tests`; the most important numerical check is the periodic entropy wave, because it tests reconstruction order without contaminating the result with a discontinuity. The Sod case is the robustness complement, not the convergence case.

### v0.15.5 GPU-residency study path

Start with `cfd/core/device_residency.hpp` and treat transfer bytes and synchronization points as first-class performance counters. Reset the counters after setup and identify every operation that forces host visibility.

Then compare the SYCL LBM initialization and diagnostics with the previous host-staged shape. Follow `refresh_device_macroscopic()` to see how derived `rho/u` remain available to downstream GPU consumers without requiring a host copy, and compare `total_mass()` with a full macroscopic download.

Next inspect `cfd/core/sycl_sparse.hpp` and `src/sycl/sparse_linalg_sycl.cpp`. The CSR matrix was already persistent; v0.15.5 also keeps Krylov work vectors allocated across solves and adds direct device-USM SpMV/CG calls. The host-span compatibility API still stages vectors, while the direct device API avoids that volume transfer. Per-iteration scalar convergence checks remain a synchronization boundary, so this is not yet a fully resident FVM pressure loop.

Finally read `FLUIDX3D_GPU_RESIDENCY_RESEARCH.md`. Apply its acceptance test to each solver family: after initialization and before output, a qualified device time loop should show zero explicit host-field transfer bytes. Keep physical hardware parity separate from source-level/SYCL syntax validation.


### v0.15.6 GPU-R1 study path

Start with `include/cfd/solvers/lbm/esoteric_pull_sycl.hpp` and trace the optional resident acceleration field from `device_local_acceleration()` into the Guo-force collision term. Then inspect `resident_multiphysics_sycl.hpp`: thermal D3Q7 reads the device macroscopic view, free-surface VOF computes flux divergence entirely on-device, particle coupling atomically scatters equal-and-opposite drag into the same force field, and Q-criterion consumes the resident velocity scratch. The intended timestep has no full-field host transfer between these stages.

When extending this path, preserve the same rule: setup/output may cross the host boundary, but ordinary no-output timesteps should not. Prefer persistent allocations and device views over convenience APIs that materialize host vectors.
