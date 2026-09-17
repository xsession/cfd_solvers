# CST-style coverage audit and expansion map

This document maps the current `cfd_solvers` implementation against the CST-style domain coverage supplied for the continuation work. It is a clean-room capability map, not an assertion of product parity. The goal is to identify which numerical families already exist, which have reference baselines, and which still need production-grade implementation.

Status meanings:

- **Implemented**: a public API plus deterministic regression/validation exists in this repository.
- **Baseline**: a useful numerical reference exists, but the dimensionality, model breadth, or validation set is not yet comparable to a mature specialist solver.
- **Open**: no solver of the requested class exists yet.

## Coverage matrix

| Domain | Current repository state | Status | Principal implementation |
| --- | --- | --- | --- |
| Electrostatics | 2-D FEM electrostatic potential/field solve | Implemented baseline | `src/fem/electromagnetics2d.cpp` |
| Magnetostatics | 2-D out-of-plane vector-potential FEM with reconstructed B | Implemented baseline | `src/fem/magnetostatics2d.cpp` |
| Low-frequency harmonic EM | 2-D harmonic eddy-current formulation | Implemented baseline | `src/fem/eddy_current2d.cpp` |
| Low-frequency transient EM / machines | Nonlinear B-H magnetostatics, stranded-coil flux-linkage/circuit seam and Maxwell-stress force/torque now exist; transient A-phi, nonlinear eddy current and moving-band machines remain open | Baseline | `src/fem/magnetostatics2d.cpp` |
| Full-wave time-domain EM | 1-D/3-D Cartesian Yee FDTD, cylindrical TM, materials, PML/Mur, ports, probes, NF2FF, SAR | Implemented baseline | `src/fdtd/` |
| Full-wave frequency-domain EM | 1-D reference, 2-D Tri3 Nedelec and sparse 3-D Tet4 Nedelec driven Maxwell with PEC boundaries and edge-current excitation | Implemented baseline | `src/em/frequency_domain.cpp`, `src/em/edge_fem2d.cpp`, `src/em/edge_fem3d.cpp` |
| Electromagnetic eigenmodes | Generalized 1-D dielectric-loaded PEC cavity modes plus analytical rectangular PEC TE/TM wave-port modes with 1 W normalization; arbitrary 2-D port and 3-D cavity eigensystems remain open | Baseline | `src/em/frequency_domain.cpp`, `src/em/wave_port.cpp` |
| Integral-equation EM | Thin-wire EFIE MoM with arbitrary orientation, feeds, loads, junctions, loss and ground image model | Baseline | `src/rf/thin_wire_mom.cpp` |
| MLFMM | No fast multipole acceleration | Open | Phase 12 |
| Asymptotic SBR / physical optics | Optical ray/physical-optics modules exist, but no electromagnetic SBR/PO/RCS solver | Open | Phase 12 |
| Planar/multilayer EM | Microstrip, stripline and CPW quasi-static models; no multilayer full-wave MoM | Baseline | `src/rf/network.cpp` |
| RF/microwave/antennas | N-port networks, Touchstone, transmission lines, waveguide baseline, thin-wire MoM, arrays, matching, efficiency and optimization | Implemented baseline | `src/rf/` |
| Circuit/system EM | SPICE-class MNA, nonlinear devices, AC/transient/noise/PSS/PZ, RF N-port bridges | Implemented baseline | `src/circuit/` |
| Signal integrity | Transmission-line/N-port/Touchstone primitives plus NRZ eye height/width and Gaussian BER post-processing | Implemented baseline | `src/rf/`, `src/circuit/`, `src/em/system.cpp` |
| Power integrity | PEEC R/L/C extraction plus PDN target-impedance, decoupling selection and loaded resistive-grid IR-drop workflow | Implemented baseline | `src/rf/peec.cpp`, `src/em/system.cpp` |
| EMC/EMI | FDTD, antennas, RF networks, cable/shield coupling and waveform/probe helpers for double-exponential and damped-sine EMC stimuli; installed co-site workflows remain open | Baseline | `src/fdtd/`, `src/rf/`, `src/em/cable.cpp`, `src/em/system.cpp` |
| Cable/harness | Frequency-dependent full-matrix multiconductor RLCG propagation, matrix loads, shield transfer impedance and effective-height field coupling | Implemented baseline | `src/em/cable.cpp` |
| Thermal steady/transient | FEM transient heat, FVM thermal transport, phase change | Implemented baseline | `src/fem/heat2d.cpp`, `src/fvm/thermal_transport.cpp` |
| Conjugate heat transfer | Shared CHT coupling layer over fluid/thermal fields | Implemented baseline | `src/multiphysics/cht.cpp` |
| Structural mechanics | Linear/axisymmetric elasticity, material laws, contact, dynamics and modal baseline | Implemented baseline | `src/fem/` |
| EM -> thermal -> mechanical | Joule/electrothermal and thermoelastic coupling paths exist | Implemented baseline | `src/multiphysics/electro_thermal.cpp`, FEM elasticity/heat |
| Charged-particle tracking | Non-relativistic/relativistic Boris and Vay pushers plus absorbing/specular walls, secondary-emission yield and Monte-Carlo neutral collision baseline | Implemented baseline | `src/particle/electromagnetic.cpp` |
| Electrostatic PIC | New periodic 1-D CIC deposition + spectral Poisson + particle push | Implemented baseline | `src/particle/electromagnetic.cpp` |
| Electromagnetic PIC | Periodic 1-D/3V, 2-D/3-V and 3-D/3-V self-consistent EM-PIC baselines plus a compact staggered 3-D Yee EM-PIC foundation with Vay push, sponge/electric-wall primitives and spectral charge-continuity current reconstruction; local deposition/PML/domain decomposition remains open | Baseline | `src/particle/electromagnetic.cpp` |
| Wakefields | Causal resonator point-charge wake, bunch wake convolution and longitudinal impedance transform | Implemented baseline | `src/particle/wakefield.cpp` |
| Plasma/discharge | Secondary-emission surface-yield, Monte-Carlo elastic/ionization collisions and ionization-source coupling into EM-PIC exist; multi-species plasma chemistry, multipactor and corona remain open | Baseline | `src/particle/electromagnetic.cpp` |
| Bioelectromagnetics | SAR post-processing, SAR material helper, heterogeneous voxel tissue ingestion, local mass-averaged SAR, voxel-to-Pennes projection and implicit 2-D Pennes bioheat coupling | Implemented baseline | `src/fdtd/postprocess.cpp`, `src/multiphysics/bioheat.cpp` |
| Photonics | Sequential/non-sequential optics, diffraction/POP, Gaussian beam, polarization and dispersive glasses; dedicated full-wave photonic eigenmode/band solver open | Baseline | `src/optics/`, `src/fdtd/` |
| Optimization / UQ | Optical optimization/tolerancing, antenna optimizer, circuit Monte Carlo/sensitivity primitives | Implemented baseline | `src/optics/optimization.cpp`, `src/rf/antenna.cpp`, `src/circuit/analysis.cpp` |

## New numerical baselines in this expansion

### Driven frequency-domain Maxwell

`cfd::em::solve_pec_driven_maxwell_1d` solves the transverse 1-D phasor equation with PEC end walls, material permittivity/permeability, conductivity and impressed current density. The implementation uses a complex finite-difference Helmholtz system and reconstructs staggered magnetic field samples. Validation uses a sinusoidal manufactured current whose continuum field amplitude is known analytically.

### Electromagnetic cavity eigenmodes

`cfd::em::pec_cavity_eigenmodes_1d` assembles the generalized dielectric-loaded 1-D curl-curl analogue and solves the symmetric transformed eigenproblem. Regression checks the first three modes and the analytical PEC-cavity fundamental frequency.

### 2-D Nedelec edge-element Maxwell

`cfd::em::solve_driven_edge_maxwell_2d` adds a true curl-conforming vector-FEM baseline on triangular meshes. It uses lowest-order first-kind Nedelec basis functions, deterministic global edge orientation, complex curl-curl/mass/conduction assembly, impressed-current loading and PEC tangential constraints on boundary edges. A manufactured in-plane field is solved on two mesh resolutions and the regression requires the refinement error to decrease.

The 2-D implementation remains useful as a compact reference, while v0.9.1 extends the same edge-element ideas to sparse 3-D Tet4 Maxwell. Production gaps now move to open/impedance/radiation boundaries, arbitrary cross-section FEM port eigensystems, 3-D cavity eigenanalysis and larger-scale preconditioning.

### 3-D sparse edge-FEM, ports, Q and adaptive refinement

`cfd::em::solve_driven_edge_maxwell_3d` assembles a lowest-order Tet4 Nedelec curl-curl system directly into the shared complex sparse layer. The solver supports PEC boundary edges, volumetric impressed current and directed interior-edge excitation, reconstructs element E/H fields, and feeds resonator energy/loss/Q post-processing. Rectangular PEC TE/TM wave-port modes provide analytical cutoff and 1 W power normalization. Face-jump indicators, Dorfler marking and conforming longest-edge star bisection provide the first local 3-D RF mesh-adaptation loop.

### Cable/harness and low-frequency magnetics

`cfd::em::solve_multiconductor_cable` adds full-matrix frequency-dependent RLCG propagation with arbitrary matrix loads, shield transfer impedance and simple field-to-cable coupling. The magnetostatic baseline now supports nonlinear B-H Picard iteration, stranded-coil current density and flux-linkage/inductance export, plus Maxwell-stress boundary force/torque integration. These are still reference implementations rather than electrical-machine or EMC certification workflows.

### Relativistic particles and wakes

The particle core now includes a proper-velocity relativistic Boris update, a Vay update for high-gamma crossed-field work, absorbing/specular walls, a secondary-emission yield seam, deterministic Monte-Carlo elastic/ionization collisions and a periodic 1-D/3V self-consistent electromagnetic PIC baseline. `cfd::particle::wakefield` adds causal resonator wakes, bunch convolution and longitudinal impedance transforms. Multi-species plasma chemistry, multipactor and gas-breakdown workflows remain open.

### Charged-particle tracking and electrostatic PIC

`cfd::particle::boris_push` implements the non-relativistic Boris rotation. A uniform magnetic-field regression checks long-run speed conservation. `cfd::particle::ElectrostaticPic1D` adds periodic cloud-in-cell deposition, a spectral periodic Poisson field solve, field interpolation and particle advance. The Poisson kernel is validated against a sinusoidal charge-density field and the PIC loop against a neutralized uniform distribution.

### Self-consistent 1-D electromagnetic PIC in v0.9.3

`cfd::particle::ElectromagneticPic1D` adds the first self-consistent EM-PIC workflow: periodic Yee-style transverse `Ey/Ez/By/Bz` field updates, optional longitudinal `Ex` from periodic Poisson, Vay particle pushes from gathered fields, CIC transverse current deposition and a spectral longitudinal-current reconstruction that satisfies charge continuity from old/new particle charge density. The same step can apply the existing Monte-Carlo neutral collision/ionization model so emitted secondary macro-particles continue in later field/particle updates.

Validation covers charge-continuity residuals, source-free periodic wave energy envelope behavior, uniform transverse beam-current response through Ampere's law and MCC ionization/secondary coupling. It is still a reference 1-D/3V periodic implementation; production gaps include 2-D/3-D charge-conserving current deposition, absorbing/metallic boundaries, domain decomposition, GPU kernels and detailed plasma chemistry.

### SAR to Pennes bioheat

`cfd::multiphysics::sar_from_rms_electric_field` converts conductivity, tissue density and RMS electric field to W/kg SAR. `PennesBioheat2D` advances the Pennes equation implicitly with conduction, blood perfusion, metabolic heat and spatial SAR loading. Periodic and fixed-temperature thermal boundaries are supported. A uniform analytical implicit step and fixed-boundary regression validate the implementation.


### SI/PI/EMC workflows and heterogeneous bioheat in v0.9.2

`cfd::em::analyze_nrz_eye` adds deterministic eye/BER post-processing for known-symbol NRZ waveforms. `pdn_parallel_impedance_ohm`, `optimize_decoupling_greedy` and `solve_pdn_ir_drop_grid` provide early power-integrity target-impedance, decoupling and loaded-rail drop checks. `sample_double_exponential_pulse`, `sample_damped_sine` and `probe_waveform` add reusable ESD/BCI/lightning-style source and probe primitives.

The bioelectromagnetic layer now accepts heterogeneous voxel tissues through `VoxelTissueGrid3D`, computes material-specific voxel SAR, estimates local mass-averaged SAR and projects 3-D SAR into the 2-D implicit Pennes solver.

## Highest-value next implementation order

1. **3-D RF FEM boundaries/eigenanalysis**: impedance/open boundaries, arbitrary cross-section wave ports, cavity edge-element eigenmodes and stronger preconditioning.
2. **Surface-current MoM + MLFMM**: RWG basis first, then acceleration and hybridization with volumetric solvers.
3. **Electrical machines**: transient A-phi, nonlinear eddy-current materials, solid-conductor circuit coils and moving/sliding interfaces.
4. **Electromagnetic PIC scale-up**: replace the nodal 3-D/3-V periodic baseline with production staggered Yee fields, local charge-conserving deposition, PML/conductor boundaries, domain decomposition and accelerator kernels.
5. **Plasma/breakdown**: multi-species plasma chemistry source terms, multipactor and gas-breakdown threshold workflows.
6. **Bioelectromagnetic feedback**: temperature-dependent EM/perfusion feedback and implant/wearable exposure validation.
7. **Full-wave photonics**: dispersive waveguide modes, Bloch boundaries and band structures.

## Explicit non-parity statement

The repository now covers a wide set of the same *physics categories* as a CST-style multiphysics environment, but it is not yet CST-equivalent in solver maturity. The largest remaining gaps are production open-boundary/eigenmode RF FEM, RWG/MLFMM/asymptotic EM, full certification-grade SI/PI/EMC workflows, transient electrical-machine formulations, production 2-D/3-D EM-PIC/plasma/breakdown, closed-loop bioelectromagnetic feedback, and full-wave photonic band/mode solvers.


## v0.9.7 audit addendum

The particle/PIC branch now includes a 3-D electrostatic PIC foundation: trilinear CIC deposition, spectral 3-D Poisson fields, trilinear gather, particle stepping and a spectral 3-D charge-conserving current reconstruction. This reduces risk for later full 3-D Yee EM-PIC because the particle geometry and continuity layer are now independently validated. Remaining production CST-class gaps are still full 3-D EM-PIC staggering, PML/conductor boundaries, domain decomposition, GPU kernels, plasma chemistry and breakdown workflows.


## v0.9.8 audit addendum

The particle/PIC branch now includes a compact periodic 3-D/3-V electromagnetic PIC baseline. `ElectromagneticPic3D` stores all electric and magnetic components on a 3-D nodal grid, advances centered Maxwell curls, gathers trilinear E/B fields, pushes `PicParticle3D` macro-particles with the Vay pusher and couples the v0.9.7 spectral Jx/Jy/Jz current reconstruction into Ampere's law. This validates full 3-D self-consistent particle-field plumbing before the production staggered-Yee/PML/GPU path. Remaining CST-class gaps are still local deposition, true 3-D Yee staggering, conductor/PML boundaries, domain decomposition, accelerator kernels, plasma chemistry and breakdown workflows.


## v0.9.9 audit addendum

The particle/PIC branch now includes `StaggeredElectromagneticPic3D`, a compact true-Yee 3-D/3-V reference. It keeps Ex/Ey/Ez and Bx/By/Bz on component-specific staggered locations, advances Faraday/Ampere curl updates with forward/backward Yee differences, gathers staggered E/B fields with component offsets, and reuses the validated 3-D spectral current reconstruction for charge-consistent current coupling. Electric-wall and absorbing-sponge field-boundary primitives plus absorbing/specular particle walls give a bridge from the periodic reference toward conductor/PML workflows. Remaining production gaps are local Esirkepov/Villasenor-Buneman current deposition, CPML/PML material boundaries, parallel particle sorting, guard-cell exchange, and accelerator kernels.

## v0.10.0 audit addendum

The PIC roadmap now has a first local, finite-volume charge-conserving 3-D current reconstruction. Earlier 3-D EM-PIC releases used spectral current reconstruction to guarantee continuity; v0.10.0 adds a production-path bridge that satisfies the periodic backward-difference continuity equation using compact cumulative face fluxes. The staggered 3-D EM-PIC solver can select either the spectral or local finite-volume current path.

This is still not a full Esirkepov/Villasenor-Buneman deposition implementation. The next CST-class gap remains true local charge-conserving current deposition with higher-order particle shapes, followed by CPML/conductor material boundaries and MPI/GPU particle-field decomposition.


## v0.10.1 audit addendum

The PIC roadmap now includes the first execution-layout layer needed for large particle-field runs: stable 3-D particle sorting by cell, cell offset/count tables for cache-friendly particle loops, original-index recovery for diagnostics/restart, and guard-halo classification for particles near sub-domain faces. `StaggeredElectromagneticPic3D` can optionally sort particles at a configurable interval and reports sort passes plus occupied cells.

This is not yet MPI domain decomposition, guard-cell field exchange or GPU particle bins. It is the validated serial data-structure foundation for those next steps, and it intentionally remains separate from the current-deposition algorithm so the existing continuity regressions remain meaningful.


## v0.10.2 audit addendum

Added the first serial execution-decomposition layer for the CST-class PIC roadmap:

- Cartesian 3-D domain descriptors split global PIC grids into uneven but non-empty subdomains.
- Particle migration bucketing reports the destination subdomain for each macro-particle, wrapping periodic positions and preserving outside-particle diagnostics for nonperiodic runs.
- Scalar guard-cell exchange constructs ghost-padded local blocks by copying neighbouring subdomain layers, providing a deterministic contract for later MPI field halo exchange.

Remaining execution gaps: actual MPI particle migration, guard-cell field exchange over communicator topologies, local high-order Esirkepov/Villasenor-Buneman current deposition, CPML/conductor material boundaries, particle sorting by tile on GPU, and distributed particle-field load balancing.


## v0.10.3 audit addendum

The PIC execution roadmap now has explicit communicator-ready message contracts:

- particle migration messages split retained local particles from cross-domain payloads and preserve source-domain/source-index metadata;
- scalar guard-cell messages carry destination padded coordinates, so the communication layer can be separated from decomposition geometry;
- applying those messages reproduces the validated serial guard-cell exchange.

This reduces risk for the next MPI implementation because the communication payload format is now tested independently of MPI. Remaining execution gaps are actual MPI send/receive plumbing, vector/staggered EM field guard exchange, particle migration integrated into a distributed time step, GPU particle bins/kernels, dynamic load balancing and production CPML/conductor boundaries.


## v0.10.4 audit addendum

The PIC execution roadmap now has a deterministic transport layer between serial pack/apply helpers and future MPI. `PicRankTopology3D` maps subdomains to logical ranks, `PicTransportEnvelope3D` carries rank-addressed particle or scalar-guard payloads, and `InMemoryPicTransport3D` validates send/receive-style delivery without an MPI runtime. This closes the architecture gap between communicator-ready payload generation and actual distributed exchange, while keeping real MPI transport and GPU bins as open production steps.


## v0.10.5 addendum - communicator serialization

The PIC/plasma execution path now has a deterministic serialized message contract for rank-addressed particle migration and scalar guard-cell exchange. This is still not a completed distributed EM-PIC engine: real MPI multi-rank execution, guard-field exchange inside the time step, GPU particle bins/kernels and local high-order Esirkepov/Villasenor-Buneman current deposition remain open production work. v0.10.5 also adds initial plasma reaction-network and RF-breakdown threshold baselines, but these are threshold/chemistry utilities rather than full self-consistent discharge solvers.


## v0.10.6 audit addendum

The PIC execution roadmap now has a validated rank-local distributed exchange round. Earlier releases separately validated particle migration messages, scalar guard-cell messages, rank-addressed envelopes and serialization. v0.10.6 composes those pieces into `run_serialized_distributed_pic_exchange_round_3d`, which packs particles and scalar guards, serializes by logical destination rank, decodes and applies the results. It also adds `exchange_electromagnetic_field_guard_cells_3d`, a six-component Ex/Ey/Ez/Bx/By/Bz guard wrapper needed by the staggered 3-D EM-PIC field arrays.

Remaining CST-class execution gaps are real MPI runtime validation, integration of this exchange round into the live staggered EM-PIC step, CPML/conductor material boundaries, dynamic particle load balancing, GPU particle bins/kernels and production local high-order charge-conserving deposition.

## v0.10.7 audit addendum

The PIC execution path now includes a distributed timestep bridge rather than only standalone exchange primitives. `run_serialized_distributed_staggered_pic_step_3d` performs one rank-local orchestration pass: particles are pushed with domain-local EM fields, particle ownership is migrated through the serialized rank-addressed message contract, and six-component EM guard fields are reconstructed for subsequent local updates. This closes the validated serial timestep-ordering gap before real MPI send/receive and GPU kernels are attached.

Limitations remain explicit: the bridge is not a production multi-rank MPI run, does not yet perform local Esirkepov/Villasenor-Buneman current deposition across domain boundaries and does not overlap communication with computation.

## v0.10.8 Geant4 addendum

Geant4 was added as a clean-room reference for particle-through-matter transport.
The new baseline does not attempt Geant4 compatibility. It imports the
architecture lesson: a track/step/process separation, geometry/material-limited
step arbitration, secondaries and scoring. This complements the CST particle/PIC
path by adding a non-PIC Monte-Carlo transport seam suitable for later dose,
shielding, detector and radiation-material interaction studies.


## v0.10.9 addendum - stochastic transport, dose scoring and campaign automation

The v0.10.9 slice improves the CST-class particle/dose/workflow path without importing upstream code. It adds stochastic sampled process lengths, physics-list bundling, BVH region lookup, sensitive-detector hit collections, dose-grid scoring and SAR projection into the Pennes bioheat solver. This is useful for early particle-through-matter, detector, dose and bioheat coupling studies.

The same release also adds SU2/csauto-inspired workflow primitives: deterministic factorial/LHS DOE generation, template placeholder and IF/ENDIF rendering, solver-adapter descriptors, campaign registry summaries and finite-difference gradient utilities. These support validation campaigns and optimization sweeps around the solvers already present in the repository.

Remaining gaps: hierarchical constructive geometry, variance reduction, production particle physics data tables, full adjoint solvers, external solver launching, web dashboards, MPI/GPU campaign execution and quantitative dosimetry validation remain future work.
