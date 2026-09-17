# Source capability map

This file tracks what to reimplement, not what to copy.

## OpenFOAM 14 -> FVM and continuum physics

Target families:

- incompressible fluid: SIMPLE/PISO/PIMPLE-style pressure/velocity coupling
- compressible fluid and energy equation
- turbulence: laminar, RANS, LES
- VOF and multiphase transport
- conjugate heat transfer
- reacting/species transport
- particles and porous/media source terms
- dynamic meshes/AMR
- function-object style diagnostics and post-processing

Phase 4A adds the clean-room `PolyMesh` owner/neighbour topology, boundary patches, cell-parallel Gauss gradient/divergence, an orthogonal two-point Laplacian, shared matrix-free CG, and a staggered conservative pressure projection. Phase 4B adds linear/upwind face schemes, conservative convective divergence and transient Taylor-Green flow. Phase 4C adds collocated unstructured momentum, non-orthogonal correction, Rhie-Chow-style flux handling and SIMPLE/PISO/PIMPLE-style coupling. v0.5.0 moves collocated momentum to shared ILU(0)-GMRES and adds bounded reconstruction plus generic scalar transport. Remaining FVM targets are turbulence, compressible/energy/species flow, VOF/multiphase, moving mesh/AMR and CHT.

## FluidX3D -> LBM performance architecture

Target families:

- D3Q19/D3Q27 fluid LBM
- thermal and multiphase distribution functions
- high-throughput device-resident timestepping
- published in-place streaming patterns
- optional compact distribution storage
- heterogeneous multi-device decomposition

Phase 2 adds D3Q19/D3Q27 and an independent implementation of the published Esoteric-Pull single-grid algorithm on CPU plus optional SYCL. The original D2Q9 two-grid implementation remains as a parity oracle. No FluidX3D source is copied.

## Elmer FEM -> general multiphysics FEM

Target families:

- unstructured element/topology library
- quadrature and basis functions
- sparse CSR assembly
- CG/GMRES/BiCGStab and preconditioners
- heat, elasticity, electrostatics/electromagnetics
- coupled multiphysics and nonlinear Newton loops
- adaptivity/error estimators

v0.6.0 extends the original minimal FEM proof into a reusable reference-element/assembly layer: Line2/Tri3/Quad4/Tet4/Hex8/Prism6/Pyramid5 topology, shape functions, quadrature, isoparametric Jacobians, mixed scalar boundary conditions, assembled/matrix-free Tri3 Laplace, 2-D heat/elasticity/electrostatics/DC conduction, axisymmetric elasticity and a convergent 3-D Tet4 Poisson baseline. Nonlinear mechanics, frequency-domain EM, adaptivity and distributed FEM remain planned.

## openEMS FDTD -> electromagnetic time domain

Target families:

- 3-D Yee cells
- PEC/PMC/periodic boundaries
- Mur/CPML/PML absorbing boundaries
- anisotropic/dispersive materials
- lumped ports, waveguide excitation, probes
- near-to-far field transforms
- MPI/SYCL domain decomposition

The 1-D Yee/CFL baseline includes heterogeneous dielectric/conductive coefficients, hard/soft sources, first-order Mur absorption and reusable time/DFT monitors. v0.6.0 adds a 3-D Cartesian Maxwell PEC baseline. v0.7.0 independently adds Debye/Drude/Lorentz ADE media, a PMC enclosure baseline, TEM-like forward/backward wave ports with S-parameter extraction, and legacy VTK E/H output. v0.7.1 adds a 1-D CPML baseline, +x TFSF injection and field-coupled parallel lumped R/L/C cells. Anisotropic media, general 3-D waveguide ports/TFSF, NF2FF/SAR and distributed SYCL FDTD remain planned.

## Optiland -> optics

Target families:

- sequential/non-sequential ray tracing
- surfaces, coordinate transforms and apertures
- materials/dispersion/coatings
- paraxial and physical-optics analyses
- aberrations/MTF/spot diagrams
- optical optimization
- large batched ray tracing on SYCL

The optics layer now includes reflection/Snell primitives, sequential spherical/plane real-ray tracing, paraxial first-order tracing, Sellmeier dispersion, Jones/Stokes polarization, dielectric Fresnel coefficients and normal-incidence multilayer thin films. Non-sequential tracing, wave/physical optics, full analysis/tolerancing/optimization and SYCL batched tracing remain planned.

## Phase 3 runtime references

The distributed runtime is framework code rather than a port from one solver project.

- MPI topology/exchange semantics follow the MPI standard's Cartesian topology and nonblocking point-to-point communication model. Phase 3A intentionally uses explicit 26-neighbour messages because D3Q19/D3Q27 require edge/corner data as well as faces.
- AdaptiveCpp/SYCL device assignment uses the devices visible to the SYCL runtime. Backend/device visibility remains controlled externally by AdaptiveCpp/vendor environment settings, so the solver does not hard-code NVIDIA-, AMD-, or Intel-specific device IDs.
- No MPI or AdaptiveCpp implementation source is copied into cfd_solvers.

Phase 3A provides the CPU/MPI correctness reference. Phase 3B implements selective persistent population packing plus staged and opt-in direct-device MPI+SYCL transport; real accelerator/MPI hardware validation remains pending.


## Phase 4C clean-room pressure coupling

The Phase-4C collocated solver was independently implemented in C++20. OpenFOAM 14's public incompressible module was consulted only for architectural concepts: momentum equation -> pressure-free predictor/inverse diagonal -> pressure/non-orthogonal correction -> face flux -> reconstructed velocity. No GPL implementation text or mechanically translated source is included.


## Coupled multiphysics -> shared field/transfer architecture

The coupling layer is framework code rather than a port from one upstream. v0.7.0 introduces explicit field metadata (units/location/topology/producer), conservative remapping/flux transfer and an Aitken-accelerated partitioned fixed-point driver. The first physical chain reuses the native FEM DC-conduction and transient-heat modules to compute `J -> |J|^2/sigma -> q'` without hidden shared globals. General FEM/FVM projections and monolithic block coupling remain planned.

## RF/antenna references -> Phase 10

Palace (Apache-2.0), OpenSEMBA FDTD (MIT) and OpenNEC (MIT) are capability/validation references for full-wave FEM, EMC/FDTD and thin-wire MoM respectively. The cfd_solvers RF code is independently implemented and currently covers network/N-port math, common analytic transmission structures, antenna references, PEEC extraction and a readable thin-wire MoM baseline. GPL or unclear-license RF projects are behavioral/scientific references only.

## SPICE references -> Phase 11

ngspice is the primary permissive behavioral reference for MNA analyses and compact-model workflows. Xyce and QucsatorRF are clean-room feature/architecture references. OpenVAF/OSDI defines the planned dynamic compact-model boundary. v0.8.0 implements an independent MNA engine and parser plus OSDI library discovery; no ngspice/Xyce/QucsatorRF/OpenVAF simulator source is copied or mechanically translated.

## CST-class system workflows -> Phase 12

The CST-style expansion is clean-room framework code rather than a port of CST or any commercial solver. v0.9.0 introduced frequency-domain Maxwell references, electrostatic PIC and SAR-to-Pennes coupling. v0.9.1 added complex sparse algebra, 3-D Tet4 Nedelec Maxwell, rectangular wave-port normalization, resonator Q, adaptive RF refinement, multiconductor cable/harness propagation, relativistic Boris particles/wakes and nonlinear low-frequency magnetics. v0.9.2 adds system-level SI/PI/EMC utilities, a Vay pusher, Monte-Carlo neutral collisions and heterogeneous voxel SAR transfer into the Pennes bioheat solver.

Relevant files:

- `include/cfd/em/system.hpp`, `src/em/system.cpp`: eye/BER metrics, PDN impedance/IR-drop/decoupling and EMC waveform/probe helpers.
- `include/cfd/particle/electromagnetic.hpp`, `src/particle/electromagnetic.cpp`: Boris/Vay pushers, particle-wall interactions, secondary emission, Monte-Carlo collisions and electrostatic PIC.
- `include/cfd/multiphysics/bioheat.hpp`, `src/multiphysics/bioheat.cpp`: SAR conversion, heterogeneous voxel tissues, mass-averaged SAR, voxel-to-Pennes projection and implicit Pennes stepping.

Remaining Phase-12 gaps are still explicit in `docs/INTEGRATION_TRACKER.md`: RWG/MLFMM/SBR, transient machine A-phi, EM-PIC, plasma chemistry/breakdown, temperature-dependent bio-EM feedback and full-wave photonic band/mode solvers.
