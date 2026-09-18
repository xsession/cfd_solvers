# Unified solver roadmap

The authoritative, machine-counted checklist is [`INTEGRATION_TRACKER.md`](INTEGRATION_TRACKER.md). This file gives the same phase numbering at a higher level. A capability is only checked in the tracker after code and a deterministic validation path exist.

Current development status is produced with:

```bash
python3 scripts/integration_status.py
python3 scripts/integration_status.py --json
```

## Phase 0 - project rules and traceability

Build/release provenance, license boundaries, pinned upstream references, compiler/sanitizer gates, reproducible source archives, patches and Git bundles. Baseline complete.

## Phase 1 - common HPC runtime

Shared aligned/SoA storage, OpenMP, optional SYCL/MPI, device assignment, distributed decomposition, checkpoint/restart and common sparse/Krylov algebra. Next priorities are NUMA/SIMD abstractions, memory pools, distributed SpMV and accelerator sparse kernels.

## Phase 2 - FluidX3D-class LBM

D2Q9/D3Q19/D3Q27, two-grid reference and single-grid in-place streaming, forcing, physical boundaries, CPU/SYCL parity, multi-rank compact halo exchange and scaling tools are implemented baselines. Remaining work includes advanced collision models, thermal/passive-scalar/free-surface physics, curved/immersed geometry, particles, voxelization and hardware-specific mixed/compressed precision.

## Phase 3 - OpenFOAM-class finite volume

Owner/neighbour `PolyMesh`, Gauss operators, bounded reconstruction, transient and collocated incompressible flow, Rhie-Chow, SIMPLE/PISO/PIMPLE, non-orthogonal correction, scalar transport and shared Krylov momentum solves are in place. Remaining groups include runtime-selectable schemes, higher-order temporal integration, turbulence, compressible flow, VOF/multiphase, reacting/species/energy equations, particles, moving mesh/AMR and production I/O.

## Phase 4 - Elmer-class finite element

Reference elements and quadrature, sparse/matrix-free assembly, scalar diffusion/heat, electrostatics/DC conduction, linear/axisymmetric elasticity, Tet4 Poisson, Newton nonlinear solve, Tri3 adaptivity, Darcy flow, 2-D magnetostatics and modal eigenanalysis are implemented baselines. Remaining work includes mixed/vector elements, nonlinear mechanics/contact, Brinkman/incompressible FEM, harmonic/edge-element electromagnetics, multiphysics couplings and distributed/GPU FEM.

## Phase 5 - openEMS-class FDTD

1-D and 3-D Yee baselines, heterogeneous/lossy media, PEC/PMC, Mur ABC, hard/soft sources, probes/DFT, Debye/Drude/Lorentz ADE media, TEM-like wave ports/S-parameters and VTK E/H output are implemented. Next priorities are CPML, TFSF, lumped RLC, cylindrical/multigrid coordinates, anisotropy, NF2FF, SAR, HDF5 and SYCL/MPI production kernels.

## Phase 6 - Optiland-class optics

Sequential spherical/plane ray tracing, paraxial analysis, Sellmeier materials, polarization, Fresnel interfaces and multilayer coatings are implemented baselines. Next groups are transforms/decenters, conics/aspheres/freeforms, non-sequential/ghost tracing, physical optics, PSF/MTF/wavefront/Zernike analysis, optimization/tolerancing, import/export and accelerator ray batches.

## Phase 7 - chemistry, electrochemistry and corrosion

Species/Arrhenius kinetics, Nernst/Butler-Volmer/Faraday relations, 1-D and `PolyMesh` Nernst-Planck, electroneutral and Poisson potentials, Scharfetter-Gummel fluxes, reacting electrode boundaries and galvanic mixed-potential calculations are implemented baselines. Next priorities are reversible/stiff chemistry, thermodynamic/activity models, aqueous equilibrium/speciation, porous electrodes, implicit electrode/transport coupling, passivation/pitting/product layers, moving/phase-field corrosion interfaces, stress coupling and EIS.

## Phase 8 - coupled multiphysics

Started in v0.7.0. Implemented baselines now include a unit/location/topology-aware field registry, conservative cell and face transfers, fixed-point partitioned coupling, Aitken relaxation, and a shared-mesh DC-conduction -> Joule-heating -> transient-FEM-thermal chain. Next are FEM<->FVM and structured<->unstructured projections, monolithic block coupling, thermal-structural-optical chains, electrochemistry-to-thermal coupling, corrosion-driven mesh motion and porous electrochemical flow/thermal cases.

## Phase 9 - interoperability, workflow and optimization infrastructure

The tracker, upstream pin manifest, benchmark history and reproducible release tooling already exist. Remaining work includes common case/material schemas, Gmsh/VTK/OpenFOAM mesh import, STL/OBJ geometry, VTK/HDF5/XDMF output, Python bindings, parameter/optimization runners, restartable workflow graphs and additional performance work such as reusable sparse assembly patterns, SIMD, mixed precision and backend autotuning.

## Phase 10 - RF, antennas and microwave networks

Started in v0.8.0. Validated baselines include generic N-port network math, Touchstone SnP, de-embedding, mixed-mode transforms, stability/gain metrics, common transmission-line/waveguide models, thin-wire dipole/array utilities, center-fed, coupled parallel-wire and arbitrary-orientation disjoint-wire MoM references, PEEC conductor extraction and circuit-port S-parameter bridges. Remaining work includes higher-fidelity MoM/NEC arbitrary connected geometry, junction basis functions, ground/image models and canonical convergence validation, broadband vector fitting/model-order reduction, PEEC capacitive/proximity effects, discontinuity/modal matching, full 3-D RF FEM wave ports/eigenmodes and hardware-backed distributed/GPU RF.

## Phase 11 - SPICE-class circuits and compact models

v0.8.0 provides a clean-room MNA baseline with passive elements, controlled sources, nonlinear diode/MOS/BJT/JFET devices, switches, mutual inductance, RF N-port elements, DC/AC/transient analysis, sweeps/noise/sensitivity/Monte Carlo/Fourier, hierarchical parameterized netlists and an OSDI/OpenVAF loader seam. v0.8.1 adds adaptive step-doubling LTE control, a fitted small-signal pole/zero baseline and periodic steady-state shooting-by-settling. v0.8.2 adds unequal-step adaptive BDF2/Gear control, an ideal-transformer convenience model, and validated lossy sampled transmission-line behavior. Next priorities are sparse MNA migration, BDF3+ order selection, direct generalized-eigenvalue pole-zero analysis, harmonic balance, richer source/device models, full OSDI descriptor evaluation/DAE state handling, RAW/data interoperability and broader vendor-model compatibility.


## Phase 13 - multibody, rigid bodies and granular dynamics

v0.16.0 starts a Project-Chrono-class clean-room mechanics family: 6-DOF quaternion rigid bodies, Jacobian/PGS constraints, distance/spherical/revolute/prismatic/fixed/gear joints, sweep-and-prune plus BVH sphere collision, smooth penalty and nonsmooth impulse contact, and explicit spherical DEM with Hertz-type contact, friction, rolling resistance and cohesion. v0.16.1 adds motors, GJK/EPA, triangle-mesh contact, persistent manifolds and Mindlin history; v0.16.2-v0.16.6 add resident GPU DEM and deterministic distributed contact/history/force-exchange infrastructure; v0.16.7 adds conservative resolved and unresolved CFD/DEM coupling. v0.16.8 adds nonlinear implicit Newmark/generalized-alpha/HHT integration, a coupled fixed-base reduced-coordinate articulated-tree solver, and conservative floating-frame FEM-interface kinematic/force transfer. Phase 13 is therefore algorithmically complete except for the intentionally open distributed-memory DEM qualification item, which still requires a real multi-rank GPU-aware MPI run.


## Planned solver-family expansion after Phase 13

The uncovered-domain audit is retained in `UNCOVERED_SOLVER_FIELDS_RESEARCH.md`. These phases are intentionally **not** added to the machine-counted integration denominator until implementation begins. Current planned order:

- Phase 14: acoustics, ultrasound and photoacoustics (k-Wave-class pseudospectral methods).
- Phase 15: semiconductor TCAD (DEVSIM-class Poisson/carrier drift-diffusion).
- Phase 16A: battery cell physics (PyBaMM-class SPM/SPMe/DFN and degradation).
- Phase 16B: nuclear neutronics (OpenMC-class Monte Carlo and depletion).
- Phase 17: electrical power systems (pandapower/VeraGrid-class load flow, OPF and dynamics).
- Phase 18: geomechanics and porous THMC (OpenGeoSys/PyLith-class coupled subsurface physics).
- Phase 19: atomistics and molecular dynamics (HOOMD-blue/LAMMPS-class particle mechanics).
- Phase 20: explicit impact/crash mechanics (OpenRadioss-class high-rate nonlinear FEM).
- Phase 21: electronic structure/DFT (Quantum ESPRESSO-class plane-wave methods).

Horizontal common-runtime references from the audit—SUNDIALS, hypre and AMReX—should be mined continuously rather than treated as separate physics phases.

## Release checkpoints

- v0.2.x: LBM performance and physical boundaries.
- v0.3.x: distributed CPU/GPU runtime baseline.
- v0.4.x: finite-volume operators and incompressible pressure-velocity coupling.
- v0.5.0: shared sparse numerics plus electrochemistry/corrosion foundation.
- v0.6.0/v0.6.1: FEM/FDTD/optics breadth, then nonlinear/adaptive FEM, Darcy, magnetostatics and modal analysis.
- v0.7.0: dispersive/port/output FDTD expansion plus first validated coupled-multiphysics infrastructure.
- v0.7.1: CPML, TFSF and field-coupled lumped R/L/C FDTD baselines.
- v0.8.0: RF/microwave networks, thin-wire/PEEC references and the first SPICE-class MNA engine.
- v0.8.1: adaptive circuit LTE control, pole-zero/PSS baselines and coupled parallel-wire MoM.
- v0.8.2: unequal-step adaptive BDF2/Gear plus arbitrary-orientation disjoint-wire MoM.
- v0.16.0: Phase-13 multibody/contact/DEM foundation.
- v0.16.1: Phase-13 contact/collision maturity (motors, GJK/EPA, triangle mesh, persistent manifolds, Mindlin history).
- v0.16.2: resident SYCL DEM neighbor/contact hot loop.
- v0.16.3: bonded-particle damage/fracture plus distributed DEM ownership/migration/ghost-exchange foundations.
- v0.16.4: distributed DEM contact ownership, persistent-history migration, reverse-force exchange and MPI timestep driver.


## Phase 14 - acoustics, ultrasound and photoacoustics [complete]

v0.17.0 introduced the project first dedicated Fourier-pseudospectral wave backend: radix-2 FFTs, spectral derivatives/fractional Laplacians, a 2-D first-order k-space acoustic solver, heterogeneous media, nonlinear B/A response, split-field PML, power-law attenuation, photoacoustic initial pressure, transducer/sensor arrays, beamforming and acoustic post-processing. v0.17.1 completes the phase with Dirichlet time reversal plus structural, CFD, Pennes bioheat, piezoelectric FEM and optimization/inverse-problem coupling baselines.

## Battery continuation v0.19.1

Series SPMe pack control and dissipative balancing now have deterministic validation. Remaining Phase 16A work is full DFN/P2D and 3-D battery thermal coupling.

## Battery continuation v0.19.2

3-D pack thermal fields and conservative SPMe electrothermal coupling are validated. Phase 16A now has one remaining item: full DFN/P2D.
