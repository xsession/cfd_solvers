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

Started in v0.8.0. Validated baselines include generic N-port network math, Touchstone SnP, de-embedding, mixed-mode transforms, stability/gain metrics, common transmission-line/waveguide models, thin-wire dipole/array utilities, a center-fed thin-wire MoM reference, PEEC conductor extraction and circuit-port S-parameter bridges. Remaining work includes higher-fidelity MoM/NEC geometry and ground models, broadband vector fitting/model-order reduction, PEEC capacitive/proximity effects, discontinuity/modal matching, full 3-D RF FEM wave ports/eigenmodes and hardware-backed distributed/GPU RF.

## Phase 11 - SPICE-class circuits and compact models

v0.8.0 provides a clean-room MNA baseline with passive elements, controlled sources, nonlinear diode/MOS/BJT/JFET devices, switches, mutual inductance, RF N-port elements, DC/AC/transient analysis, sweeps/noise/sensitivity/Monte Carlo/Fourier, hierarchical parameterized netlists and an OSDI/OpenVAF loader seam. Next priorities are sparse MNA migration, adaptive LTE/time-step control, pole-zero/PSS/harmonic-balance, richer source/device models, full OSDI descriptor evaluation/DAE state handling, RAW/data interoperability and broader vendor-model compatibility.

## Release checkpoints

- v0.2.x: LBM performance and physical boundaries.
- v0.3.x: distributed CPU/GPU runtime baseline.
- v0.4.x: finite-volume operators and incompressible pressure-velocity coupling.
- v0.5.0: shared sparse numerics plus electrochemistry/corrosion foundation.
- v0.6.0/v0.6.1: FEM/FDTD/optics breadth, then nonlinear/adaptive FEM, Darcy, magnetostatics and modal analysis.
- v0.7.0: dispersive/port/output FDTD expansion plus first validated coupled-multiphysics infrastructure.
- v0.7.1: CPML, TFSF and field-coupled lumped R/L/C FDTD baselines.
