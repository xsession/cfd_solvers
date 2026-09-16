# Roadmap

The authoritative, checkable feature-completion matrix is now [`INTEGRATION_TRACKER.md`](INTEGRATION_TRACKER.md).
The chemistry/corrosion research basis is [`CHEMISTRY_CORROSION_RESEARCH.md`](CHEMISTRY_CORROSION_RESEARCH.md), and cross-cutting performance work is tracked in [`OPTIMIZATION_PLAN.md`](OPTIMIZATION_PLAN.md).

This file keeps the chronological implementation checkpoints.

The current unreleased CPU continuation is documented in [`CONTINUATION_CPU.md`](CONTINUATION_CPU.md).
It advances LBM, FVM, FEM/coupling, FDTD, optics and chemistry together; the unchecked entries below remain open work.

## Phase 1 - unified core and proof solvers [complete baseline]

- [x] C++20 build and shared namespace/layout
- [x] OpenMP multicore backend
- [x] AdaptiveCpp/SYCL optional backend wiring
- [x] aligned SoA storage
- [x] D2Q9 CPU LBM
- [x] D2Q9 SYCL baseline kernel
- [x] finite-volume diffusion
- [x] finite-element Poisson
- [x] 1-D Maxwell FDTD
- [x] optical reflection/refraction
- [x] CPU regression tests and GCC/Clang CI

## Phase 2A - LBM performance foundation [implemented]

- [x] D3Q19 descriptor and solver
- [x] D3Q27 descriptor and solver
- [x] single-population-grid Esoteric-Pull-style streaming from the published algorithm
- [x] two-grid One-Step-Pull reference implementation
- [x] odd/even-layout D2Q9/D3Q19/D3Q27 numerical parity tests
- [x] reusable benchmark harness: MLUPS, estimated DDF GB/s, memory and mass drift
- [x] SYCL single-grid D2Q9/D3Q19/D3Q27 implementation
- [x] CPU/SYCL parity tests compiled when SYCL is enabled
- [x] GCC/OpenMP, Clang, serial and ASan+UBSan CPU verification
- [ ] accelerator CI runner that compiles and executes the SYCL path

## Phase 2B - physical LBM cases [implemented on CPU reference path]

- [x] body-force formulation
- [x] halfway bounce-back walls
- [x] moving-wall boundary
- [x] inlet/outlet boundary set
- [x] lid-driven cavity regression
- [x] Poiseuille/channel manufactured regression
- [x] Taylor-Green convergence/decay study
- [x] FP32 versus FP64 reference comparison
- [ ] optional compressed population storage with explicit error gates (moved to the next performance pass)

## Phase 3A - distributed CPU/hybrid runtime [implemented baseline]

- [x] MPI wrapper and Cartesian/node-local topology discovery
- [x] aspect-aware 2-D/3-D brick decomposition with uneven extents
- [x] 26-neighbour face/edge/corner halo representation
- [x] nonblocking MPI halo exchange baseline
- [x] MPI + OpenMP hybrid execution
- [x] overlap strict-interior kernels with halo transfer
- [x] virtual-rank backend and D3Q19/D3Q27 monolithic-parity tests
- [x] rank-local binary checkpoint/restart
- [x] local-rank -> visible-SYCL-device assignment
- [x] four-rank OpenMPI CI job
- [x] local four-rank runtime execution with MPICH 4.0 in WSL (2026-09-16)

## Phase 3B - distributed accelerator/runtime optimization [implemented baseline]

- [x] device-resident distributed D3Q19/D3Q27 block kernels
- [x] one-GPU-per-rank execution using the Phase-3A local-rank mapping
- [x] pinned-host halo staging for non-GPU-aware MPI
- [x] direct device-buffer exchange for GPU-aware MPI
- [x] pack only populations that cross each process boundary
- [x] persistent communication buffers/requests and exact communication-volume profiling
- [x] CSV strong/weak-scaling benchmark output with messages/bytes per step
- [ ] real accelerator CI/runtime validation on NVIDIA, AMD and Intel devices
- [ ] multiple local devices per rank after one-GPU-per-rank is hardware-validated
- [ ] decomposition-independent checkpoint manifest / portable case format

## Phase 4A - OpenFOAM-class FVM foundation [implemented baseline]

- [x] owner/neighbour polyhedral cell/face primitives and boundary patches
- [x] validated Cartesian-hexa mesh generator
- [x] cell-to-face adjacency and oriented area-vector convention
- [x] OpenMP cell-parallel Gauss scalar gradient
- [x] OpenMP cell-parallel Gauss vector divergence
- [x] orthogonal two-point scalar Laplacian
- [x] shared matrix-free conjugate-gradient solver
- [x] staggered conservative pressure projection with periodic face fluxes
- [x] manufactured gradient/divergence/Laplacian regressions

## Phase 4B/4C - incompressible Navier-Stokes [collocated baseline implemented]

- [x] linear and first-order upwind face interpolation
- [x] conservative scalar/vector convective divergence on `PolyMesh`
- [x] transient periodic staggered momentum + pressure projection correctness solver
- [x] Taylor-Green analytical error/refinement validation
- [x] bounded higher-order convection limiter baseline
- [x] transient collocated momentum equation on arbitrary `PolyMesh`
- [x] collocated flux path with Rhie-Chow-style pressure/velocity coupling
- [x] non-orthogonal/skewness correction
- [x] SIMPLE algorithm
- [x] PISO algorithm
- [x] PIMPLE-style outer coupling
- [x] collocated cavity/channel validation cases
- [x] periodic staggered Taylor-Green FVM validation case

## Phase 4D - advanced continuum physics

- turbulence models
- dynamic mesh and AMR
- VOF/multiphase
- compressible flow
- heat/species/reacting flow
- conjugate heat transfer

## Phase 5 - Elmer-class FEM [advanced baseline in v0.6.1]

- [x] reference element catalogue and Gaussian quadrature
- [x] isoparametric mapping/Jacobians
- [x] sparse CSR and matrix-free Tri3 operator APIs
- [x] 2-D Poisson/heat/linear elasticity
- [x] mixed Dirichlet/Neumann/Robin scalar boundaries
- [x] electrostatics and DC conduction
- [x] axisymmetric linear elasticity
- [x] convergent 3-D Tet4 Poisson baseline
- [x] generic Newton nonlinear solve + nonlinear scalar FEM baseline
- [x] residual/jump error estimator + conforming Tri3 h-adaptivity baseline
- [x] saturated Darcy porous-flow baseline
- [x] 2-D magnetostatic vector-potential baseline
- [x] generalized sparse eigen/modal baseline
- [ ] nonlinear elasticity/material laws
- [ ] Brinkman/mixed incompressible FEM
- [ ] frequency-domain EM
- [ ] distributed FEM assembly/solve

## Phase 6 - openEMS-class FDTD [baseline expanded in v0.6.0]

- [x] 1-D Yee/CFL baseline
- [x] heterogeneous dielectric/conductive material coefficients
- [x] hard/soft sources
- [x] first-order Mur absorption
- [x] time/DFT monitors
- [x] 3-D Cartesian Maxwell PEC baseline
- [x] 1-D polynomial matched electric/magnetic absorbing layer
- [x] 1-D Debye/Drude/Lorentz ADE with analytical response/refinement gates
- [x] 1-D wave decomposition and legacy VTK export baseline
- [ ] 3-D SYCL production update
- [ ] multidimensional UPML/CPML
- [ ] 3-D dispersive/anisotropic materials
- [ ] lumped/waveguide ports and TFSF
- [ ] near-to-far/SAR
- [ ] multi-GPU decomposition

## Phase 7 - Optiland-class optics [baseline expanded in v0.6.0]

- [x] sequential spherical/plane real-ray tracing and apertures
- [x] paraxial first-order tracing
- [x] Sellmeier glass/material dispersion baseline
- [x] Jones/Stokes polarization primitives
- [x] Fresnel dielectric interfaces
- [x] normal-incidence multilayer thin films
- [x] conic/even-asphere sag and sequential intersection baseline
- [x] Gaussian-beam ABCD propagation and thin-lens focusing
- [ ] non-sequential/ghost/stray-light tracing
- [ ] PSF/MTF/wavefront/Zernike analyses
- [ ] physical optics propagation
- [ ] tolerancing and optimization
- [ ] SYCL batched ray/surface kernels

## Phase 8 - coupled multiphysics

Shared coupling graph for combinations such as EM heating -> thermal -> structural deformation -> fluid/optical response. Field transfer is explicit and testable rather than hidden global state.

- [x] field metadata/units, conservative 1-D cell remap and PolyMesh face-flux transfer baseline
- [x] fixed-point driver with Aitken relaxation
- [x] DC conduction -> Joule heat -> transient heat -> small-strain thermal expansion on a shared Tri3 mesh
- [ ] general coupling graph, cross-mesh FEM/FVM projection and optical deformation

## Phase 4D - numerics/performance foundation [started]

- [x] bounded limited-linear unstructured FVM reconstruction
- [x] CSR sparse matrix builder/storage
- [x] Jacobi-preconditioned CG
- [x] Jacobi-preconditioned BiCGStab
- [x] replace collocated momentum Jacobi sweeps with assembled ILU(0)-GMRES solve
- [x] reusable scalar FVM gradient/reconstruction/divergence workspaces
- [x] scalar-transport operator/ILU cache with invalidation
- [x] least-squares gradients and minmod/van-Leer reconstruction baseline
- [ ] workspace adoption throughout collocated pressure/momentum loops and additional NVD schemes

## Phase 7A - chemistry/electrochemistry/corrosion foundation [started]

- [x] species + elementary Arrhenius/mass-action kinetics
- [x] reversible reactions, van't Hoff equilibrium constants and implicit isothermal reaction integration
- [x] Nernst, Butler-Volmer and Faraday helpers
- [x] conservative 1-D Nernst-Planck transport
- [x] nonlinear 1-D ohmic-electrolyte corrosion cell
- [x] multi-species PolyMesh Nernst-Planck
- [x] electroneutral and Poisson potential solves
- [x] Scharfetter-Gummel electrochemical fluxes
- [x] Butler-Volmer/Faradaic PolyMesh electrode flux boundary
- [x] galvanic mixed-potential solver
- [ ] fully implicit nonlinear electrode/transport coupling
- [ ] aqueous speciation/equilibrium and pH coupling
- [ ] galvanic multi-electrode corrosion case
- [x] ideal dilute acid/base equilibrium and pH calculation (standalone baseline)
- [ ] moving-interface/phase-field corrosion
