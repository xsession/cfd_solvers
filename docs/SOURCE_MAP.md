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

Phase 4A adds the clean-room `PolyMesh` owner/neighbour topology, boundary patches, cell-parallel Gauss gradient/divergence, an orthogonal two-point Laplacian, shared matrix-free CG, and a staggered conservative pressure projection. The Phase-4B baseline adds linear/upwind face schemes, conservative convective divergence, and a transient periodic staggered Navier-Stokes Taylor-Green solver. The next FVM milestone is collocated unstructured momentum transport, non-orthogonal correction, Rhie-Chow-style flux handling, and SIMPLE/PISO/PIMPLE-style coupling.

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

Phase 1 implements a minimal assembled linear Poisson FEM problem to establish the API.

## openEMS FDTD -> electromagnetic time domain

Target families:

- 3-D Yee cells
- PEC/PMC/periodic boundaries
- Mur/CPML/PML absorbing boundaries
- anisotropic/dispersive materials
- lumped ports, waveguide excitation, probes
- near-to-far field transforms
- MPI/SYCL domain decomposition

Phase 1 implements the 1-D Yee update and CFL-safe timestep.

## Optiland -> optics

Target families:

- sequential/non-sequential ray tracing
- surfaces, coordinate transforms and apertures
- materials/dispersion/coatings
- paraxial and physical-optics analyses
- aberrations/MTF/spot diagrams
- optical optimization
- large batched ray tracing on SYCL

Phase 1 implements vector, reflection and Snell refraction primitives.

## Phase 3 runtime references

The distributed runtime is framework code rather than a port from one solver project.

- MPI topology/exchange semantics follow the MPI standard's Cartesian topology and nonblocking point-to-point communication model. Phase 3A intentionally uses explicit 26-neighbour messages because D3Q19/D3Q27 require edge/corner data as well as faces.
- AdaptiveCpp/SYCL device assignment uses the devices visible to the SYCL runtime. Backend/device visibility remains controlled externally by AdaptiveCpp/vendor environment settings, so the solver does not hard-code NVIDIA-, AMD-, or Intel-specific device IDs.
- No MPI or AdaptiveCpp implementation source is copied into cfd_solvers.

Phase 3A provides the CPU/MPI correctness reference. Phase 3B implements selective persistent population packing plus staged and opt-in direct-device MPI+SYCL transport; real accelerator/MPI hardware validation remains pending.


## Phase 4C clean-room pressure coupling

The Phase-4C collocated solver was independently implemented in C++20. OpenFOAM 14's public incompressible module was consulted only for architectural concepts: momentum equation -> pressure-free predictor/inverse diagonal -> pressure/non-orthogonal correction -> face flux -> reconstructed velocity. No GPL implementation text or mechanically translated source is included.
