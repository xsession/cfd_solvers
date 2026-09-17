# GPU backend strategy

The accelerator language remains C++20 and SYCL is the portability API. AdaptiveCpp is the first implementation target.

## Memory model

Hot solver state uses device USM (`sycl::malloc_device`) so population arrays remain resident on the selected accelerator across time steps. Host USM (`sycl::malloc_host`) is used only for the default MPI staging path and behaves like pinned/device-accessible host memory.

The distributed runtime never copies the entire local lattice merely to exchange a halo. Phase 3B first packs only boundary-crossing populations into compact device buffers.

## Distributed transports

### Staged host - portable default

Compact device send buffers are copied into pinned host buffers before MPI. Received pinned-host buffers are copied back into compact device buffers and unpacked by a device kernel. This works with an ordinary host-buffer MPI implementation.

### Direct device MPI - explicit opt-in

When the selected MPI stack is accelerator-aware, `--gpu-aware-mpi` passes the compact device-USM buffers directly to persistent MPI point-to-point requests. This removes the host staging copies, but support depends on the MPI transport and the accelerator backend behind SYCL.

Do not infer support only because the code compiled. Validate the MPI stack and backend on the target machine.

## Overlap

The distributed device time step uses separate compact send/receive buffers, allowing this order:

1. device halo pack;
2. start MPI;
3. device strict-interior collision while MPI is active;
4. finish MPI and device halo unpack;
5. device boundary-shell collision;
6. population-grid swap.

The current implementation uses an in-order SYCL queue so staging copies, unpack and boundary kernels have deterministic dependencies without accidental implementation-specific scheduling.

## Device mapping

Node-local MPI rank selects a visible accelerator ordinal deterministically. External visibility controls can still be used to map ranks to NUMA-local GPUs. Oversubscription is reported rather than silently hidden.

## Numerical policy

- D3Q19/D3Q27 distributed GPU kernels retain the conventional two-grid pull algorithm as the distributed correctness reference.
- The single-grid in-place scheme remains the local high-performance path and will be introduced into distributed execution only after the two-grid GPU/MPI parity matrix is hardware-validated.
- CPU reference solvers remain available for field-level comparisons.
- Fast-math and reduced-precision storage stay opt-in.

## Current validation gap

The development container has GCC/Clang but no AdaptiveCpp or physical accelerator runtime. The MPI+SYCL code is source/syntax validated locally, and staged distributed parity tests are compiled whenever both features are enabled, but execution still requires GPU-capable CI or target hardware.

## Device-resident execution contract (v0.15.5)

A solver is not considered GPU-resident merely because its dominant operator has a SYCL kernel. The hot numerical state must remain on the accelerator across the time loop and host transfers must be explicit boundaries.

The v0.15.5 contract is:

1. initialize or upload state before the measured time loop;
2. reset `cfd::core::DeviceTransferStats`;
3. execute one or more time steps without full-field host transfers;
4. perform reductions/derived-field operations on device where practical;
5. transfer only the data explicitly requested for host-visible API/I/O;
6. count synchronization points separately from bytes moved;
7. validate the contract on physical accelerators before claiming production qualification.

### LBM changes in v0.15.5

- Uniform and Taylor-Green populations are generated directly in device USM.
- The in-place Esoteric-Pull path reconstructs `rho`, `ux`, `uy`, `uz` into persistent device scratch.
- `download_macroscopic()` copies four macroscopic arrays rather than q population arrays.
- `total_mass()` performs a device reduction and returns one scalar instead of downloading the lattice.
- CLI/benchmark mass checks use the reduction path.
- Resident-byte and transfer-stat inspection APIs expose the current memory/transfer boundary.

For D3Q19, a macroscopic host extraction now transfers 4 values/cell instead of 19 population values/cell, a 78.9% reduction in this explicit readback. For D3Q27 the corresponding reduction is 85.2%. This comparison concerns diagnostic readback volume, not the stream-collide kernel traffic itself.

### Sparse linear algebra changes in v0.15.5

CSR matrix storage was already persistent on-device. The SYCL sparse backend now also owns persistent RHS/solution/residual/search/SpMV work vectors and reuses them between calls, removing per-call device allocation/free churn. Direct device-USM `multiply_device()` and `conjugate_gradient_device()` entry points allow resident callers to avoid full-vector host staging.

This is still a bridge rather than a fully resident FVM solve: the current CG convergence path observes shared scalar reductions on the host, and today's FVM pressure/momentum fields/operators are still CPU-owned. A future GPU FVM field/operator layer must remove those remaining boundaries from pressure/momentum coupling.

### Clean-room performance model

The public FluidX3D architecture is used only as a requirements reference. Its license is source-available non-commercial, so no source translation is used here. The independent rule carried into `cfd_solvers` is bandwidth-first: measure bytes/cell/update, keep hot fields resident, communicate only halos, reuse scratch allocations, and avoid whole-domain transfers for diagnostics.

See `FLUIDX3D_GPU_RESIDENCY_RESEARCH.md` for the detailed source review, cross-solver mapping and future GPU-R1..R4 plan.

## Resident LBM multiphysics GPU-R1 (v0.15.6)

v0.15.6 extends the device-resident contract beyond the core population lattice. `EsotericPullSyclSolver` now exposes an optional three-component device acceleration field. Thermal buoyancy, capillary forces and immersed-boundary particle reaction can accumulate into this field on the same in-order queue, and the next stream/collide kernel consumes it directly.

New device-resident components:

- `ResidentThermalD3Q7Sycl`: D3Q7 passive-scalar/temperature transport, resident scalar extraction/reduction, and Boussinesq force accumulation.
- `ResidentFreeSurfaceSycl3D`: periodic first-order conservative VOF advection, on-device interface normals/curvature, volume reduction and CSF capillary acceleration.
- `ResidentParticlesSycl3D`: SoA particle state, trilinear fluid interpolation, drag/gravity push and atomic equal-and-opposite two-way force spreading.
- `ResidentFlowDiagnosticsSycl3D`: central-difference 3-D Q-criterion generated directly from resident macroscopic velocity.

The intended no-output step order is: clear the device force field -> refresh resident macroscopic fields -> advance/couple thermal, free-surface and particles -> execute the D3Q19 step. No complete lattice, velocity, temperature, fill or particle arrays need to cross the host memory boundary inside that block. Initial particle upload, explicit downloads, checkpoints and visualization exports remain boundary operations by design.

The current implementation uses an in-order SYCL queue to preserve coupling dependencies without host waits. Particle reaction spreading uses device atomics because several particles may contribute to the same Eulerian cell. Future profiling should decide whether spatial particle binning and block-local accumulation outperform direct atomics on each backend.

Real-device numerical and performance qualification is not claimed by v0.15.6. The development environment has no physical SYCL accelerator, so the release combines default CPU/HDF5/Python regression with strict fake-SYCL compilation of the resident paths.

## Resident FDTD GPU-R2 (v0.15.7)

`cfd::fdtd::ResidentMaxwell3DSycl` is the first structured electromagnetic solver in the project that follows the same strict device-ownership contract as the resident LBM path. It owns six Yee field arrays, five per-cell material arrays, six electric CPML memory arrays, six magnetic CPML memory arrays, and compact per-axis CPML coefficient profiles in SYCL allocations.

The hot path alternates device magnetic and electric curl kernels and then enforces only the configured physical outer faces. Anisotropic permittivity, permeability and electric conductivity are applied directly from resident coefficients. For CPML boundaries, every derivative participating in a curl has its own convolution-memory term; the compact per-axis `kappa/b/c` profiles are uploaded once during setup.

Device-facing APIs cover:

- Gaussian field initialization and soft-source injection;
- device-side material-box updates;
- E/H stepping with CPML memory retained across steps;
- a scalar energy reduction without downloading the six full field arrays;
- direct probe capture into caller-owned device sample buffers;
- six-face pack/unpack into caller-owned device halo buffers;
- per-face physical-boundary masks so internal decomposition faces do not receive PEC/PMC/CPML treatment.

The halo seam is deliberately not labelled `MPI + SYCL domain decomposition`: no distributed Maxwell scheduler, neighbor exchange orchestration or overlap timing has been implemented yet. The intended next distributed step is to connect these device buffers to the existing GPU-aware MPI exchange infrastructure and then split update work into communication-independent interior and halo-dependent boundary regions.

Physical-GPU qualification is also still open. v0.15.7 is validated by the default CPU project matrix plus strict fake-SYCL compilation; real devices must still measure CPML reflection, energy evolution, field parity, bytes/update and multi-device scaling.

## Resident FVM GPU-R3 (v0.15.8)

`cfd::fvm::ResidentPolyMeshSycl` mirrors the existing general owner/neighbour finite-volume topology rather than introducing a separate Cartesian-only GPU mesh. Cell-to-face adjacency is flattened once, while cell centres/volumes, face centres/areas, interpolation weights, face-normal metrics and orthogonal pressure metrics remain in persistent device allocations.

The first resident operator set covers scalar face interpolation, Gauss gradients, vector divergence, orthogonal Laplacians, predictor face fluxes, continuity RHS generation and pressure/velocity correction. These kernels are cell-owned or face-owned and therefore avoid global atomics in the baseline.

`SyclCsrLinearAlgebra` now accepts an existing `sycl::queue`. This is a context-ownership change rather than a convenience overload: USM fields allocated by the FVM resident context can be passed directly into `conjugate_gradient_device()` without crossing into a second queue/context or being staged through host vectors.

`ResidentPressureProjectionSycl` demonstrates that seam with a constant positive mobility and a symmetric pinned pressure-reference cell. The pressure solve still synchronizes scalar Krylov reductions to make convergence decisions, but the predictor, pressure, RHS, pressure gradient, corrected velocity and face flux remain device-resident.

The production CPU SIMPLE/PISO/PIMPLE implementation is not yet replaced or claimed fully resident. GPU-R3 still needs resident momentum assembly, variable coefficient refresh, corrected/non-orthogonal operators, device-native boundary descriptors, turbulence/species/energy assembly and distributed FVM halo scheduling. Physical accelerator qualification remains open.

