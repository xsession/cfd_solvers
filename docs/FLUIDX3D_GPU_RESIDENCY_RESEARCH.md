# FluidX3D GPU-residency research and clean-room mapping

Research refresh: 2026-09-17

## Scope and clean-room rule

This note studies public architectural ideas visible in FluidX3D documentation, its public OpenCL wrapper, published papers and public issue reports. It is not a source-code port. FluidX3D explicitly describes its license as **source-available, no-cost, non-commercial**, not open source, and does not permit commercial use under its current license. `cfd_solvers` therefore uses the public material only to derive independent requirements and performance patterns; implementation in this repository remains independently written.

Primary public references:

- FluidX3D repository and README: <https://github.com/ProjectPhysX/FluidX3D>
- FluidX3D OpenCL wrapper: <https://github.com/ProjectPhysX/FluidX3D/blob/master/src/opencl.hpp>
- FluidX3D documentation: <https://github.com/ProjectPhysX/FluidX3D/blob/master/DOCUMENTATION.md>
- M. Lehmann, *Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs*, Computation 10(6), 92 (2022): <https://doi.org/10.3390/computation10060092>
- Public portability failure example, Apple OpenCL compiler issue #349: <https://github.com/ProjectPhysX/FluidX3D/issues/349>

## The coding mindset

### 1. The accelerator owns hot state

A high-performance GPU solver should not be structured as a CPU solver that periodically calls an accelerator. The dominant state for the time loop should be allocated once on the device and remain there across many time steps. Host memory is a control, I/O and interoperability boundary, not the canonical copy of every evolving field.

FluidX3D's OpenCL memory wrapper makes host and device allocations independent, including the option to create a device allocation without a host buffer. Transfers are explicit and can cover subranges rather than the whole allocation. This is the key architectural lesson to preserve across solver families.

**`cfd_solvers` contract:** after initialization, a device solver should be able to execute `N` time steps with zero explicit host-transfer bytes. A host transfer is allowed only when an API call explicitly requests host-visible data, checkpoint/output requires it, or a non-device-aware external library forms a deliberate boundary.

### 2. Optimize bytes per cell per time step before FLOPs

FluidX3D treats the LBM update as memory-bandwidth bound. Its published/public performance analysis reports D3Q19 traffic in bytes per cell per time step and uses in-place streaming and reduced storage precision to lower the memory wall. The Esoteric Pull/Push paper shows why eliminating a second population lattice is valuable even when arithmetic count does not shrink.

For bandwidth-dominated kernels, the first questions should therefore be:

1. How many bytes are read and written per cell/update?
2. Are the accesses coalesced/contiguous for adjacent work items?
3. Can a temporary field be recomputed more cheaply than stored?
4. Can two passes be fused without excessive register pressure?
5. Can storage precision be reduced with an explicit accuracy budget?

Kernel-count reduction by itself is not the goal. Fusion is useful only when it eliminates meaningful memory traffic or launch/synchronization cost without destroying occupancy, cache behavior or maintainability.

### 3. Perform initialization, reductions and derived fields on the device

A device-resident time loop can still be undermined by apparently small convenience operations:

- constructing a full initial population lattice on the CPU and uploading it;
- downloading all populations merely to compute density/velocity;
- downloading a full force field to sum a scalar diagnostic;
- downloading volume fields just to render them.

These operations should instead become device kernels/reductions, with only the requested result transferred. Scalar reductions ideally return one scalar; visualization should consume device fields directly whenever possible.

### 4. Treat synchronization as a performance cost

A tiny transfer can still serialize a queue. The runtime should track synchronization points separately from transfer bytes. Convergence checks in an iterative solver can therefore be expensive even if only one residual scalar reaches the CPU per iteration.

Longer-term pressure/FEM Krylov work should use device-side convergence or bounded batches of iterations when practical, with explicit accuracy/termination semantics.

### 5. Persistent workspaces beat repeated allocation

Sparse solvers, particle systems and multiphysics operators should allocate their hot scratch space once and reuse it. Allocation is both overhead and an implicit lifetime/synchronization hazard. A persistent matrix with per-call GPU vector allocations is only partially resident.

### 6. Communicate halos, not domains

Multi-GPU execution should exchange only boundary-crossing data. Pack/unpack should run on the accelerator. With GPU-aware MPI or another validated device transport, device buffers should be passed directly; otherwise pinned-host staging is an explicit portability fallback. Interior computation should overlap communication where dependencies permit.

### 7. Unified-memory zero-copy is an optimization, not an assumption

FluidX3D can fuse host/device storage on CPUs/iGPUs using host-backed OpenCL buffers and imposes alignment requirements for efficient transfers. This is attractive on unified-memory systems, but driver behavior matters. `cfd_solvers` should keep device-USM as the conservative discrete-GPU baseline and add zero-copy/shared-USM policies only behind runtime capability checks and hardware regression tests.

### 8. Hardware validation is part of numerical correctness

Public FluidX3D history contains several vendor/driver workarounds, and a 2026 Apple OpenCL issue reports a compiler/runtime combination that silently produced incorrect physics for particular feature combinations. Accelerator qualification therefore cannot stop at "compiled successfully." Each supported backend/device family needs physical numerical parity tests.

`cfd_solvers` must keep NVIDIA/AMD/Intel/other hardware qualification checkboxes open until real hardware runners execute the corresponding validation matrix.

## Mapping onto cfd_solvers

### LBM: closest to full GPU residency

LBM is the best fit because collision/streaming is explicit and local/stencil based.

Already present before v0.15.5:

- persistent SYCL population storage;
- in-place Esoteric-Pull baseline;
- D2Q9/D3Q19/D3Q27 GPU kernels;
- compact distributed halo exchange;
- device-aware/direct-device MPI seam;
- optional reduced population storage and autotuning infrastructure.

v0.15.5 closes several accidental host round-trips:

- common uniform and Taylor-Green initial conditions are generated directly on device;
- macroscopic `rho/u` are reconstructed into persistent device scratch;
- host macroscopic extraction transfers only four fields instead of the entire q-population lattice;
- mass is reduced on-device and only the final scalar is observed by the host;
- transfer/synchronization accounting provides an executable residency contract.

Remaining work for a genuinely feature-complete GPU-resident LBM path:

- port thermal/free-surface extension updates to SYCL-resident state;
- keep particle integration, force spreading and compaction entirely on device;
- make Q-criterion, interface geometry and visualization consume resident fields directly;
- qualify multi-GPU parity on physical NVIDIA/AMD/Intel hardware.

### FDTD: very suitable

FDTD is also naturally GPU-resident because E/H updates are explicit structured stencils. The target architecture is resident E/H/material/PML arrays, alternating update kernels, on-device monitors/reductions and compact halo exchange. Port/PML/postprocessing data should cross to the host only when explicitly requested.

### PIC and Lagrangian particles: suitable with data-oriented redesign

Particle state should be SoA and resident. Push, field gather, deposition, boundary handling, sorting/binning, compaction and collision candidate generation should happen on device. Avoid rebuilding host-side particle lists each step. Output should sample or compact the requested subset before transfer.

### Finite volume: possible, but much more work

A fully GPU-resident FVM loop is feasible but requires more than a GPU SpMV:

1. device representation/mirror of `PolyMesh` connectivity and geometry;
2. resident `GeometricField`-style state;
3. device face interpolation, gradient/divergence and flux kernels;
4. device matrix assembly or matrix-free pressure/momentum operators;
5. persistent device linear-solver state and preconditioners;
6. pressure-velocity correction without downloading full fields;
7. device-side residual/norm decisions or deliberately batched host checks;
8. GPU-resident boundary/source-model execution.

v0.15.5 improves the sparse-solver bridge by reusing persistent Krylov work buffers and exposing direct device-USM SpMV/CG entry points. The direct CG path keeps volume vectors on-device but still synchronizes scalar reductions for convergence decisions. It does **not** claim that the current FVM SIMPLE/PISO/PIMPLE loop is fully GPU-resident.

### FEM: possible with matrix-free/element kernels

Element integration, residual/Jacobian-vector application and iterative linear algebra map well to GPUs, especially matrix-free formulations. Sparse global assembly is harder and should use coloring/atomics or element-by-element operators deliberately. Current FEM GPU-element tracker items remain future work.

### Chemistry/reacting flow: good for batching, difficult for divergence

Independent per-cell thermochemistry is highly parallel. GPU batching can be effective for equation-of-state evaluation, source evaluation and groups of similarly stiff reactors. Highly stiff adaptive ODE integration introduces branch divergence and variable work, so solver selection/binning matters more than simply offloading each cell.

### Optics/ray workloads: naturally parallel

Rays, wavelengths and many analysis samples are independent enough for GPU batching. Device-resident material/geometry acceleration structures should be shared by repeated analyses, with only final images/statistics transferred.

### Workflow and I/O should remain host-side

There is no benefit in forcing orchestration, filesystem access, case parsing, controller/API logic or HDF5 metadata management onto a GPU. "Whole simulation on GPU" means the **hot numerical loop and derived diagnostics** remain device-resident; it does not mean eliminating the CPU from the application.

## v0.15.5 executable residency contract

`cfd::core::DeviceTransferStats` records:

- host -> device bytes;
- device -> host bytes;
- device -> device explicit bytes;
- explicit synchronization points.

A device solver can reset the counters immediately after setup, execute a time-step block, and assert that `host_transfer_bytes() == 0`. Host extraction then has a measurable transfer budget.

The new focused LBM regression uses exactly this pattern. It intentionally treats a scalar mass reduction as a synchronization point but not as a bulk lattice transfer, because the reduction output uses shared scalar storage in the current SYCL implementation.

## Concrete next implementation phases

### GPU-R1 - resident LBM extensions

Implemented in v0.15.6:

- [x] D3Q7 thermal/passive-scalar populations and resident temperature extraction;
- [x] free-surface fill transport plus resident normals/curvature/capillary forcing;
- [x] immersed-boundary particle state, interpolation and atomic two-way reaction spreading;
- [x] on-device Q-criterion;
- [x] shared persistent device acceleration field consumed directly by D3Q19;
- [ ] device-native visualization/export pipeline for all new fields;
- [ ] spatial particle binning/local accumulation after hardware profiling;
- [ ] physical NVIDIA/AMD/Intel regression and bandwidth profiling.

### GPU-R2 - structured FDTD residency

Implemented in v0.15.7:

- [x] persistent 3-D E/H fields;
- [x] persistent anisotropic/lossy material coefficients;
- [x] persistent CPML coefficient profiles and twelve convolution-memory fields;
- [x] alternating device E/H update kernels;
- [x] scalar device reduction and probe capture into device buffers;
- [x] six-face device halo pack/unpack plus physical-face masks;
- [ ] interior/boundary kernel split and communication overlap;
- [ ] complete MPI + SYCL Maxwell domain-decomposition driver;
- [ ] kernel fusion only if hardware profiling shows a benefit;
- [ ] physical NVIDIA/AMD/Intel CPML/parity/bandwidth qualification.

### GPU-R3 - FVM mesh/field mirror

- [x] immutable general `PolyMesh` topology/geometry mirror;
- [x] flattened resident cell/face connectivity and persistent pressure/projection work fields;
- [x] device gradient/divergence/orthogonal-laplacian/interpolation kernels;
- [x] direct queue/context bridge into device CSR/CG;
- [x] resident constant-mobility pressure projection with pressure-reference handling;
- [x] resident first-order momentum coefficient/source assembly;
- [x] direct device-USM BiCGStab for nonsymmetric momentum/scalar CSR systems;
- [x] resident fixed/zero-gradient/slip velocity descriptors and non-orthogonal pressure correction;
- [x] variable-coefficient pressure matrix numeric refresh on-device;
- [x] resident first-order SIMPLE/PISO/PIMPLE baseline;
- [x] fixed-pressure/outlet pressure patch baseline;
- [x] resident generic implicit scalar advection-diffusion seam for temperature/species fields;
- [ ] full CPU boundary-condition and higher-order scheme parity;
- [x] shared resident field registry plus SA/k-epsilon/SST turbulence transport;
- [x] resident temperature/species scalar equations with one-step Arrhenius heat/species source coupling;
- [ ] full wall-function, thermophysical-property and detailed chemistry-mechanism parity inside the coupled loop;
- [ ] real-device parity/performance qualification.

### GPU-R4 - resident multiphysics task graph

- express solver coupling as data dependencies between device kernels;
- remove unnecessary `wait()` calls;
- optionally use command graphs only after backend support is mature and qualified;
- expose transfer/synchronization counters in benchmark JSON.

## Acceptance criteria

A solver family can be called **GPU-resident** only when all of the following hold for its qualified feature set:

1. setup transfers are completed before the measured time loop;
2. `N` time steps execute without full-field host round-trips;
3. mandatory residual/diagnostic traffic is scalar or intentionally compact;
4. checkpoint/output transfers happen only at requested output boundaries;
5. host-transfer byte counters remain zero inside a no-output time-loop regression;
6. physical hardware parity passes on every backend/device family claimed as supported;
7. memory capacity and bytes/update are reported alongside throughput.

This definition is intentionally stricter than "there is a GPU kernel" and gives the project a measurable route toward FluidX3D-class accelerator efficiency without importing FluidX3D implementation code.

## v0.15.7 GPU-R2 result

The structured FDTD family now has a concrete resident accelerator implementation rather than only a roadmap. `ResidentMaxwell3DSycl` keeps the Yee fields, material coefficients and CPML state on the accelerator and exposes device-native probes and halo buffers. In contrast to a per-operation offload design, an ordinary no-output timestep does not require a field-volume CPU transfer.

The remaining distributed gap is scheduling, not data ownership: the new physical-face mask allows a local brick to distinguish physical boundaries from MPI interfaces, while pack/unpack writes directly to device buffers. The next MPI step should reuse the existing GPU-aware communication layer, overlap interior E/H work with exchange, and update boundary slabs after halo completion.
## v0.15.8 GPU-R3 result

The unstructured/polyhedral FVM family now has a concrete device-ownership layer. `ResidentPolyMeshSycl` preserves the project's existing owner/neighbour face representation and uploads only a flattened adjacency/geometry view needed by GPU kernels. This avoids creating a second solver-specific mesh model and lets CPU and accelerator paths share the same case/mesh semantics.

The resident pressure path also removes a subtle USM ownership barrier: `SyclCsrLinearAlgebra` can reuse the FVM queue/context, so pressure/RHS vectors do not need a host or cross-context copy before device CG. A constant-mobility projection baseline demonstrates predictor flux, continuity RHS, pressure solve and flux/velocity correction as one no-bulk-host-transfer chain.

This is not yet the end state described by the FluidX3D-inspired residency rule. Momentum assembly, turbulence/species/energy equations, variable pressure coefficients, general boundary conditions, non-orthogonal corrections and MPI halo overlap still need GPU-native implementations. The release therefore calls the result a GPU-R3 baseline, not a fully resident production Navier-Stokes solver.


## v0.15.9 GPU-R3.1 result

The FVM residency work now crosses the control-loop boundary. Momentum coefficients and sources are formed from resident face fluxes and pressure gradients, `V/aP` mobility is refreshed on-device, and the fixed pressure sparsity pattern receives new numeric coefficients through a device-to-device update. Explicit non-orthogonal pressure fluxes are included in the pressure RHS and final corrected face flux.

`ResidentIncompressibleSycl` wires those pieces into SIMPLE, PISO and PIMPLE baselines while retaining the existing CPU `CollocatedIncompressible` implementation untouched. The hot-loop residency target is therefore demonstrated architecturally, but the resident path is deliberately narrower numerically: Euler time stepping, upwind convection, Jacobi momentum solution and pinned all-Neumann pressure are the current qualified design surface. Full boundary/turbulence/species/energy parity and real accelerator testing are the next requirements before this can be called a production-equivalent GPU Navier-Stokes path.


## v0.15.10 GPU-R3.2 result

The resident FVM stack now distinguishes the linear-algebra character of its equations instead of forcing every block through the same iteration. Pressure retains the SPD CG path; first-order upwind momentum and scalar transport use a new persistent device-USM BiCGStab implementation. Momentum CSR numeric values are regenerated on-device from resident face fluxes while connectivity remains fixed.

Fixed-pressure patches are also supported without a host sparsity rebuild. If no fixed pressure is present, the pressure operator uses the symmetric pinned reference-cell formulation. If a fixed-pressure patch exists, the reference pin is disabled and Dirichlet face coefficients/RHS terms are assembled on-device. Boundary pressure values are included in the Gauss pressure gradient and non-orthogonal correction path.

`ResidentScalarTransportSycl` adds the first reusable non-momentum transported equation to GPU-R3. It accepts a resident face-flux pointer, fixed/zero-gradient scalar boundary conditions and an optional device volumetric source. This is the execution seam needed to move temperature, passive species and later turbulence variables into the same no-bulk-host-transfer timestep. The standalone class currently mirrors mesh geometry separately; a future field-registry/task-graph layer should share one resident mesh allocation across all equations.


## v0.15.11 GPU-R3.3 result

The FVM accelerator path now shares one mesh/queue/context across flow, turbulence and coupled scalar equations. `ResidentFvmFieldRegistrySycl` owns only the evolving cell fields, while `ResidentScalarEquationSycl` supplies a common variable-diffusivity/source/sink transport operator. This removes the per-equation mesh duplication noted in v0.15.10 and gives SA, k-epsilon, SST, temperature and species the same residency contract.

The reacting seam intentionally starts with a compact one-step Arrhenius source kernel: rate, fuel sink and heat source are produced on-device and copied device-to-device into the resident scalar source fields. Detailed CPU chemistry mechanisms, adaptive stiff ODE integration, wall functions and higher-order convection remain separate follow-on work rather than being hidden behind a generic 'GPU chemistry' label.


## v0.16.2 GPU DEM residency result

Phase 13 now has a resident spherical-DEM path that follows the same ownership rule as the LBM/FDTD/FVM accelerator work: evolving particle state remains device-owned through neighbor construction, contact evaluation and integration. A uniform grid is formed entirely on-device by atomic cell counts, a device prefix pass and compact bucket placement. Pair kernels inspect only 27 neighboring cells and process each pair once (`j > i`).

Persistent Mindlin state is also kept on-device. Instead of requiring a globally concurrent pair hash, each lower-index particle owns a bounded history table keyed by its higher-index contact neighbors. This provides deterministic tangential spring memory without a host reconciliation pass. The current baseline is intentionally spherical/single-device and uses single-precision particle/contact storage; distributed halos, bonded fracture, non-spherical GPU collision and physical backend qualification remain open.
