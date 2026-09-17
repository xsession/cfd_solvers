# Unified Solver Integration Tracker

Last research refresh: 2026-09-17

This is the authoritative completion checklist for `cfd_solvers`. A feature is not considered complete merely because code exists. Unless explicitly marked as an exploratory baseline, completion requires:

1. a documented public API or case format;
2. a deterministic correctness/regression test;
3. at least one physical or manufactured validation case;
4. CPU execution and profiling;
5. accelerator execution when the feature is expected to run on GPUs;
6. distributed-memory validation when the feature is expected to scale across nodes;
7. documentation of assumptions and unsupported regimes;
8. no hidden dependency on upstream source code whose license is incompatible with this repository's clean-room strategy.

Status notation uses normal Markdown checkboxes so completion can be counted automatically.

## Source families being unified

- OpenFOAM 14: finite-volume CFD, mesh/case/runtime selection, multiphase, turbulence, reacting flow, particles, moving mesh, AMR, post-processing.
- FluidX3D: bandwidth-oriented LBM, cross-vendor accelerator execution, compact/in-place populations, multi-device domain decomposition, geometry voxelization, particles/free-surface/thermal extensions, force/torque and visualization.
- Elmer FEM: general multiphysics FEM spanning structural mechanics, heat transfer, fluid dynamics, electromagnetics, adaptivity and massively parallel execution.
- openEMS FDTD (`thliebig/openEMS`): 3-D Cartesian/cylindrical FDTD, absorbing boundaries, dispersive media, ports/excitations, NF2FF, SAR and parallel execution.
- Optiland: sequential/non-sequential optical design, paraxial/real ray tracing, polarization, coatings, optical analysis, optimization, tolerancing, materials, import/export and differentiable backends.
- Electrochemistry/reactive chemistry references: EchemFEM, echemAMR, Cantera, Reaktoro, PHREEQC/PhreeqcRM, and published corrosion-specific Poisson-Nernst-Planck work.
- RF/antenna references: Palace, OpenSEMBA FDTD, OpenNEC, plus clean-room behavioral reference to openEMS/Meep/scuff-em.
- SPICE/RF-circuit references: ngspice, Xyce, QucsatorRF/Qucs-S, OpenVAF/OSDI and GNUCAP.

## Phase 0 - project rules and traceability

- [x] C++20 primary implementation language.
- [x] CMake build.
- [x] OpenMP CPU backend.
- [x] optional SYCL accelerator backend.
- [x] optional MPI distributed backend.
- [x] GCC and Clang regression builds.
- [x] sanitizer smoke configuration.
- [x] third-party/license research document.
- [x] upstream feature tracker.
- [x] machine-readable upstream source manifest with commit pin and license fingerprint.
- [x] script that reports checklist completion by upstream family and phase.
- [x] benchmark result JSON schema and historical performance dashboard.
- [x] reproducible release workflow that creates source archive, patch, bundle and checksums in CI.

## Phase 1 - common HPC runtime [baseline complete]

### Data and execution
- [x] aligned allocations.
- [x] structure-of-arrays storage.
- [x] OpenMP parallel loops/reductions.
- [x] CPU thread control.
- [x] SYCL queue/device path.
- [x] rank-to-device assignment.
- [x] Cartesian MPI decomposition.
- [x] nonblocking halo exchange.
- [x] selective LBM halo exchange.
- [x] persistent MPI request path.
- [x] pinned-host and direct-device halo modes.
- [x] checkpoint/restart baseline.
- [x] asynchronous task graph shared by FVM/FEM/FDTD/LBM.
- [x] OS first-touch NUMA placement baseline for caller-owned CPU buffers.
- [x] explicit NUMA node/thread affinity and placement policy.
- [x] SIMD abstraction for small fixed-size kernels.
- [x] device-memory pool / scratch allocator.
- [x] runtime backend selection without recompiling the application.
- [x] profiling ranges/counters for every solver family.

### Shared linear algebra
- [x] matrix-free conjugate gradient.
- [x] CSR sparse matrix storage and assembly builder.
- [x] Jacobi preconditioner.
- [x] preconditioned conjugate gradient.
- [x] BiCGStab for nonsymmetric systems.
- [x] restarted GMRES.
- [x] complex sparse systems through real-block CSR + ILU(0)-GMRES with sparse fallback.
- [x] block CSR / vector-valued matrices.
- [x] ILU(0) preconditioner.
- [x] algebraic multigrid preconditioner adapter with externally supplied hierarchy/cycle callback.
- [x] geometric multigrid for structured grids.
- [x] SYCL SpMV and Krylov kernels.
- [x] distributed local-row CSR with explicit halo-column values and halo-aware SpMV baseline.
- [x] MPI-owned distributed sparse-vector exchange integrated directly into Krylov iterations.
- [x] mixed-precision iterative refinement.

**v0.15.6 accelerator-residency policy:** hot device state should remain accelerator-owned across the time loop; initialization/derived fields/reductions should execute on-device where practical; host transfers are explicit API/I/O boundaries and are measurable through `DeviceTransferStats`. The SYCL LBM path now extends the no-bulk-host-transfer contract through resident D3Q7 thermal transport, free-surface VOF/capillary coupling, immersed-boundary particle reaction forces and Q-criterion diagnostics. SYCL sparse linear algebra reuses persistent Krylov work buffers, but full FVM/FEM loops are not yet claimed device-resident. Physical GPU qualification remains open until real hardware runners execute the parity matrix.

**v0.15.7 FDTD residency:** the 3-D Yee accelerator path now keeps E/H, material coefficients and CPML memory on-device; probes and halo packing can target device buffers directly. This closes the implementation-level `SYCL FDTD kernels` item. `MPI + SYCL domain decomposition` remains open because v0.15.7 supplies only the device halo/physical-face seam, not a complete distributed Maxwell driver or hardware-qualified GPU-aware MPI execution.

**v0.15.8 FVM residency:** `ResidentPolyMeshSycl` now mirrors general owner/neighbour `PolyMesh` geometry/connectivity, and resident interpolation/gradient/divergence/Laplacian plus pressure-correction kernels can share a queue/context with direct-device CSR/CG. A pinned-Neumann resident pressure-projection baseline keeps its volume fields on-device. The optimization gate `SYCL FVM/FEM sparse algebra` remains open because resident FEM assembly/operators and physical-GPU qualification are still missing.

**v0.15.10 resident FVM GPU-R3.2:** the FVM accelerator path now assembles nonsymmetric momentum CSR values on-device, solves them through direct device-USM BiCGStab, supports fixed-pressure/outlet pressure patches, and provides a resident implicit scalar advection-diffusion equation that can consume the flow face flux directly. This still does not close `SYCL FVM/FEM sparse algebra`: full turbulence/thermophysical/reacting parity, FEM device assembly/operators, distributed FVM execution and physical-GPU qualification remain open.

## Phase 2 - FluidX3D-class LBM [advanced baseline]

### Core numerics
- [x] D2Q9 BGK.
- [x] D3Q19 BGK.
- [x] D3Q27 BGK.
- [x] two-grid reference streaming.
- [x] single-grid in-place/pull streaming.
- [x] Guo forcing.
- [x] FP32/FP64 validation.
- [x] MRT collision.
- [x] TRT collision.
- [x] regularized collision option.
- [x] cumulant collision option.
- [x] Smagorinsky/subgrid LES option.
- [x] thermal DDF lattice.
- [x] passive scalar DDF lattice.
- [x] multiphase/free-surface extension.

### Boundaries and geometry
- [x] halfway bounce-back.
- [x] moving wall.
- [x] velocity inlet.
- [x] pressure/density outlet.
- [x] channel and cavity regressions.
- [x] interpolated curved-wall bounce-back.
- [x] immersed-boundary particles.
- [x] two-way particle coupling.
- [x] STL/triangle-mesh voxelization.
- [x] GPU voxelization.
- [x] moving/rotating geometry re-voxelization.
- [x] force/torque integration on solids.
- [x] porous-media drag models.

### Performance/portability
- [x] SYCL D2Q9/D3Q19/D3Q27 baseline.
- [x] multi-rank decomposition.
- [x] selective halo payloads.
- [ ] hardware CI on NVIDIA.
- [ ] hardware CI on AMD.
- [ ] hardware CI on Intel GPU.
- [x] compressed population storage with explicit conservation/error gates.
- [x] FP16/BF16/custom packed storage experiments.
- [x] kernel-fusion/autotuning database by device.
- [x] multiple GPUs per rank.
- [x] adaptive domain repartitioning.

### Visualization/output concepts worth carrying forward
- [x] legacy VTK structured-grid export for 2-D LBM scalar/vector fields.
- [x] LBM vorticity derived-field baseline.
- [x] LBM Q-criterion derived field.
- [x] 2-D LBM streamline extraction baseline.
- [x] LBM horizontal slice extraction baseline.
- [x] headless PGM scalar-field image output baseline.
- [x] optional lightweight interactive viewer kept separate from solver core.

## Phase 3 - OpenFOAM-class finite-volume framework [in progress]

### Mesh/fields/operator framework
- [x] owner/neighbour polyhedral face representation.
- [x] boundary patches.
- [x] oriented face-area convention.
- [x] cell-to-face adjacency.
- [x] Cartesian-hexa test meshes.
- [x] deterministic skew/sheared mesh.
- [x] Gauss scalar gradient.
- [x] vector divergence.
- [x] orthogonal scalar Laplacian.
- [x] linear interpolation.
- [x] first-order upwind interpolation.
- [x] bounded limited-linear MUSCL/Barth-Jespersen reconstruction.
- [x] least-squares gradient.
- [x] corrected/snGrad operator.
- [x] deferred non-orthogonal Laplacian.
- [x] generic dimensioned field classes.
- [x] generic finite-volume matrix assembly.
- [x] run-time selectable discretization schemes.
- [x] reusable old-time field history plus validated BDF2 derivative/Crank-Nicolson update primitives.
- [x] second-order temporal schemes integrated into production FVM equation solvers.

### Pressure-velocity coupling
- [x] staggered projection reference.
- [x] collocated momentum baseline.
- [x] Rhie-Chow-style face flux.
- [x] non-orthogonal pressure correction.
- [x] SIMPLE.
- [x] PISO.
- [x] PIMPLE-style outer loops.
- [x] channel/cavity/skew regressions.
- [x] momentum solve through shared CSR/Krylov stack instead of Jacobi sweeps.
- [x] pressure AMG/multigrid path.
- [x] adaptive timestep from Courant number.
- [x] production pressure-velocity solver integration of second-order backward/CN time stepping.

### Transport and thermo
- [x] scalar advection-diffusion equation.
- [x] energy/enthalpy equation.
- [x] temperature-dependent material properties.
- [x] multi-species mass transport.
- [x] buoyancy/Boussinesq.
- [x] compressible equation-of-state framework.
- [x] compressible pressure-energy coupling.
- [x] MUSCL-minmod + Rusanov shock-capturing compressible baseline.
- [x] characteristic high-order/WENO compressible schemes.

### Turbulence
- [x] laminar/turbulence runtime interface.
- [x] mixing-length reference.
- [x] Spalart-Allmaras eddy-viscosity constitutive helper.
- [x] Spalart-Allmaras transport equation.
- [x] k-epsilon eddy-viscosity constitutive helper.
- [x] k-epsilon transport equations.
- [x] k-omega/SST eddy-viscosity constitutive helper.
- [x] k-omega/SST transport equations.
- [x] Reynolds-stress Boussinesq reconstruction baseline.
- [x] Reynolds-stress transport model framework.
- [x] LES filters and Smagorinsky/WALE.
- [x] wall functions and y+ diagnostics.
- [x] DES length-scale switching baseline.
- [x] full DES/hybrid RANS-LES transport framework.

### Multiphase/interface physics
- [x] VOF volume fraction transport.
- [x] bounded interface compression.
- [x] surface tension/CSF.
- [x] contact angle.
- [x] enthalpy/mushy-zone phase-change material baseline.
- [ ] compressible VOF.
- [x] Euler-Euler interphase-drag source baseline.
- [ ] coupled Euler-Euler phase transport framework.
- [x] drift-flux model.
- [x] lubrication thin-film flux baseline.
- [ ] coupled thin-film transport solver.

### Lagrangian/dispersed phase
- [x] particle cloud container.
- [x] particle tracking through polyhedral cells.
- [x] drag force baseline for Lagrangian particles.
- [x] lift/virtual-mass force models.
- [x] one-way/two-way particle momentum coupling baseline.
- [x] binary collision/coalescence/breakup baselines.
- [ ] many-particle collision/coalescence/breakup framework.
- [x] D2-law droplet evaporation baseline.
- [ ] coupled spray injection/evaporation solver.
- [x] cell particle-volume-fraction accumulation baseline.
- [ ] dense-particle rheology/coupled momentum framework.

### Reacting CFD
- [x] chemistry mechanism abstraction shared with Phase 8.
- [x] finite-rate species source coupling.
- [x] ideal-gas combustion thermo-state baseline with mixture cp/cv/gamma/density/enthalpy inversion.
- [x] bounded 1-D premixed reaction-diffusion laminar-flame front baseline.
- [ ] premixed/non-premixed model interfaces.
- [ ] soot/radiation coupling hooks.

### Heat/radiation/CHT
- [x] solid heat conduction.
- [x] conjugate fluid-solid heat transfer.
- [x] gray two-surface exchange baseline.
- [x] general diffuse-gray surface-to-surface radiosity/view-factor network baseline.
- [x] optically-thin participating-media source baseline.
- [x] participating-media radiation source-model interface with optically-thin implementation.
- [x] radiation/chemistry energy coupling.

### Mesh motion/adaptation
- [x] mesh-motion field.
- [x] ALE flux correction.
- [x] rigid-body motion.
- [x] topology change interface.
- [x] local refinement/coarsening.
- [x] error-indicator-driven AMR.
- [x] conservative field remap after topology changes.

### Case/IO/post-processing
- [x] unit-aware native YAML/JSON case-parameter parsing baseline.
- [x] OpenFOAM mesh importer.
- [x] OpenFOAM field importer/exporter where licensing/interoperability permits.
- [x] VTK/VTU output.
- [x] HDF5 checkpoint/field output.
- [x] probes/sampling.
- [x] forces/coefficients.
- [x] residual-history and solver function-object callback framework baseline.
- [x] derived fields: vorticity, Q, Lambda2, enstrophy, y+.

## Phase 4 - Elmer-class finite-element framework [in progress]

### FEM core
- [x] 1-D linear Poisson proof solver.
- [x] 3-D Tet4 Poisson manufactured-solution and mesh-convergence baseline.
- [x] shared CSR/Krylov foundation.
- [x] Tri3 unstructured topology + structured triangle-mesh generator baseline.
- [x] P1 affine shape-gradient/Jacobian evaluation.
- [x] one-point triangle quadrature baseline.
- [x] component/node Dirichlet constraint baseline.
- [x] global sparse element assembly into shared CSR.
- [x] element topology catalogue: line/tri/quad/tet/hex/prism/pyramid.
- [x] reference-element shape functions.
- [x] Gaussian quadrature.
- [x] isoparametric mapping/Jacobians.
- [x] Dirichlet/Neumann/Robin BC framework.
- [x] global sparse assembly.
- [x] matrix-free element operator API.
- [ ] mixed/vector finite elements.
- [x] nonlinear Newton solve with backtracking line search and shared ILU(0)-GMRES Jacobian solves.
- [x] consistent P1 transient mass matrix + implicit-Euler time integration baseline.
- [x] residual/flux-jump error estimator, Dorfler marking and conforming longest-edge Tri3 refinement.

### Structural mechanics
- [x] 2-D Tri3 linear elasticity baseline.
- [x] plane stress/strain constitutive matrices.
- [x] axisymmetric elasticity.
- [x] Green-Lagrange finite-strain truss geometry baseline.
- [ ] general nonlinear-geometry FEM assembly.
- [x] hyperelasticity.
- [x] plasticity/material-law interface.
- [x] contact.
- [x] generalized sparse eigenfrequency baseline with consistent mass and fixed-free bar validation.
- [x] structural dynamics.

### Thermal/fluid/porous
- [x] 2-D transient heat equation baseline (steady Poisson/diffusion covered separately).
- [x] phase change/latent heat.
- [x] porous heat/mass transport.
- [x] penalty FEM incompressible Stokes reference.
- [ ] mixed velocity-pressure Navier-Stokes FEM.
- [x] saturated Darcy porous-flow baseline.
- [x] Brinkman porous flow.

### Electromagnetics
- [x] electrostatics.
- [x] DC conduction.
- [x] 2-D magnetostatics through out-of-plane vector potential A_z and reconstructed B field.
- [x] eddy-current harmonic formulation.
- [x] 2-D driven frequency-domain Maxwell baseline with lowest-order Nedelec edge elements.
- [ ] 3-D frequency-domain Maxwell with curl-conforming vector elements.
- [x] electro-thermal coupling.

### Multiphysics/HPC
- [x] generic sparse monolithic block linear-system assembly and Krylov solve baseline.
- [ ] physics-specific mixed FEM monolithic block formulations.
- [x] partitioned coupling.
- [ ] MPI mesh partitioning.
- [ ] distributed sparse assembly.
- [ ] GPU element kernels.
- [ ] restart/checkpoint and adaptive repartition.

## Phase 5 - openEMS-class FDTD [in progress]

- [x] 1-D Yee proof solver.
- [x] 3-D Cartesian Yee grid baseline.
- [x] axisymmetric cylindrical TM(r,z) Maxwell FDTD baseline with explicit r=0 update.
- [ ] cylindrical multi-grid refinement.
- [x] material coefficient preprocessing.
- [x] PEC enclosure boundary baseline.
- [x] PMC boundary baseline on the 3-D Cartesian engine.
- [x] Mur absorbing boundary.
- [x] 1-D CPML absorbing-boundary baseline with graded conductivity/kappa and convolution memory.
- [x] 1-D +x total-field/scattered-field (TFSF) incident-wave injection baseline.
- [x] soft electric-field source baseline.
- [x] hard field source.
- [x] field-coupled parallel lumped R/L/C cell element baseline.
- [x] Drude ADE dispersion baseline.
- [x] Lorentz ADE dispersion baseline.
- [x] Debye ADE dispersion baseline.
- [x] anisotropic materials.
- [x] 1-D TEM-like wave-port decomposition and complex S-parameter baseline.
- [x] probes/time signals/DFT monitors.
- [x] NF2FF transform.
- [x] SAR 1g/10g calculation.
- [ ] HDF5 output.
- [x] legacy ASCII VTK E/H field output baseline.
- [x] SYCL FDTD kernels.
- [ ] MPI + SYCL domain decomposition.
- [x] geometry/material bridge shared with optics/FEM.

## Phase 6 - Optiland-class optics [in progress]

### System/modeling
- [x] vector reflection/refraction proof functions.
- [x] sequential optical system/surface list baseline.
- [x] coordinate transforms/decenters/tilts.
- [x] spherical surfaces.
- [x] conics/aspheres.
- [x] polynomial/freeform surfaces.
- [x] circular surface apertures/vignetting baseline.
- [x] explicit aperture stops and general aperture shapes.
- [x] wavelength and field models.
- [x] glass/material dispersion database.

### Ray tracing and physics
- [x] paraxial trace.
- [x] sequential real-ray trace with CPU batch execution.
- [x] polarization/Jones/Stokes representation.
- [x] birefringence.
- [x] Fresnel coatings.
- [x] multilayer coating transfer matrices.
- [x] non-sequential ray tracing.
- [x] ghost/multi-sequence tracing.
- [x] absorption/scatter.
- [x] physical-optics propagation / diffraction.
- [x] Gaussian beam propagation.

### Analysis
- [x] first-order/paraxial properties.
- [x] RMS spot-radius analysis baseline.
- [x] spot-diagram sampling/reporting.
- [x] ray fans.
- [x] wavefront error.
- [x] Zernike decomposition.
- [x] PSF.
- [x] MTF.
- [x] distortion/chromatic analysis.
- [x] non-sequential stray-light power-accounting baseline.
- [x] non-sequential multi-bounce ghost-path enumeration/ranking by residual power baseline.

### Optimization/tolerance/interoperability
- [x] variable/operand merit function.
- [x] local optimizers.
- [x] global optimizers.
- [x] automatic differentiation path.
- [x] Monte Carlo tolerancing.
- [x] sensitivity analysis.
- [x] material/glass search.
- [ ] Zemax ZMX import/export.
- [ ] CODE V SEQ import/export.
- [ ] OSLO LEN import/export.
- [x] JSON native format.
- [ ] SYCL batched ray-surface kernels.

## Phase 7 - chemistry, electrochemistry and corrosion [started]

Research basis is documented in `CHEMISTRY_CORROSION_RESEARCH.md`.

### Chemistry/thermodynamics inspired by Cantera/Reaktoro/PHREEQC
- [x] species metadata: name, molar mass, charge.
- [x] elementary mass-action reaction network.
- [x] Arrhenius forward rates.
- [x] reversible reactions/equilibrium constants.
- [x] third-body/falloff reactions.
- [x] general reaction orders.
- [x] thermodynamic phase/property model interface (ideal-gas/ideal-solution chemical potentials).
- [x] NASA-7 polynomial thermochemistry.
- [x] NIST-form Shomate cp/h/s polynomial thermochemistry.
- [x] ideal gas/ideal solution phases.
- [x] ideal dilute acid/base speciation with charge-balanced pH.
- [x] non-ideal activity model interface.
- [x] aqueous activity models.
- [x] binary 1:1 Pitzer mean-activity baseline.
- [ ] full multicomponent Pitzer interaction model and database parameterization.
- [x] SIT model.
- [x] chemical-potential/Gibbs minimization equilibrium.
- [x] kinetic ODE integrator for stiff reaction networks.
- [x] multiphase equilibrium/precipitation/dissolution.
- [x] surface complexation.
- [x] ion exchange.
- [x] thermochemical database import layer.
- [ ] PHREEQC database reader or optional PhreeqcRM adapter.
- [ ] Cantera YAML mechanism importer or optional Cantera adapter.
- [ ] Reaktoro dynamic adapter kept separate from clean core.
- [x] reaction-path analysis.
- [x] reactor-network abstractions.

### Electrochemical transport inspired by EchemFEM/echemAMR
- [x] Nernst equation helper.
- [x] Butler-Volmer current-density helper.
- [x] Faraday dissolution/penetration helper.
- [x] conservative 1-D Nernst-Planck diffusion/advection/electromigration baseline.
- [x] nonlinear 1-D ohmic-electrolyte + Butler-Volmer corrosion cell.
- [x] multi-species Nernst-Planck on PolyMesh.
- [x] electroneutral potential formulation.
- [x] Poisson-Nernst-Planck potential formulation.
- [x] Scharfetter-Gummel/exponential-fitting flux option.
- [ ] SUPG/DG FEM transport option.
- [x] porous-electrode effective transport.
- [x] electronic solid-phase conduction.
- [x] smooth signed-distance immersed electrode/interface weighting baseline.
- [ ] cut-cell/embedded-boundary conservative electrode interface.
- [x] Butler-Volmer interfacial species-flux boundary on PolyMesh.
- [x] implicit Butler-Volmer + mass-transfer boundary coupling baseline.
- [ ] fully implicit spatial PNP/electrode monolithic coupling.
- [x] classical Marcus electron-transfer kinetics baseline.
- [x] numerical Marcus-Hush-Chidsey integral kinetics baseline.
- [x] double-layer capacitance.
- [x] Bikerman finite-size activity correction baseline.
- [ ] generalized modified-PNP transport solver.
- [ ] GPU kernels for electrochemical transport.
- [ ] MPI decomposition and AMR.

### Corrosion-specific physics
- [x] galvanic/mixed-potential solver with multiple reactions/metals.
- [x] anodic metal dissolution species source via Faradaic Butler-Volmer boundary flux.
- [x] cathodic oxygen/hydrogen reaction models.
- [x] concentration-dependent Nernst equilibrium potentials.
- [x] passive-film/oxide growth state variable.
- [x] passivation/transpassive kinetics.
- [x] precipitation/product-layer coupling.
- [x] pH and aqueous speciation coupling.
- [x] localized/pitting initiation model.
- [x] level-set recession update baseline.
- [ ] moving-interface transport/remesh solver.
- [x] local Allen-Cahn phase-field corrosion update baseline.
- [x] bounded 1-D spatial phase-field corrosion PDE baseline.
- [ ] multidimensional phase-field corrosion coupled to transport.
- [x] stress-assisted exchange-current mechanochemical coupling baseline.
- [ ] fully coupled stress-corrosion transport/solid mechanics.
- [x] Faradaic geometry recession-distance helper.
- [ ] Faradaic recession coupled to remeshing.
- [x] cathodic-protection/anode models.
- [x] electrochemical impedance spectroscopy small-signal solver.
- [ ] corrosion benchmarks against published PNP/FEM cases.

## Phase 8 - coupled multiphysics [started]

- [x] explicit field registry with units/location/topology metadata.
- [x] conservative cell-to-cell transfer baseline using exact-overlap 1-D remap.
- [x] conservative face-to-cell transfer on `PolyMesh` integrated fluxes.
- [x] FEM<->FVM projection/interpolation.
- [x] structured-grid<->unstructured transfer.
- [x] fixed-point partitioned coupler.
- [x] Aitken relaxation.
- [x] monolithic sparse block coupling interface and linear solve baseline.
- [x] DC conduction -> Joule heating -> transient FEM thermal coupling on a shared Tri3 mesh.
- [x] thermal -> structural deformation.
- [x] axial structural-deformation -> sequential optical-surface vertex transfer baseline.
- [ ] CFD -> thermal -> structural -> optics chain.
- [x] electrochemistry -> heat generation -> thermal.
- [ ] corrosion recession -> mesh motion -> CFD/structural update.
- [ ] battery/electrolyzer porous electrochemistry + thermal + flow.

## Phase 9 - interoperability, workflow and UX [complete]

- [x] common case schema with units and validation.
- [x] material database schema shared by CFD/FEM/FDTD/optics/electrochemistry.
- [x] mesh import: Gmsh.
- [x] mesh import: VTK/VTU.
- [x] mesh import: OpenFOAM polyMesh.
- [x] geometry import: STL/OBJ.
- [x] output: VTK/VTU.
- [x] output: HDF5/XDMF.
- [x] Python bindings.
- [x] parameter sweep runner.
- [x] optimization/inverse-problem runner.
- [x] restartable workflow graph.
- [x] provenance metadata in every result file.
- [x] benchmark/validation catalog CLI.
- [x] factorial and Latin-hypercube DOE case generation with deterministic seed support.
- [x] campaign template placeholder plus IF/ENDIF conditional renderer.
- [x] solver-adapter descriptor and capability boundary for workflow front ends.
- [x] campaign registry summary and finite-difference gradient utility for design loops.
- [x] persistent campaign folder writer with rendered case files, DOE row snapshots, template manifest and registry roundtrip.
- [x] external solver command-plan boundary for native, Docker, Apptainer/Singularity and Slurm-style runtime wrappers.
- [x] residual/performance log parser interfaces plus a gradient-descent campaign optimization loop.
- [x] runtime doctor checks for native, container and scheduler command boundaries.
- [x] safe native argv process launcher with captured stdout/stderr logs.
- [x] local campaign runner with output discovery, residual/performance parsing and registry updates.
- [x] live-control directive files for stop/extend/checkpoint/flush through the solver-adapter capability boundary.
- [x] Slurm-style scheduler queue parser and registry status application for queued/running/done/failed/cancelled jobs.
- [x] campaign output refresh that scans logs/residual/performance files and updates registry objective/iteration/status without relaunching cases.
- [x] multi-server campaign placement planner with server capacity, online filtering, tag filtering and deterministic assignment.
- [x] generated SSH/rsync/native/Docker remote command scripts for distributed campaign launch planning.
- [x] Docker/Compose deployment scaffold with worker replicas, manager service, environment, resource limits and server inventory template.
- [x] supervised multi-server execution plan with remote health checks, launch/status/cancel/fetch-log command scripts and job metadata.
- [x] remote job status parser with registry updates for running/done/failed/stopped distributed cases.
- [x] Docker deployment healthcheck and worker entrypoint scripts for containerized campaign workers.
- [x] multi-server supervision access probes for SSH/local execution, Docker daemon availability and campaign-root writability.
- [x] retry-wrapped remote launch plans plus stdout/stderr log-tail command generation.
- [x] dashboard-ready multi-server JSON/TSV summaries with sensitive command-display redaction.
- [x] dependency-free multi-server controller API scaffold with health/status/cases/log-tail routes.
- [x] token-gated controller mutation route that writes durable campaign control directives.
- [x] generated controller OpenAPI, status JSON, routes table, env example, README and launch script artifacts.
- [x] generated browser dashboard assets with polling summary cards, case table, log preview and control buttons.
- [x] controller static asset routes for `/`, `dashboard.js` and `dashboard.css` without external frontend dependencies.
- [x] newline-delimited controller event snapshot route and persisted `events.ndjson` for downstream dashboards.
- [x] dependency-free SSE controller status stream with event IDs and heartbeat comments.
- [x] browser dashboard EventSource updates with automatic reconnect and polling fallback.
- [x] TLS 1.2+ controller socket configuration with operator-supplied certificate and key paths.
- [x] hashed bearer-token identities with viewer/operator/admin role enforcement.
- [x] bounded durable NDJSON controller audit history with rotation.
- [x] authenticated live SSH/Docker/filesystem access-probe endpoint with execution timeout.
- [x] hardened systemd service artifact with restart policy and filesystem restrictions.
- [x] TLS reverse-proxy example with SSE flushing and security headers.


## Phase 10 - RF, antennas and microwave networks [started]

Research basis is documented in `RF_CIRCUIT_RESEARCH.md`.

### RF network mathematics and interchange
- [x] complex two-port Z <-> S conversion with configurable real reference impedance.
- [x] S <-> ABCD conversion.
- [x] ABCD/two-port cascade.
- [x] ideal/lossy transmission-line ABCD and S-parameter baseline.
- [x] Hammerstad/Jensen-style quasi-static microstrip impedance/effective-permittivity baseline.
- [x] Touchstone v1 S2P RI/MA/DB reader and RI writer.
- [x] generic N-port matrix representation.
- [x] Touchstone SnP v2 import/export.
- [x] arbitrary positive-real per-port impedance renormalization.
- [x] complex-reference power-wave impedance renormalization.
- [x] de-embedding and reference-plane shifts.
- [x] mixed-mode/differential S-parameters.
- [x] N-port passivity and reciprocity checking baseline.
- [x] passivity enforcement and broadband causality checking.
- [x] two-port Rollett K and mu stability metrics.
- [x] source/load gain and stability circles.
- [x] low-pass L-match synthesis utility baseline.
- [x] adaptive broadband network sweep with midpoint interpolation-error refinement.
- [ ] rational/vector-fitting model-order reduction.
- [x] transducer-gain and load-pull sampling baseline.

### Quasi-static conductor extraction
- [x] slender-filament PEEC DC resistance extraction.
- [x] PEEC self/mutual partial-inductance matrix with numerical Neumann integration.
- [x] complex PEEC branch-impedance/current solve.
- [x] PEEC coefficient-of-potential/capacitance extraction.
- [x] round-wire PEEC skin-effect resistance baseline.
- [ ] proximity-effect conductor models.
- [ ] surface/volume PEEC discretization for PCB/package conductors.

### Antenna analysis
- [x] sinusoidal thin-wire dipole far-field radiation baseline.
- [x] numerical radiated power and radiation resistance integration.
- [x] directivity and normalized pattern extraction.
- [x] thin-wire EFIE Method-of-Moments current/input-impedance reference baseline with symmetry/passivity gates.
- [x] coupled parallel-wire MoM baseline with independent segmentation, delta-gap feeds and lumped series loads.
- [x] arbitrary-orientation disjoint straight-wire MoM with dyadic free-space Green kernel and rotational-invariance validation.
- [ ] NEC-grade thin-wire MoM accuracy/convergence across canonical antenna benchmarks.
- [x] endpoint-connected arbitrary wire geometry with explicit junction KCL constraints.
- [ ] general interior wire junctions/meshed conductor junction topology.
- [x] wire lumped loads and lossy two-port transmission-line sections.
- [x] infinite PEC ground-plane image-kernel baseline.
- [ ] finite/conductive ground models.
- [x] feed impedance/reactance and VSWR/return-loss extraction.
- [x] radiation efficiency including round-wire skin loss and equivalent distributed dielectric loss.
- [x] polarization/axial-ratio analysis.
- [x] phased-array/array-factor framework.
- [x] multi-feed impedance/S-matrix extraction for mutual impedance/coupling.
- [x] bounded derivative-free antenna optimization driver.

### Full-wave RF / antenna bridge
- [x] reuse 3-D Cartesian FDTD full-wave engine from Phase 5.
- [x] reuse TEM-like port/S-parameter extraction from Phase 5.
- [x] reuse NF2FF and SAR post-processing from Phase 5.
- [x] 2-D driven frequency-domain edge-element FEM Maxwell solver baseline.
- [ ] production 3-D driven edge-element RF FEM solver.
- [x] lowest-order first-kind Nedelec Tri3 edge-element basis with global edge orientation.
- [ ] tetrahedral Nedelec/Whitney edge-element basis for 3-D RF FEM.
- [ ] wave-port eigenmode solver.
- [ ] lumped-port frequency-domain excitation/termination.
- [ ] resonance/eigenmode/Q solver for RF cavities/antennas.
- [ ] adaptive frequency sweep and reduced-order model.
- [ ] RF full-wave adaptive mesh refinement.
- [ ] conductor/surface-impedance and skin-effect boundary model.
- [ ] frequency-dependent anisotropic tensor materials in RF FEM.

### EMC, cables and field/circuit coupling
- [x] quasi-static PEEC filament partial R/L extraction and conductor impedance solve baseline.
- [ ] graded rectilinear RF mesh generation.
- [ ] thin-wire/subcell wire model embedded in FDTD.
- [ ] multiconductor transmission-line network embedded in FDTD.
- [ ] Huygens/equivalent surface excitation.
- [ ] Hertzian dipole field source.
- [ ] thin-slot and impedance-sheet models.
- [x] shielding-effectiveness magnitude utility baseline.
- [ ] broadband EMC transfer-function workflow.
- [ ] circuit load/termination coupling to FDTD/FEM ports.
- [ ] bidirectional EM <-> SPICE co-simulation.

## Phase 11 - SPICE-class circuit and compact-model simulation [started]

Research basis is documented in `RF_CIRCUIT_RESEARCH.md`.

### Circuit kernel and analyses
- [x] clean-room modified nodal analysis correctness baseline.
- [x] resistor stamping.
- [x] capacitor stamping.
- [x] inductor/branch-current stamping.
- [x] independent voltage sources.
- [x] independent current sources.
- [x] Newton-Raphson nonlinear DC operating point.
- [x] complex small-signal AC analysis.
- [x] backward-Euler transient analysis.
- [x] sparse MNA stamping through shared CSR + ILU(0)-GMRES backend with dense fallback.
- [x] memoryless multi-terminal compact-device residual/Jacobian callback seam.
- [x] dynamic-state charge-based DAE residual/Jacobian device abstraction.
- [x] memoryless nonlinear multi-terminal residual/Jacobian compact-device seam.
- [x] PN-junction voltage limiting primitive.
- [x] source stepping.
- [x] gmin stepping/homotopy schedule.
- [x] adaptive transient timestep with local truncation error control.
- [x] trapezoidal integration.
- [x] BDF2 transient integration baseline.
- [x] adaptive variable-step BDF2/Gear integration with exact unequal-step coefficients and Richardson LTE control.
- [x] DC voltage-source sweep.
- [x] general parameter sweep.
- [x] temperature sweep.

- [x] logarithmic AC sweep framework.
- [x] multi-source noise sweep/integration framework with per-source contributions.
- [x] resistor thermal-noise output-referred baseline.
- [x] pole-zero analysis.
- [x] numerical DC sensitivity baseline.
- [x] small-signal second/third-order distortion baseline from nonlinear transfer derivatives.
- [x] Fourier/THD measurements.
- [x] periodic steady state baseline / harmonic balance remains open.
- [x] Monte-Carlo resistor tolerancing baseline.

### Devices and compact models
- [x] Shockley diode DC/small-signal/transient baseline.
- [x] MOSFET Level-1/Shichman-Hodges DC/small-signal baseline with channel-length modulation.
- [x] BJT Ebers-Moll baseline.
- [ ] Gummel-Poon BJT model.
- [x] JFET baseline.
- [ ] MOS level 2/3 baseline.
- [ ] VDMOS/power MOS model.
- [x] voltage-controlled switch baseline.
- [x] VCVS/VCCS/CCVS/CCCS controlled sources.
- [x] mutual-inductor coupling baseline.
- [x] ideal lossless transformer convenience model using power-conserving controlled-source MNA stamps.
- [ ] nonlinear magnetic core/hysteresis transformer model.
- [x] lossless/lossy sampled TEM transmission-line circuit elements for AC MNA.
- [x] programmatic nonlinear static-device callback.
- [x] SPICE B-source/equation-defined voltage/current source syntax with numerical Jacobian linearization.
- [x] coupled electrothermal/self-heating device terminal framework.
- [x] resistor thermal and diode shot-noise contributions.
- [x] BJT/MOS/JFET compact-model noise contributions baseline.

### SPICE/netlist/model interchange
- [x] SPICE engineering suffix parser (`meg`, `k`, `m`, `u`, `n`, `p`, `f`).
- [x] basic R/C/L/V/I/D/M netlist parser.
- [x] `.MODEL` parser baseline for diode, NMOS/PMOS, NPN/PNP, NJF/PJF and voltage-controlled switch models.
- [x] `.SUBCKT` hierarchy and parameter passing.
- [x] `.PARAM` and arithmetic expression evaluator.
- [x] `.FUNC` user functions.
- [x] `.INCLUDE` and `.LIB` sections with relative-path and named-section handling.
- [x] initial conditions / `.IC` / `.NODESET`.
- [x] dependent-source `POLY(1)` / `TABLE` syntax baseline.
- [ ] PSpice/LTspice/HSPICE compatibility modes.
- [x] ASCII SPICE RAW transient/AC writer/reader.
- [ ] optional external ngspice adapter for cross-validation.
- [ ] optional external Xyce adapter for cross-validation/scaling.
- [x] OSDI dynamic-library/version/descriptor discovery seam.
- [ ] OSDI descriptor evaluation/stamping adapter.
- [ ] OpenVAF Verilog-A compilation workflow.
- [ ] validation with BSIM-CMG/PSP/HICUM-class public Verilog-A models.

### RF circuit integration
- [x] two-port S/Z/ABCD mathematics shared with Phase 10.
- [x] Touchstone S2P parser shared with Phase 10.
- [x] generic sampled N-port Touchstone device in AC MNA.
- [x] microstrip/stripline/coplanar sampled circuit elements.
- [x] rectangular-waveguide TE10 sampled circuit element baseline.
- [ ] equivalent discontinuity/modal circuit elements.
- [x] S-parameter N-port to admittance stamping for AC circuit analysis.
- [x] direct circuit port S-parameter extraction.
- [x] shared two-port K/mu stability metrics.
- [x] source/load stability circles.
- [x] constant-noise circle calculation from two-port noise parameters.
- [x] shared L-match synthesis baseline.
- [x] sampled field/network N-port -> AC circuit component bridge.
- [ ] circuit termination -> full-wave field-solver port bridge.

## Optimization gates applied continuously

A feature is optimized only after a correctness baseline exists.

- [x] SoA LBM storage.
- [x] in-place LBM streaming.
- [x] compact distributed LBM halo payload.
- [x] OpenMP cell-owned FVM accumulation avoiding atomics.
- [x] bounded reconstruction implemented without global locks.
- [x] shared CSR/Krylov layer to replace solver-specific dense/Jacobi paths.
- [x] eliminate transient heap allocations from hot FVM loops with reusable workspaces.
- [x] matrix assembly reuse when sparsity pattern is unchanged.
- [x] Morton-order cell permutation utility baseline.
- [ ] connectivity-preserving face/cell storage reorder integrated into solver meshes.
- [x] SIMD fixed-width vector math.
- [x] mixed precision with iterative refinement.
- [ ] SYCL FVM/FEM sparse algebra. (FVM resident operators, pressure CG, nonsymmetric BiCGStab momentum/scalar transport, fixed-pressure outlets and first-order SIMPLE/PISO/PIMPLE are implemented through v0.15.10; FEM integration and physical-GPU qualification remain open.)
- [ ] kernel fusion guided by profiler measurements.
- [ ] communication/computation overlap beyond LBM.
- [ ] automatic backend/autotuning profiles saved per device.
- [ ] regression performance thresholds for representative CPU/GPU hardware.

## Phase 12 - CST-class electromagnetic, particle and bioelectromagnetic coverage [started]

This phase maps the additional solver domains identified in the CST-style coverage audit in `CST_COVERAGE_AUDIT.md`.
Existing capabilities remain owned by their original phases; this phase tracks only the missing cross-domain expansion work and newly added reference baselines.

### Frequency-domain and eigenmode electromagnetics
- [x] driven 1-D complex frequency-domain Maxwell/Helmholtz PEC reference solver with conductive loss.
- [x] generalized 1-D PEC cavity electromagnetic eigenmode solver with spatial dielectric loading.
- [x] 2-D Tri3 curl-conforming Nedelec edge-element Maxwell assembly.
- [x] 3-D tetrahedral curl-conforming Nedelec/Whitney Maxwell assembly.
- [x] driven 3-D frequency-domain sparse complex Maxwell solve.
- [x] rectangular PEC wave-port eigenmode extraction and 1 W modal normalization baseline.
- [x] cavity/resonator Q extraction including dielectric and conductor loss.
- [x] adaptive RF frequency sweep + stable common-pole rational reduced-order model.
- [x] RF FEM adaptive mesh refinement driven by face-jump indicators, Dorfler marking and conforming longest-edge star bisection.

### Integral-equation and asymptotic electromagnetics
- [ ] RWG surface-current MoM for PEC surfaces.
- [ ] dielectric surface integral equations.
- [ ] MLFMM acceleration for electrically large integral-equation systems.
- [ ] physical-optics current approximation for large smooth conductors.
- [ ] shooting-and-bouncing-rays high-frequency propagation/RCS solver.
- [ ] hybrid MoM/FEM/FDTD/asymptotic domain coupling.

### Low-frequency machines and conductors
- [ ] transient magneto-quasistatic A-phi formulation.
- [x] nonlinear B-H curve material integration in magnetostatic FEM.
- [ ] nonlinear B-H curve material integration in harmonic eddy-current FEM.
- [x] stranded-coil excitation with flux-linkage/inductance export for circuit coupling.
- [ ] solid-conductor coil excitation with external circuit coupling.
- [ ] moving-band/sliding-interface electrical-machine formulation.
- [x] force/torque extraction from Maxwell-stress boundary integration.

### Signal/power integrity, cables and EMC
- [x] multiconductor transmission-line cable/harness solver with frequency-dependent RLCG and full matrix loads.
- [x] cable shield/transfer-impedance model and effective-height field-to-cable coupling.
- [x] eye-diagram/BER-oriented signal-integrity post-processing.
- [x] PDN impedance/IR-drop/decoupling optimization workflow.
- [x] ESD/BCI/lightning waveform source library and standardized EMC probes.
- [ ] installed-antenna/co-site hybrid coupling workflow.

### Charged particles, PIC, plasma and wakefields
- [x] non-relativistic 3-D Lorentz-force charged-particle tracker using the Boris pusher.
- [x] periodic 1-D electrostatic PIC baseline with CIC deposition and spectral Poisson field solve.
- [x] periodic 2-D electrostatic PIC baseline with bilinear CIC deposition, spectral Poisson field solve and field gather.
- [x] periodic 3-D electrostatic PIC baseline with trilinear CIC deposition, spectral Poisson field solve and field gather.
- [x] 2-D spectral charge-conserving current reconstruction from old/new particle charge density.
- [x] 3-D spectral charge-conserving current reconstruction from old/new particle charge density.
- [x] local finite-volume 3-D current reconstruction satisfying periodic discrete continuity.
- [x] 3-D particle geometry/deposition validation foundation for future full Yee EM-PIC.
- [x] relativistic Boris particle pusher using proper velocity.
- [x] Vay pusher for ultra-relativistic crossed-field regimes.
- [x] periodic 1-D/3V self-consistent electromagnetic Yee-grid PIC baseline.
- [x] periodic 2-D/3-V electromagnetic PIC baseline with centered Maxwell curl update and Vay particle push.
- [x] periodic 3-D/3-V electromagnetic PIC baseline with centered Maxwell curl update and Vay particle push.
- [x] true staggered 2-D Yee electromagnetic PIC baseline with CFL guard and component-specific field locations.
- [x] integrated electric-wall and absorbing-sponge field boundaries in the staggered EM-PIC update.
- [x] integrated absorbing/specular particle walls with secondary-yield reporting in the staggered EM-PIC step.
- [x] transverse 2-D/3-V current deposition (`Jz`) coupled into the electromagnetic field update.
- [x] full 3-D spectral current coupling (`Jx/Jy/Jz`) into the electromagnetic field update.
- [x] selectable local finite-volume current coupling in the staggered 3-D EM-PIC update.
- [x] stable 3-D PIC particle sorting by cell with cell offsets and original-index recovery.
- [x] 3-D particle guard-halo classification for face/edge/corner domain exchange.
- [x] optional sorted particle storage in the staggered 3-D EM-PIC update with diagnostics.
- [x] Cartesian 3-D PIC domain decomposition descriptors with uneven block coverage.
- [x] serial particle migration bucketing between PIC subdomains.
- [x] serial scalar guard-cell exchange producing ghost-padded decomposed blocks.
- [x] logical-rank topology mapping for PIC subdomains before MPI transport.
- [x] rank-addressed particle/guard transport envelopes with diagnostics.
- [x] deterministic in-memory PIC transport reproducing migration and guard exchange through send/receive-style inboxes.
- [x] serialized rank-addressed PIC transport envelope contract with deterministic roundtrip and payload diagnostics.
- [x] destination-rank serialized exchange planning for particle/guard messages.
- [x] six-component electromagnetic field-vector guard exchange for Ex/Ey/Ez/Bx/By/Bz blocks.
- [x] serialized rank-local distributed PIC exchange round combining particle migration and scalar field guards.
- [x] distributed PIC exchange-round and EM field-guard CLI/regression coverage.
- [x] rank-local distributed staggered 3-D PIC timestep bridge combining local field push, serialized migration and EM guard exchange.
- [x] Geant4-inspired track/step/process transport baseline with geometry-limited steps and process-limited GPIL selection.
- [x] continuous energy-loss and discrete secondary-production process hooks with per-step scoring.
- [x] slab-region material lookup, boundary crossing and production-cut handling for particle-through-matter smoke tests.
- [x] stochastic sampled physical interaction lengths for discrete transport processes.
- [x] physics-list bundle with clean production-cut override.
- [x] BVH-accelerated region lookup for axis-aligned transport regions.
- [x] sensitive-detector hit collection and dose-grid/SAR/Pennes projection coupling.
- [ ] real MPI PIC envelope send/receive validation across ranks.
- [x] trilinear 3-D E/B gather with Monte-Carlo collision coupling inside the 3-D EM-PIC step.
- [x] true staggered 3-D Yee electromagnetic PIC baseline with CFL guard and component-specific Ex/Ey/Ez/Bx/By/Bz locations.
- [x] component-offset trilinear 3-D E/B gather for staggered Yee particles.
- [x] reusable 3-D electric-wall and absorbing-sponge EM-PIC field-boundary primitives.
- [x] configurable polynomial-order 3-D absorbing-sponge profile.
- [x] integrated 3-D absorbing/specular particle walls with secondary-yield reporting in the staggered EM-PIC step.
- [x] particle boundary interaction with absorption, specular reflection and secondary-emission yield baseline.
- [x] reusable 2-D electric-wall and absorbing-sponge PIC field-boundary primitives.
- [x] Monte-Carlo collision model for neutral gas/plasma interactions.
- [x] MCC ionization source coupling into the PIC particle population.
- [x] multi-species plasma chemistry reaction-network coupling baseline with charge-conservation diagnostics.
- [x] causal beam wake-potential convolution and wake-impedance solver.
- [x] Paschen gas-breakdown and parallel-plate multipactor threshold workflow baseline.

### Bioelectromagnetics
- [x] RMS electric-field -> local SAR material conversion helper.
- [x] implicit 2-D Pennes bioheat solver with perfusion, metabolic heat and spatial SAR source.
- [x] SAR -> bioheat coupling baseline with periodic/fixed-temperature thermal boundaries.
- [x] heterogeneous anatomical voxel material ingestion.
- [x] 1 g / 10 g SAR map -> Pennes mesh transfer on heterogeneous tissues.
- [ ] temperature-dependent dielectric/perfusion feedback to the EM solve.
- [ ] implant/wearable exposure validation cases.

### Photonics and specialized EM
- [ ] frequency-domain dispersive photonic waveguide mode solver.
- [ ] periodic/Bloch boundary conditions and photonic band-structure solver.
- [ ] metasurface/generalized sheet transition-condition model.
- [ ] nonlinear optical material polarization models.
- [ ] dedicated optical full-wave validation beyond ray/POP/FDTD baselines.

## Phase 13 - multibody, rigid-body and granular dynamics [started]

Clean-room capability/reference family: Project Chrono. This phase adds constrained 6-DOF mechanics and granular/contact dynamics that are distinct from the existing FEM structural and transported-particle paths.

### Rigid-body state and integration
- [x] 6-DOF rigid-body state with quaternion orientation, force/torque accumulators and world-space inertia application.
- [x] rigid-body system orchestration with gravity, external forces, attached collision shapes and deterministic stepping.
- [x] semi-implicit Euler rigid-body integrator.
- [x] velocity-Verlet rigid-body integrator baseline.
- [ ] implicit Newmark multibody integrator.
- [ ] generalized-alpha/HHT multibody integrator.

### Constraint and joint framework
- [x] reusable Jacobian-row constraint representation with projected Gauss-Seidel impulse solve.
- [x] distance constraint.
- [x] spherical joint.
- [x] revolute joint with one free angular DOF.
- [x] prismatic joint with one free translational DOF.
- [x] fixed joint.
- [x] gear-ratio angular constraint.
- [x] motor/actuator constraint family.
- [ ] articulated reduced-coordinate solver.

### Collision and contact
- [x] sphere AABB generation and sweep-and-prune broad phase.
- [x] BVH broad-phase collision candidate generation.
- [x] sphere-sphere narrow-phase contact geometry.
- [x] sphere-plane narrow-phase contact geometry.
- [x] smooth penalty contact with damping.
- [x] non-smooth unilateral impulse/contact complementarity baseline.
- [x] Coulomb friction impulse/force limiting.
- [x] convex GJK/EPA narrow phase.
- [x] mesh/triangle collision and persistent contact manifolds.

### Granular / DEM
- [x] explicit spherical DEM system.
- [x] Hertz-type nonlinear normal contact with tangential damping/friction.
- [x] rolling-resistance torque.
- [x] cohesive normal contact force baseline.
- [x] history-dependent Mindlin tangential spring.
- [x] bonded particles with progressive tensile/shear damage and fracture/bond failure.
- [x] GPU/SYCL particle neighbor search and contact kernels.
- [ ] distributed-memory DEM domain decomposition (slab ownership, migration/ghost planning and MPI exchange implemented; end-to-end multi-rank contact/integration runtime validation pending).

### Multiphysics coupling
- [x] spherical Stokes-drag/reaction-force primitive for conservative CFD/DEM coupling.
- [ ] resolved CFD <-> rigid-body surface traction/force/torque coupling.
- [ ] unresolved CFD <-> many-particle drag/void-fraction coupling.
- [ ] FEM flexible-body <-> multibody coupling.
- [x] particle <-> sphere/plane contact-geometry bridge through the explicit DEM system.

