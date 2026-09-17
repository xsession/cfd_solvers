# Unified Solver Integration Tracker

Last research refresh: 2026-09-16

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
- [ ] explicit NUMA node/thread affinity and placement policy.
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
- [x] block CSR / vector-valued matrices.
- [x] ILU(0) preconditioner.
- [x] algebraic multigrid preconditioner adapter with externally supplied hierarchy/cycle callback.
- [x] geometric multigrid for structured grids.
- [ ] SYCL SpMV and Krylov kernels.
- [x] distributed local-row CSR with explicit halo-column values and halo-aware SpMV baseline.
- [ ] MPI-owned distributed sparse-vector exchange integrated directly into Krylov iterations.
- [x] mixed-precision iterative refinement.

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
- [ ] cumulant collision option.
- [ ] Smagorinsky/subgrid LES option.
- [x] thermal DDF lattice.
- [x] passive scalar DDF lattice.
- [ ] multiphase/free-surface extension.

### Boundaries and geometry
- [x] halfway bounce-back.
- [x] moving wall.
- [x] velocity inlet.
- [x] pressure/density outlet.
- [x] channel and cavity regressions.
- [ ] interpolated curved-wall bounce-back.
- [ ] immersed-boundary particles.
- [ ] two-way particle coupling.
- [x] STL/triangle-mesh voxelization.
- [ ] GPU voxelization.
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
- [ ] compressed population storage with explicit conservation/error gates.
- [ ] FP16/BF16/custom packed storage experiments.
- [ ] kernel-fusion/autotuning database by device.
- [ ] multiple GPUs per rank.
- [ ] adaptive domain repartitioning.

### Visualization/output concepts worth carrying forward
- [x] legacy VTK structured-grid export for 2-D LBM scalar/vector fields.
- [x] LBM vorticity derived-field baseline.
- [ ] LBM Q-criterion derived field.
- [x] 2-D LBM streamline extraction baseline.
- [x] LBM horizontal slice extraction baseline.
- [x] headless PGM scalar-field image output baseline.
- [ ] optional lightweight interactive viewer kept separate from solver core.

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
- [ ] second-order temporal schemes integrated into production FVM equation solvers.

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
- [ ] pressure AMG/multigrid path.
- [x] adaptive timestep from Courant number.
- [ ] production pressure-velocity solver integration of second-order backward/CN time stepping.

### Transport and thermo
- [x] scalar advection-diffusion equation.
- [x] energy/enthalpy equation.
- [x] temperature-dependent material properties.
- [x] multi-species mass transport.
- [x] buoyancy/Boussinesq.
- [x] compressible equation-of-state framework.
- [x] compressible pressure-energy coupling.
- [x] MUSCL-minmod + Rusanov shock-capturing compressible baseline.
- [ ] characteristic high-order/WENO compressible schemes.

### Turbulence
- [x] laminar/turbulence runtime interface.
- [x] mixing-length reference.
- [x] Spalart-Allmaras eddy-viscosity constitutive helper.
- [ ] Spalart-Allmaras transport equation.
- [x] k-epsilon eddy-viscosity constitutive helper.
- [ ] k-epsilon transport equations.
- [x] k-omega/SST eddy-viscosity constitutive helper.
- [ ] k-omega/SST transport equations.
- [x] Reynolds-stress Boussinesq reconstruction baseline.
- [ ] Reynolds-stress transport model framework.
- [x] LES filters and Smagorinsky/WALE.
- [x] wall functions and y+ diagnostics.
- [x] DES length-scale switching baseline.
- [ ] full DES/hybrid RANS-LES transport framework.

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
- [ ] solid heat conduction.
- [x] conjugate fluid-solid heat transfer.
- [x] gray two-surface exchange baseline.
- [x] general diffuse-gray surface-to-surface radiosity/view-factor network baseline.
- [x] optically-thin participating-media source baseline.
- [x] participating-media radiation source-model interface with optically-thin implementation.
- [ ] radiation/chemistry energy coupling.

### Mesh motion/adaptation
- [x] mesh-motion field.
- [x] ALE flux correction.
- [x] rigid-body motion.
- [ ] topology change interface.
- [ ] local refinement/coarsening.
- [ ] error-indicator-driven AMR.
- [x] conservative field remap after topology changes.

### Case/IO/post-processing
- [x] unit-aware native YAML/JSON case-parameter parsing baseline.
- [x] OpenFOAM mesh importer.
- [x] OpenFOAM field importer/exporter where licensing/interoperability permits.
- [x] VTK/VTU output.
- [ ] HDF5 checkpoint/field output.
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
- [ ] SYCL FDTD kernels.
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

## Phase 9 - interoperability, workflow and UX [planned]

- [x] common case schema with units and validation.
- [x] material database schema shared by CFD/FEM/FDTD/optics/electrochemistry.
- [x] mesh import: Gmsh.
- [x] mesh import: VTK/VTU.
- [x] mesh import: OpenFOAM polyMesh.
- [x] geometry import: STL/OBJ.
- [x] output: VTK/VTU.
- [ ] output: HDF5/XDMF.
- [ ] Python bindings.
- [x] parameter sweep runner.
- [x] optimization/inverse-problem runner.
- [x] restartable workflow graph.
- [x] provenance metadata in every result file.
- [x] benchmark/validation catalog CLI.


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
- [ ] SYCL FVM/FEM sparse algebra.
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
- [ ] 3-D tetrahedral curl-conforming Nedelec/Whitney Maxwell assembly.
- [ ] driven 3-D frequency-domain sparse complex Maxwell solve.
- [ ] wave-port eigenmode extraction and modal normalization.
- [ ] cavity/resonator Q extraction including dielectric and conductor loss.
- [ ] adaptive RF frequency sweep + rational reduced-order model.
- [ ] RF FEM adaptive mesh refinement driven by field/error indicators.

### Integral-equation and asymptotic electromagnetics
- [ ] RWG surface-current MoM for PEC surfaces.
- [ ] dielectric surface integral equations.
- [ ] MLFMM acceleration for electrically large integral-equation systems.
- [ ] physical-optics current approximation for large smooth conductors.
- [ ] shooting-and-bouncing-rays high-frequency propagation/RCS solver.
- [ ] hybrid MoM/FEM/FDTD/asymptotic domain coupling.

### Low-frequency machines and conductors
- [ ] transient magneto-quasistatic A-phi formulation.
- [ ] nonlinear B-H curve material integration in magnetostatic/eddy-current FEM.
- [ ] stranded/solid coil excitation and circuit coupling.
- [ ] moving-band/sliding-interface electrical-machine formulation.
- [ ] force/torque extraction from Maxwell stress/virtual work.

### Signal/power integrity, cables and EMC
- [ ] multiconductor transmission-line cable/harness solver with frequency-dependent RLCG.
- [ ] cable shield/transfer-impedance model and field-to-cable coupling.
- [ ] eye-diagram/BER-oriented signal-integrity post-processing.
- [ ] PDN impedance/IR-drop/decoupling optimization workflow.
- [ ] ESD/BCI/lightning waveform source library and standardized EMC probes.
- [ ] installed-antenna/co-site hybrid coupling workflow.

### Charged particles, PIC, plasma and wakefields
- [x] non-relativistic 3-D Lorentz-force charged-particle tracker using the Boris pusher.
- [x] periodic 1-D electrostatic PIC baseline with CIC deposition and spectral Poisson field solve.
- [ ] relativistic Boris/Vay particle pusher.
- [ ] self-consistent electromagnetic Yee-grid PIC.
- [ ] particle boundary interaction/absorption/reflection/secondary-emission models.
- [ ] Monte-Carlo collision model for neutral gas/plasma interactions.
- [ ] plasma chemistry/ionization source coupling.
- [ ] beam wake-potential and wake-impedance solver.
- [ ] multipactor and gas-breakdown threshold workflow.

### Bioelectromagnetics
- [x] RMS electric-field -> local SAR material conversion helper.
- [x] implicit 2-D Pennes bioheat solver with perfusion, metabolic heat and spatial SAR source.
- [x] SAR -> bioheat coupling baseline with periodic/fixed-temperature thermal boundaries.
- [ ] heterogeneous anatomical voxel material ingestion.
- [ ] 1 g / 10 g SAR map -> Pennes mesh transfer on heterogeneous tissues.
- [ ] temperature-dependent dielectric/perfusion feedback to the EM solve.
- [ ] implant/wearable exposure validation cases.

### Photonics and specialized EM
- [ ] frequency-domain dispersive photonic waveguide mode solver.
- [ ] periodic/Bloch boundary conditions and photonic band-structure solver.
- [ ] metasurface/generalized sheet transition-condition model.
- [ ] nonlinear optical material polarization models.
- [ ] dedicated optical full-wave validation beyond ray/POP/FDTD baselines.

