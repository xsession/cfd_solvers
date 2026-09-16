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
- [ ] asynchronous task graph shared by FVM/FEM/FDTD/LBM.
- [ ] NUMA-aware CPU memory placement.
- [ ] SIMD abstraction for small fixed-size kernels.
- [ ] device-memory pool / scratch allocator.
- [ ] runtime backend selection without recompiling the application.
- [ ] profiling ranges/counters for every solver family.

### Shared linear algebra
- [x] matrix-free conjugate gradient.
- [x] CSR sparse matrix storage and assembly builder.
- [x] Jacobi preconditioner.
- [x] preconditioned conjugate gradient.
- [x] BiCGStab for nonsymmetric systems.
- [x] restarted GMRES.
- [ ] block CSR / vector-valued matrices.
- [x] ILU(0) preconditioner.
- [ ] algebraic multigrid adapter.
- [ ] geometric multigrid for structured grids.
- [ ] SYCL SpMV and Krylov kernels.
- [ ] distributed CSR and halo-aware SpMV.
- [ ] mixed-precision iterative refinement.

## Phase 2 - FluidX3D-class LBM [advanced baseline]

### Core numerics
- [x] D2Q9 BGK.
- [x] D3Q19 BGK.
- [x] D3Q27 BGK.
- [x] two-grid reference streaming.
- [x] single-grid in-place/pull streaming.
- [x] Guo forcing.
- [x] FP32/FP64 validation.
- [ ] MRT collision.
- [x] TRT collision on the periodic CPU two-grid reference path, including Guo forcing.
- [ ] regularized/cumulant collision option.
- [ ] Smagorinsky/subgrid LES option.
- [ ] thermal DDF lattice.
- [ ] passive scalar DDF lattice.
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
- [ ] STL/triangle-mesh voxelization.
- [ ] GPU voxelization.
- [ ] moving/rotating geometry re-voxelization.
- [ ] force/torque integration on solids.
- [ ] porous-media drag models.

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
- [ ] VTK structured-grid export for LBM fields.
- [ ] Q-criterion/vorticity derived fields.
- [ ] streamline extraction.
- [ ] slice extraction.
- [ ] headless image output.
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
- [x] minmod/van-Leer scalar and componentwise vector reconstruction baseline.
- [x] least-squares gradient with sheared-mesh manufactured validation.
- [x] corrected/snGrad operator with linear-field exactness validation.
- [x] deferred non-orthogonal Laplacian baseline.
- [ ] generic dimensioned field classes.
- [ ] generic finite-volume matrix assembly.
- [ ] run-time selectable discretization schemes.
- [ ] field old-time history and multi-step temporal schemes.

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
- [ ] adaptive timestep from Courant number.
- [ ] second-order backward/CN time integration.

### Transport and thermo
- [x] scalar advection-diffusion equation.
- [ ] energy/enthalpy equation.
- [ ] temperature-dependent material properties.
- [ ] multi-species mass transport.
- [ ] buoyancy/Boussinesq.
- [ ] compressible equation-of-state framework.
- [ ] compressible pressure-energy coupling.
- [ ] shock-capturing/high-resolution compressible schemes.

### Turbulence
- [ ] laminar/turbulence runtime interface.
- [ ] mixing-length reference.
- [ ] Spalart-Allmaras.
- [ ] k-epsilon family.
- [ ] k-omega/SST.
- [ ] Reynolds-stress model framework.
- [ ] LES filters and Smagorinsky/WALE.
- [ ] wall functions and y+ diagnostics.
- [ ] DES/hybrid RANS-LES framework.

### Multiphase/interface physics
- [ ] VOF volume fraction transport.
- [ ] bounded interface compression.
- [ ] surface tension/CSF.
- [ ] contact angle.
- [ ] phase change.
- [ ] compressible VOF.
- [ ] Euler-Euler multiphase framework.
- [ ] drift-flux model.
- [ ] thin-film model.

### Lagrangian/dispersed phase
- [ ] particle cloud container.
- [ ] particle tracking through polyhedral cells.
- [ ] drag/lift/virtual-mass force models.
- [ ] one-way/two-way coupling.
- [ ] collisions/coalescence/breakup framework.
- [ ] sprays/droplets/evaporation.
- [ ] dense-particle coupling.

### Reacting CFD
- [ ] chemistry mechanism abstraction shared with Phase 8.
- [ ] finite-rate species source coupling.
- [ ] combustion thermo state.
- [ ] laminar flame baseline.
- [ ] premixed/non-premixed model interfaces.
- [ ] soot/radiation coupling hooks.

### Heat/radiation/CHT
- [ ] solid heat conduction.
- [ ] conjugate fluid-solid heat transfer.
- [ ] surface-to-surface radiation.
- [ ] participating-media radiation model interface.
- [ ] radiation/chemistry energy coupling.

### Mesh motion/adaptation
- [ ] mesh-motion field.
- [ ] ALE flux correction.
- [ ] rigid-body motion.
- [ ] topology change interface.
- [ ] local refinement/coarsening.
- [ ] error-indicator-driven AMR.
- [ ] conservative field remap after topology changes.

### Case/IO/post-processing
- [ ] unified YAML/JSON case format.
- [ ] OpenFOAM mesh importer.
- [ ] OpenFOAM field importer/exporter where licensing/interoperability permits.
- [ ] VTK/VTU output.
- [ ] HDF5 checkpoint/field output.
- [ ] probes/sampling.
- [ ] forces/coefficients.
- [ ] residual/function-object framework.
- [ ] derived fields: vorticity, Q, Lambda2, enstrophy, y+.

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
- [ ] nonlinear geometry.
- [ ] hyperelasticity.
- [ ] plasticity/material-law interface.
- [ ] contact.
- [x] generalized sparse eigenfrequency baseline with consistent mass and fixed-free bar validation.
- [ ] structural dynamics.

### Thermal/fluid/porous
- [x] 2-D transient heat equation baseline (steady Poisson/diffusion covered separately).
- [ ] phase change/latent heat.
- [ ] porous heat/mass transport.
- [ ] FEM incompressible flow reference.
- [x] saturated Darcy porous-flow baseline.
- [ ] Brinkman porous flow.

### Electromagnetics
- [x] electrostatics.
- [x] DC conduction.
- [x] 2-D magnetostatics through out-of-plane vector potential A_z and reconstructed B field.
- [ ] eddy-current harmonic formulation.
- [ ] frequency-domain Maxwell with appropriate vector elements.
- [x] one-way DC conduction/Joule heat/transient thermal coupling on a shared Tri3 mesh.

### Multiphysics/HPC
- [ ] monolithic block systems.
- [x] partitioned coupling baseline: one-way DC/heat/elasticity and a shared fixed-point driver.
- [ ] MPI mesh partitioning.
- [ ] distributed sparse assembly.
- [ ] GPU element kernels.
- [ ] restart/checkpoint and adaptive repartition.

## Phase 5 - openEMS-class FDTD [in progress]

- [x] 1-D Yee proof solver.
- [x] 3-D Cartesian Yee grid baseline.
- [ ] cylindrical coordinate formulation.
- [ ] cylindrical multi-grid refinement.
- [x] material coefficient preprocessing.
- [x] PEC enclosure boundary baseline.
- [x] PMC boundary baseline on the 3-D Cartesian engine.
- [x] Mur absorbing boundary.
- [x] 1-D polynomial matched electric/magnetic absorbing layer with pulse-energy validation.
- [ ] UPML/CPML absorbing boundary.
- [ ] TFSF excitation.
- [x] soft electric-field source baseline.
- [x] hard field source.
- [ ] lumped R/L/C elements.
- [x] Drude ADE dispersion baseline.
- [x] Lorentz ADE dispersion baseline.
- [x] Debye ADE dispersion baseline.
- [ ] anisotropic materials.
- [x] 1-D TEM-like wave-port decomposition and complex S-parameter baseline.
- [x] probes/time signals/DFT monitors.
- [ ] NF2FF transform.
- [ ] SAR 1g/10g calculation.
- [ ] HDF5 output.
- [x] legacy ASCII VTK E/H field output baseline.
- [ ] SYCL FDTD kernels.
- [ ] MPI + SYCL domain decomposition.
- [ ] geometry/material bridge shared with optics/FEM.

## Phase 6 - Optiland-class optics [in progress]

### System/modeling
- [x] vector reflection/refraction proof functions.
- [x] sequential optical system/surface list baseline.
- [ ] coordinate transforms/decenters/tilts.
- [x] spherical surfaces.
- [x] rotational conics/even aspheres with sag, normal and sequential intersection baseline.
- [ ] polynomial/freeform surfaces.
- [x] circular surface apertures/vignetting baseline.
- [ ] explicit aperture stops and general aperture shapes.
- [ ] wavelength and field models.
- [x] glass/material dispersion database.

### Ray tracing and physics
- [x] paraxial trace.
- [x] sequential real-ray trace with CPU batch execution.
- [x] polarization/Jones/Stokes representation.
- [ ] birefringence.
- [x] Fresnel coatings.
- [x] multilayer coating transfer matrices.
- [ ] non-sequential ray tracing.
- [ ] ghost/multi-sequence tracing.
- [ ] absorption/scatter.
- [ ] physical-optics propagation / diffraction.
- [x] Gaussian beam ABCD propagation, refraction-index convention and thin-lens focusing.

### Analysis
- [x] first-order/paraxial properties.
- [x] RMS spot-radius analysis baseline.
- [ ] spot-diagram sampling/reporting.
- [ ] ray fans.
- [ ] wavefront error.
- [ ] Zernike decomposition.
- [ ] PSF.
- [ ] MTF.
- [ ] distortion/chromatic analysis.
- [ ] stray-light/ghost analysis.

### Optimization/tolerance/interoperability
- [ ] variable/operand merit function.
- [ ] local optimizers.
- [ ] global optimizers.
- [ ] automatic differentiation path.
- [ ] Monte Carlo tolerancing.
- [ ] sensitivity analysis.
- [ ] material/glass search.
- [ ] Zemax ZMX import/export.
- [ ] CODE V SEQ import/export.
- [ ] OSLO LEN import/export.
- [ ] JSON native format.
- [ ] SYCL batched ray-surface kernels.

## Phase 7 - chemistry, electrochemistry and corrosion [started]

Research basis is documented in `CHEMISTRY_CORROSION_RESEARCH.md`.

### Chemistry/thermodynamics inspired by Cantera/Reaktoro/PHREEQC
- [x] species metadata: name, molar mass, charge.
- [x] elementary mass-action reaction network.
- [x] Arrhenius forward rates.
- [x] reversible mass-action reactions with concentration-form van't Hoff equilibrium constants.
- [ ] third-body/falloff reactions.
- [ ] general reaction orders.
- [ ] thermodynamic polynomial/property model interface.
- [ ] ideal gas/ideal solution phases.
- [x] ideal dilute acid/base speciation with family mass balance, electroneutrality and pH.
- [ ] non-ideal activity model interface.
- [ ] aqueous activity models.
- [ ] Pitzer model.
- [ ] SIT model.
- [ ] chemical-potential/Gibbs minimization equilibrium.
- [x] small isothermal stiff reaction networks: adaptive implicit Euler, damped Newton and positivity checks.
- [ ] multiphase equilibrium/precipitation/dissolution.
- [ ] surface complexation.
- [ ] ion exchange.
- [ ] thermochemical database import layer.
- [ ] PHREEQC database reader or optional PhreeqcRM adapter.
- [ ] Cantera YAML mechanism importer or optional Cantera adapter.
- [ ] Reaktoro dynamic adapter kept separate from clean core.
- [ ] reaction-path analysis.
- [ ] reactor-network abstractions.

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
- [ ] porous-electrode effective transport.
- [ ] electronic solid-phase conduction.
- [ ] immersed electrode/electrolyte interfaces.
- [x] Butler-Volmer interfacial species-flux boundary on PolyMesh.
- [ ] fully implicit nonlinear electrode/transport coupling.
- [ ] Marcus/Marcus-Hush-Chidsey kinetics.
- [ ] double-layer capacitance.
- [ ] generalized modified PNP / finite-size effects.
- [ ] GPU kernels for electrochemical transport.
- [ ] MPI decomposition and AMR.

### Corrosion-specific physics
- [x] galvanic/mixed-potential solver with multiple reactions/metals.
- [x] anodic metal dissolution species source via Faradaic Butler-Volmer boundary flux.
- [ ] cathodic oxygen/hydrogen reaction models.
- [ ] concentration-dependent Nernst equilibrium potentials.
- [ ] passive-film/oxide growth state variable.
- [ ] passivation/transpassive kinetics.
- [ ] precipitation/product-layer coupling.
- [ ] pH and aqueous speciation coupling.
- [ ] localized/pitting initiation model.
- [ ] moving metal/electrolyte interface.
- [ ] phase-field corrosion interface.
- [ ] stress-assisted corrosion coupling.
- [ ] Faradaic geometry recession/remeshing.
- [ ] cathodic-protection/anode models.
- [ ] electrochemical impedance spectroscopy small-signal solver.
- [ ] corrosion benchmarks against published PNP/FEM cases.

## Phase 8 - coupled multiphysics [started]

- [x] explicit field registry with units/location/topology metadata.
- [x] conservative cell-to-cell transfer baseline using exact-overlap 1-D remap.
- [x] conservative face-to-cell transfer on `PolyMesh` integrated fluxes.
- [ ] FEM<->FVM projection/interpolation.
- [ ] structured-grid<->unstructured transfer.
- [x] fixed-point partitioned coupler.
- [x] Aitken relaxation.
- [ ] monolithic block coupling interface.
- [x] DC conduction -> Joule heating -> transient FEM thermal coupling on a shared Tri3 mesh.
- [x] thermal -> small-strain plane-stress/plane-strain deformation on a shared Tri3 mesh.
- [ ] deformation -> optical surfaces.
- [ ] CFD -> thermal -> structural -> optics chain.
- [ ] electrochemistry -> heat generation -> thermal.
- [ ] corrosion recession -> mesh motion -> CFD/structural update.
- [ ] battery/electrolyzer porous electrochemistry + thermal + flow.

## Phase 9 - interoperability, workflow and UX [planned]

- [ ] common case schema with units and validation.
- [ ] material database schema shared by CFD/FEM/FDTD/optics/electrochemistry.
- [ ] mesh import: Gmsh.
- [ ] mesh import: VTK/VTU.
- [ ] mesh import: OpenFOAM polyMesh.
- [ ] geometry import: STL/OBJ.
- [ ] output: VTK/VTU.
- [ ] output: HDF5/XDMF.
- [ ] Python bindings.
- [ ] parameter sweep runner.
- [ ] optimization/inverse-problem runner.
- [ ] restartable workflow graph.
- [ ] provenance metadata in every result file.
- [ ] benchmark/validation catalog CLI.

## Optimization gates applied continuously

A feature is optimized only after a correctness baseline exists.

- [x] SoA LBM storage.
- [x] in-place LBM streaming.
- [x] compact distributed LBM halo payload.
- [x] OpenMP cell-owned FVM accumulation avoiding atomics.
- [x] bounded reconstruction implemented without global locks.
- [x] shared CSR/Krylov layer to replace solver-specific dense/Jacobi paths.
- [ ] eliminate transient heap allocations from hot FVM loops with reusable workspaces.
- [x] reusable scalar FVM gradient/reconstruction/divergence buffers with reference-parity profiling.
- [x] scalar-transport matrix and ILU reuse while coefficients remain unchanged, with boundary/flux invalidation.
- [ ] cache-friendly face/cell reordering.
- [ ] SIMD fixed-width vector math.
- [ ] mixed precision with iterative refinement.
- [ ] SYCL FVM/FEM sparse algebra.
- [ ] kernel fusion guided by profiler measurements.
- [ ] communication/computation overlap beyond LBM.
- [ ] automatic backend/autotuning profiles saved per device.
- [ ] regression performance thresholds for representative CPU/GPU hardware.
