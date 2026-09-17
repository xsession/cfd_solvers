# Changelog

## 0.8.0 - RF/antenna and SPICE-class solver expansion

- Integrated the latest remote CPU continuation with the locally validated advanced FDTD/multiphysics work without dropping either feature set.
- Added generic two-port/N-port RF conversions, Touchstone SnP, real-reference renormalization, mixed-mode conversion, reference-plane shifts, passivity/reciprocity, K/mu stability, stability circles, transducer gain and load-pull sampling.
- Added microstrip, stripline, CPW, coax, TEM-line and rectangular-waveguide TE10 analytical/network baselines plus antenna-array/polarization utilities.
- Added thin-wire center-fed MoM and PEEC filament extraction baselines, including skin-depth/AC resistance and SPICE coupled-inductor export.
- Added a clean-room MNA circuit engine with R/C/L, independent/controlled sources, switches, mutual inductance, diode, MOS Level-1, BJT, JFET and sampled RF N-port devices.
- Added DC/AC/transient analyses, backward-Euler/trapezoidal/BDF2 integration, homotopy/PN limiting, sweeps, thermal noise integration, sensitivity, Monte Carlo, Fourier/THD and circuit S-parameter extraction.
- Extended SPICE parsing with engineering suffixes, `.MODEL`, `.PARAM`, expressions, nested `.SUBCKT`, scoped local models, `.FUNC`, `.INCLUDE` and `.LIB`.
- Added an OSDI/OpenVAF dynamic-library discovery/ABI-version seam while deliberately leaving descriptor evaluation and full DAE integration for the next compact-model phase.
- Added `docs/RF_CIRCUIT_RESEARCH.md` and expanded the authoritative tracker with Palace/OpenSEMBA/OpenNEC and ngspice/Xyce/QucsatorRF capability families.
- Final tracker state for this checkpoint is 429/576 validated capabilities (74.5%).
- GCC/OpenMP clean matrix passes 43/43 tests; focused Clang and serial RF/SPICE tests pass; focused ASan+UBSan runs pass for SPICE hierarchy/file IO, OSDI loading, PEEC/MoM and nonlinear circuit paths.

## 0.7.1 - CPML, TFSF and lumped FDTD elements

- Added a clean-room 1-D convolutional PML with graded sigma/kappa profiles and convolution memory variables.
- Added a +x 1-D total-field/scattered-field source with Yee temporal and half-cell spatial staggering.
- Added field-coupled parallel lumped R/L/C elements; resistor and capacitor contributions modify local Ampere coefficients while inductor current is advanced as an auxiliary state.
- Added CPML attenuation, homogeneous TFSF leakage, analytical inductor-ramp, resistor-dissipation and capacitor-stability regressions.
- Added CLI/CTest and sanitizer-smoke coverage for the new FDTD paths.
- Integration tracker advanced to 152/398 validated capabilities (38.2%), with openEMS-class FDTD at 17/26 (65.4%).

## 0.7.0 - dispersive FDTD and coupled-multiphysics foundation

- Added Debye, Drude and Lorentz 1-D FDTD auxiliary-differential-equation material models.
- Added 3-D PMC boundary baseline, TEM-like 1-D wave-port decomposition and complex S-parameter extraction.
- Added portable legacy ASCII VTK export for 3-D electric and magnetic fields.
- Added unit/location/topology-aware multiphysics field registry.
- Added conservative exact-overlap cell remapping and conservative face-flux-to-cell transfer.
- Added partitioned fixed-point coupling with Aitken relaxation.
- Added shared-mesh DC-conduction -> Joule-heating -> transient-FEM thermal coupling.
- Added user-visible FDTD/multiphysics CLI smoke cases and sanitizer coverage.
- Reworked the high-level roadmap to use the same phase numbering as the authoritative integration tracker.
- Upgraded `scripts/integration_status.py` with explicit tracker/manifest selection and JSON output.
- Integration tracker advanced to 149/398 validated capabilities (37.4%), with openEMS-class FDTD at 14/26 (53.8%) and coupled multiphysics at 6/15 (40.0%).

## 0.6.1 - nonlinear/adaptive FEM, porous flow, magnetostatics and modal analysis

- Added a reusable Newton solver with backtracking line search and ILU(0)-GMRES Jacobian solves.
- Added a nonlinear Tri3 Poisson/reaction manufactured-solution regression using the shared Newton stack.
- Added residual/flux-jump Poisson error estimation, Dorfler marking and conforming longest-edge Tri3 refinement.
- Added saturated isotropic Darcy porous-flow FEM with analytical channel validation.
- Added 2-D magnetostatics through the out-of-plane magnetic vector potential and reconstructed magnetic flux density.
- Added a sparse generalized eigenvalue solver using M-orthogonal deflated inverse iteration and shared PCG inner solves.
- Added a fixed-free Line2 bar modal/eigenfrequency application with consistent mass and analytical mode validation.
- Added CLI/CTest smoke cases and sanitizer coverage for the new FEM capabilities.
- Integration tracker advanced to 137/398 validated capabilities (34.4%), with Elmer-class FEM at 27/46 (58.7%).

## 0.6.0 - FEM/FDTD/optics breadth and release traceability

- Completed the Phase-0 release/benchmark traceability checklist with machine-readable benchmark NDJSON, history summarization, and a tag-driven artifact workflow.
- Added a reusable FEM reference-element catalogue for Line2, Tri3, Quad4, Tet4, Hex8, Prism6 and Pyramid5 with shape functions, gradients, quadrature, isoparametric mapping and Jacobians.
- Added reusable Tri3 CSR/matrix-free Laplace assembly and mixed Dirichlet/Neumann/Robin scalar-diffusion boundaries.
- Added 2-D Poisson, heat, linear elasticity, electrostatics and DC conduction applications using shared FEM infrastructure.
- Added axisymmetric linear elasticity with cylindrical weighting and hoop strain.
- Added a convergent 3-D Tet4 Poisson baseline.
- Added a 3-D Maxwell FDTD baseline plus 1-D material preprocessing, hard/soft sources, first-order Mur boundaries, time probes and DFT monitors.
- Added sequential and paraxial optics, Sellmeier materials, Jones/Stokes polarization, Fresnel coefficients and normal-incidence multilayer transfer matrices.
- Added CLI/CTest smoke cases for 3-D FEM, Mur FDTD and sequential/paraxial optics.
- Expanded sanitizer smoke coverage across the new FEM/FDTD/optics paths.
- Integration tracker advanced to 132/397 validated capabilities (33.2%).

## 0.5.0 - tracked integration, sparse numerics and electrochemistry foundation

- Added a 386-item checkable integration tracker spanning the OpenFOAM, FluidX3D, Elmer, openEMS, Optiland and chemistry/electrochemistry target capabilities.
- Added pinned upstream provenance manifest and integration-status reporter.
- Added CSR sparse matrix assembly/SpMV, Jacobi-PCG, BiCGStab, ILU(0) and restarted GMRES.
- Replaced the default collocated momentum fixed-point sweeps with an ILU(0)-GMRES full-momentum solve while retaining the legacy fallback for A/B validation.
- Added bounded limited-linear MUSCL/Barth-Jespersen-style finite-volume reconstruction.
- Added elementary Arrhenius/mass-action reaction networks.
- Added Nernst, Butler-Volmer and Faraday electrochemistry helpers.
- Added nonlinear 1-D corrosion current/recession solver with electrolyte ohmic feedback.
- Added conservative 1-D multi-species Nernst-Planck transport with diffusion, advection and electromigration.
- Added arbitrary-`PolyMesh` Nernst-Planck transport with electroneutral and Poisson electrostatic formulations.
- Added Scharfetter-Gummel exponential-fit drift-diffusion fluxes.
- Added Butler-Volmer/Faradaic reactive electrode species-flux boundaries and anodic dissolution coupling.
- Added a multi-reaction galvanic/mixed-potential solver.
- Added implicit generic `PolyMesh` scalar advection-diffusion transport through the shared CSR/ILU-GMRES layer.
- Added chemistry/electrochemistry research and validation documentation plus CLI/CTest smoke cases.


## 0.4.2 - Phase 4C collocated pressure/velocity coupling

### Added

- patch-based fixed-value/zero-gradient/slip velocity boundary conditions;
- fixed-value/zero-gradient pressure boundary conditions;
- orthogonal/non-orthogonal face-area decomposition;
- clean-room Rhie-Chow-style face-flux reconstruction using direct pressure differences;
- arbitrary-`PolyMesh` collocated momentum predictor retaining `HbyA` and `V/aP`;
- matrix-free non-orthogonal pressure correction using shared CG;
- SIMPLE, PISO and PIMPLE-style coupling loops;
- deterministic affine-sheared Cartesian-hexa mesh generator;
- checkerboard-pressure, sheared-mesh, Poiseuille-channel and lid-cavity regressions;
- `fvm-collocated-channel`, `fvm-collocated-cavity`, and `fvm-collocated-skew` CLI cases.

### Validation

- Rhie-Chow face flux responds to an alternating collocated pressure mode;
- sheared-mesh PISO reduces continuity L2 by more than three orders of magnitude and PIMPLE improves it further;
- 32x16 pressure-driven channel agrees with analytical Poiseuille profile within the 2% regression gate (about 0.5% in the development-container run);
- 20x20 moving-lid cavity develops the expected clockwise primary recirculation and conservative face flux;
- GCC/OpenMP local CTest matrix expanded to nine tests.

## 0.4.1 - Phase 4B transient incompressible baseline

### Added

- reusable linear and first-order upwind scalar/vector face interpolation on `PolyMesh`;
- conservative face-flux scalar/vector convective-divergence operators;
- periodic staggered transient incompressible Navier-Stokes solver with conservative momentum convection, centered viscosity and shared-CG pressure projection;
- Taylor-Green analytical velocity-error measurement and refinement regression;
- `fvm-taylor-green2d` CLI case;
- sanitizer coverage for face schemes and the transient solver.

### Validation

- constant upwind transport and uniform momentum conservation regressions pass;
- linear scalar advection is exact to floating-point tolerance on the Cartesian-hexa validation mesh using linear face interpolation;
- 24x24 -> 48x48 Taylor-Green refinement reduces velocity error by the required regression factor;
- the 64x64 CLI case at `t=0.01` reports divergence near machine precision and velocity RMS error around `1.5e-7` in the development-container run.

## 0.4.0 - Phase 4A polyhedral FVM foundation

### Added

- clean-room owner/neighbour `PolyMesh` topology with boundary patches and oriented face-area vectors;
- Cartesian hexahedral mesh generator for validation;
- OpenMP cell-parallel Gauss scalar gradient and vector divergence;
- orthogonal two-point scalar Laplacian;
- reusable matrix-free conjugate-gradient core solver;
- periodic staggered finite-volume pressure projection using CG;
- manufactured FVM operator regressions and projection regressions;
- `fvm-operators3d` and `fvm-projection2d` CLI cases;
- sanitizer smoke coverage for the new mesh/operator/projection code.

### Validation

- linear Gauss gradient and linear-vector divergence are exact to floating-point tolerance on the Cartesian hexa validation mesh;
- interior quadratic Laplacian is exact to floating-point tolerance;
- the 96x80 projection CLI reduces divergence RMS from approximately `3.14e-1` to `2.85e-11` and converges its pressure solve in 7 CG iterations in the development-container case.

## 0.3.1 - Phase 3B device-resident distributed transport

### Added

- D3Q19/D3Q27 selective face/edge/corner halo plans that transmit only crossing populations.
- Shared `constexpr` host/device halo linearization.
- Persistent compact MPI point-to-point requests.
- Device-resident distributed D3Q19/D3Q27 SYCL two-grid pull reference block.
- Device halo pack/unpack kernels.
- Pinned-host SYCL USM staging transport for ordinary MPI.
- Explicit direct device-buffer `--gpu-aware-mpi` transport.
- MPI+SYCL staged distributed parity tests when both backends are enabled.
- `cfd-halo-plan` exact halo-volume reporting CLI.
- `scripts/benchmark_distributed.sh` strong/weak-scaling driver with machine-readable CSV output.
- aggregate halo message/byte counters in `cfd-distributed` benchmark reporting.
- D3Q19/D3Q27 selection and CPU/SYCL backend selection in `cfd-distributed`.
- Distributed SYCL checkpoint/restart using the Phase-3 rank-local format.

### Changed

- The CPU MPI CLI now uses the selective persistent halo exchange rather than the Phase-3A full-cell nonpersistent baseline.
- CMake project version is 0.3.1.

### Validation

- Selective virtual D3Q19/D3Q27 exchange is bitwise-identical to full-cell halo exchange on the uneven 17x15x13 2x2x2 decomposition.
- A 64^3 brick reduces exact FP32 halo payload by 74.33% for D3Q19 and 67.35% for D3Q27.
- GCC/OpenMP regressions pass.
- MPI-only and combined MPI+SYCL Phase-3B translation units pass local API syntax checks.
- Real MPI+GPU execution remains a hardware/CI validation item because this development environment has neither MPI nor AdaptiveCpp/GPU runtimes.
## 0.3.0 - Phase 3A distributed CPU/hybrid runtime

### Added

- Aspect-aware Cartesian process-grid selection and uneven brick decomposition.
- 26-neighbour face/edge/corner haloed SoA storage and virtual-rank exchange.
- Distributed D3Q19/D3Q27 two-grid pull blocks with Guo forcing.
- OpenMP collision work inside each distributed rank.
- Strict-interior / boundary-shell split for communication overlap.
- MPI Cartesian runtime, node-local rank discovery, and nonblocking `MPI_Irecv`/`MPI_Isend` halo exchange.
- Deterministic local-rank -> visible-SYCL-accelerator mapping.
- Rank-local versioned binary checkpoint/restart.
- `cfd-distributed` MPI CLI.
- Four-rank OpenMPI CI job.
- Distributed runtime and Phase-3 result documentation.

### Validation

- D3Q19 and D3Q27 2x2x2 virtual-distributed runs match the monolithic reference on an uneven 17x15x13 periodic domain.
- Face, edge, and corner ghost regions are exercised by the parity tests.
- Checkpoint reload is exact and continuation remains bitwise-identical.
- Distributed virtual halo code is covered by the ASan+UBSan smoke target.
- MPI translation units are syntax-checked locally; real MPI execution remains a CI/hardware validation item because MPI is not installed in this development container.

## 0.2.1 - Phase 2B physical LBM validation

### Added

- Guo-style body acceleration with half-step velocity correction in D2Q9, generic two-grid pull, single-grid CPU, and optional single-grid SYCL kernels.
- Halfway bounce-back stationary walls for the D2Q9 CPU reference solver.
- Moving-wall bounce-back correction and lid-driven cavity case.
- Zou-He left velocity inlet and right pressure/density outlet.
- Force-driven Poiseuille analytical regression.
- Velocity-inlet/pressure-outlet channel regression.
- Taylor-Green analytical decay/convergence regression.
- Generic templated FP32/FP64 periodic pull reference solver.
- FP32/FP64 solution-parity regression.
- `lbm-poiseuille`, `lbm-cavity`, and `lbm-channel-io` CLI cases.
- `docs/LBM_BOUNDARIES.md` and `docs/PHASE2B_RESULTS.md`.

### Validation

- D2Q9/D3Q19/D3Q27 forced in-place/two-grid parity: pass.
- Poiseuille analytical velocity profile: pass.
- Lid-driven cavity vortex-sign regression: pass.
- Zou-He inlet/outlet throughput and outlet-density checks: pass.
- Taylor-Green refinement trend: pass.
- FP32/FP64 maximum field difference gate: pass.
- Stationary-wall API reset semantics regression: pass.
- Dedicated ASan+UBSan boundary/precision smoke target: pass.
- Optional SYCL force path and forced CPU/SYCL parity tests implemented; hardware validation remains pending.

## 0.2.0 - Phase 2A fluid performance baseline

### Added

- D2Q9, D3Q19 and D3Q27 compile-time lattice descriptors with opposite-direction pairing.
- Clean-room Esoteric-Pull single-grid BGK solver for D2Q9/D3Q19/D3Q27.
- Generic two-grid pull solver for numerical/performance comparison.
- Optional AdaptiveCpp/SYCL single-grid kernels for D2Q9/D3Q19/D3Q27.
- CPU/SYCL parity tests gated by `CFD_ENABLE_SYCL`.
- D3Q19/D3Q27 descriptor-isotropy, mass, uniform-flow and Taylor-Green regression tests.
- `cfd-bench` with MLUPS, estimated DDF traffic, allocated population memory and relative mass drift.
- Side-by-side `--streaming both` benchmark mode.
- Esoteric-Pull implementation notes, GPU backend strategy and Phase-2 validation report.

### Changed

- Existing SYCL D2Q9 queue is explicitly in-order so USM kernel dependencies do not rely on implementation scheduling.
- CLI exposes D2Q9 single-grid, D3Q19 and D3Q27 examples.
- CMake project version is now 0.2.0.
