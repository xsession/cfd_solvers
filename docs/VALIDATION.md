# Validation

The regression suite preserves the Phase-2/3 LBM and distributed gates and now also validates the Phase-4A finite-volume foundation.

## Numerical regression

`cfd-tests` checks:

- D2Q9, D3Q19 and D3Q27 descriptor weights, opposite-pair symmetry and second-order isotropy moments;
- original two-grid D2Q9 uniform-flow invariance and mass conservation;
- single-grid Esoteric-Pull D2Q9 against the original two-grid solver after the same Taylor-Green initialization and time steps;
- D3Q19 uniform 3-D flow invariance and mass conservation;
- D3Q27 Taylor-Green viscous energy decay and mass conservation;
- Guo-force single-grid/two-grid macroscopic parity for D2Q9, D3Q19 and D3Q27;
- Poiseuille analytical velocity-profile agreement;
- lid-driven cavity vortex-direction and stability checks;
- Zou-He velocity-inlet / pressure-outlet throughput and outlet-density checks;
- Taylor-Green analytical decay/refinement convergence;
- FP32 versus FP64 low-Mach field parity;
- stationary-wall API reset semantics after a cell had previously been a moving wall;
- CPU/SYCL forced macroscopic parity for D2Q9/D3Q19/D3Q27 whenever `CFD_ENABLE_SYCL=ON`;
- the Phase-1 FVM, FEM, FDTD and optics regressions;
- aspect-aware Cartesian decomposition coverage and periodic/nonperiodic neighbour mapping;
- D3Q19/D3Q27 2x2x2 virtual-distributed parity against the monolithic reference on an uneven 17x15x13 domain;
- exact checkpoint reload and bitwise-identical continuation after restart;
- selective D3Q19/D3Q27 halo counts and exact byte accounting;
- selective virtual exchange bitwise parity against the Phase-3A full-cell exchange.
- shared matrix-free CG on an independent SPD diagonal system;
- `PolyMesh` owner/neighbour topology and boundary counts on a generated hexahedral mesh;
- exact linear scalar Gauss gradient and exact linear-vector divergence to floating-point tolerance;
- exact interior quadratic Laplacian on the orthogonal validation mesh;
- periodic staggered pressure projection reducing deliberately injected discrete divergence;
- projected Taylor-Green sampled-face divergence reduction with bounded kinetic-energy change.
- `PolyMesh` linear/upwind face interpolation and conservative convection invariants;
- transient Taylor-Green incompressible integration, pressure convergence and mesh/time-step refinement trend;
- collocated Rhie-Chow response to an alternating checkerboard pressure field;
- non-orthogonal PISO/PIMPLE continuity correction on a deterministic sheared mesh;
- collocated SIMPLE pressure-driven channel agreement with the analytical Poiseuille profile;
- collocated SIMPLE lid-driven-cavity recirculation and continuity gates.
- FEM reference-element partition-of-unity, gradient-sum and reference-measure quadrature checks for Line2/Tri3/Quad4/Tet4/Hex8/Prism6/Pyramid5;
- Tri3 assembled-CSR versus matrix-free Laplace action parity;
- exact mixed Dirichlet/Neumann/Robin scalar-FEM manufactured solution;
- electrostatic parallel-plate and DC-conduction uniform-field/current checks;
- axisymmetric-elasticity uniform-dilatation patch test including hoop strain;
- 3-D Tet4 Poisson manufactured-solution mesh-refinement convergence;
- heterogeneous/lossy 1-D FDTD material update, hard/soft source semantics and first-order Mur reflection suppression;
- time-probe sample counts and DFT amplitude recovery;
- 3-D Maxwell PEC stability/finite-field baseline;
- Sellmeier dispersion, Jones/Stokes conversion, Fresnel normal-incidence reflectance and quarter-wave thin-film antireflection checks;
- sequential real-ray and paraxial focal-distance/spot regressions;
- benchmark NDJSON schema smoke validation.

CPU/SYCL tests are deliberately part of the normal test executable rather than a separate benchmark so accelerator changes cannot bypass correctness gates.

## Build/toolchain validation performed for this phase

Validated in the development container on 2026-09-16:

- GCC 14.2 Release + OpenMP build + all **18** local CTest targets: pass;
- Clang 17 Release serial-fallback build + all **18** local CTest targets: pass;
- GCC 14.2 Release serial build (`CFD_ENABLE_OPENMP=OFF`) + all **18** local CTest targets: pass;
- GCC AddressSanitizer + UndefinedBehaviorSanitizer serial build: dedicated `cfd-sanitize-smoke` passes with FEM, FDTD, optics, FVM and electrochemistry coverage. The long physical-regression executable remains outside the sanitizer gate because the O0-instrumented workload is disproportionately slow; Release builds carry the full numerical suite.

Neither MPI nor an AdaptiveCpp/SYCL runtime is installed in the development container. MPI/SYCL source paths remain covered by their existing syntax/CI wiring, while real multi-rank and accelerator execution remains a hardware validation item.

Each release patch is verified with `git apply --check`, `git diff --check`, a clean rebuild, and CTest after application to the previous source release. v0.6.0 uses the sealed v0.5.0 source release as its patch/rebuild baseline.
The strong/weak scaling shell driver is also smoke-tested with mocked MPI/solver executables to validate its CSV schema and weak-axis dimension scaling without requiring MPI locally.

## Development-container benchmark snapshot

Hardware visible to the container: 5 vCPUs from an Intel Xeon Platinum 8370C VM. GCC 14.2, Release, `-march=native`, OpenMP enabled. These numbers are a smoke/performance-regression baseline, not a hardware comparison claim.

| Lattice | Grid | Threads | MLUPS | Estimated DDF GB/s | Single-grid DDF | Two-grid DDF |
|---|---:|---:|---:|---:|---:|---:|
| D3Q19 | 64^3 | 1 | 6.83 | 1.04 | 19.0 MiB | 38.0 MiB |
| D3Q19 | 64^3 | 2 | 11.15 | 1.70 | 19.0 MiB | 38.0 MiB |
| D3Q19 | 64^3 | 5 | 16.15 | 2.45 | 19.0 MiB | 38.0 MiB |
| D3Q27 | 48^3 | 1 | 4.12 | 0.89 | 11.39 MiB | 22.78 MiB |
| D3Q27 | 48^3 | 5 | 8.60 | 1.86 | 11.39 MiB | 22.78 MiB |

`estimated DDF GB/s` counts one load and one store of every FP32 population per lattice update; it is a model of DDF traffic, not a hardware-counter measurement.

Reproduce with, for example:

```bash
./build/cfd-bench --lattice d3q19 --nx 64 --ny 64 --nz 64 \
  --warmup 10 --steps 80 --threads 5
```

## Phase 4A FVM regression snapshot

On the generated orthogonal hexahedral regression mesh, the linear scalar gradient, linear-vector divergence and interior quadratic Laplacian match their manufactured exact values to the configured floating-point tolerances.

For the 96x80 periodic projection CLI case in the development container:

- divergence RMS before projection: about `3.14091e-1`;
- divergence RMS after projection: about `2.85e-11`;
- matrix-free pressure CG iterations: `7`;
- reported pressure residual RMS: about `2.85e-8`.

These values verify the discrete projection implementation on that case; they are not a claim of general Navier-Stokes accuracy or a portable performance result.

## Phase 4B transient-flow regression snapshot

At physical time `t=0.005`, the Taylor-Green refinement test observed:

- 24x24, `dt=1e-4`: velocity RMS error `5.56e-7`, divergence RMS `7.33e-15`;
- 48x48, `dt=2.5e-5`: velocity RMS error `1.39e-7`, divergence RMS `1.67e-15`;
- fine/coarse velocity-error ratio: about `0.250`.

The normal 64x64 CLI case at `t=0.01` reports velocity RMS error about `1.50e-7` and divergence near machine precision. Exact pressure-CG iteration counts can vary slightly with reduction ordering and thread count, so convergence is gated rather than a fixed iteration count.


## Phase 4C collocated-flow regression snapshot

Development-container examples on the small regression cases:

- 32x16 pressure-driven channel: profile relative L2 error about `4.7e-3`, continuity L2 about `2.2e-8`, pressure CG converged;
- 20x20 moving-lid cavity: centre `u_x` about `-0.191`, continuity L2 order `1e-7`;
- 12x10x1 xy-sheared mesh: continuity L2 falls from about `3.01e-1` to `1.20e-4` after one PISO step and about `2.00e-5` with two PIMPLE outer correctors.

These cases validate the current pressure/velocity coupling and mesh correction. They are not high-Re benchmark claims or full OpenFOAM-equivalence claims.

## v0.6.1 FEM advanced gates

- nonlinear Tri3 Poisson/reaction manufactured solution converges through shared Newton + ILU(0)-GMRES and remains below the FEM error gate;
- residual/jump estimator + Dorfler marking selects a strict subset on the manufactured Poisson case, conforming longest-edge refinement validates, and refined error is lower than the coarse error;
- saturated Darcy channel reproduces the analytical streamwise velocity with negligible transverse component;
- 2-D magnetostatic vector-potential manufactured solution and reconstructed magnetic flux are finite and convergent;
- fixed-free Line2 generalized eigenmodes agree with analytical axial-bar eigenfrequencies;
- dedicated CLI smoke tests cover nonlinear FEM, Darcy, magnetostatics and modal analysis;
- the sanitizer smoke target exercises nonlinear solve, adaptive mesh generation, Darcy, magnetostatics and the generalized eigen path on small meshes.

## v0.7.0 FDTD and multiphysics gates

Clean v0.7.0 source configurations define **26 CTest tests** with MPI/SYCL disabled. The complete CPU suite passes under GCC 14 + OpenMP, Clang 17 without OpenMP, and serial GCC; the dedicated ASan+UBSan smoke target also passes.

- Debye, Drude and Lorentz ADE material regions remain finite under the regression pulses and produce nonzero polarization response; replacing a dispersive region with `set_material()` clears its ADE state.
- 3-D PMC baseline clamps the magnetic boundary state and remains finite.
- Synthetic TEM-like port data recovers `|S11|=0.20` and `|S21|=0.75` in the CLI smoke case (unit test also checks complex amplitudes against known values).
- Maxwell3D legacy ASCII VTK output is written and parsed back for expected geometry/E/H field sections.
- Exact-overlap 1-D cell transfer preserves the integrated scalar to floating-point tolerance.
- `PolyMesh` face-flux-to-cell transfer cancels internal flux exactly in the volume-weighted global sum.
- Aitken partitioned coupling reaches the cosine fixed point in six iterations in the development-container smoke case.
- Shared-mesh DC conduction with 4 V across 2 m and `sigma=5 S/m` reproduces `q=20 W/m^3`; the transient heat solve remains converged and raises the interior temperature above the fixed 300 K boundaries.
- `scripts/integration_status.py --json` emits a machine-readable report with the same completion counts as the Markdown tracker.

## v0.7.0 build matrix

Validated on 2026-09-16 in the development container:

- GCC 14.2 Release + OpenMP: full regression executable passes; 26 CTest entries are configured.
- Clang 17 Release with `CFD_ENABLE_OPENMP=OFF`: full regression executable and sanitizer-smoke executable pass.
- GCC 14.2 Release serial (`CFD_ENABLE_OPENMP=OFF`): full regression executable, sanitizer smoke, dispersive-FDTD and electrothermal CLI cases pass.
- GCC ASan + UBSan Debug serial build: `cfd-sanitize-smoke` passes with the new dispersive-FDTD and multiphysics paths.
- `git diff --check`: clean before release freeze.

The OpenMP small-loop cutoff was separately microbenchmarked on the 5-vCPU development runner. Forced OpenMP was slower than serial execution below roughly 2k trivial iterations in that environment, while becoming beneficial at and above the selected 2048-iteration threshold. With the cutoff/nested-team guard enabled, the full OpenMP regression executable completed in about 2.4 s in the measured run, versus the earlier order-of-13--15 s baseline. This is a development-runner observation, not a portable performance guarantee.
## v0.7.1 advanced FDTD gates

A clean v0.7.1 CPU configuration defines **28 CTest tests** with MPI/SYCL disabled.

- 1-D CPML reduces the residual Gaussian-pulse field energy below `1e-4` of its initial value and more than four orders of magnitude below the PEC residual in the regression case.
- The homogeneous +x TFSF regression launches a transmitted peak above `0.5` while keeping the sampled scattered-field-side leakage below `1e-4` of that peak; the development run is approximately `5.2e-6`.
- A parallel lumped inductor driven by a repeatedly imposed constant cell voltage reproduces the analytical current-density ramp to floating-point tolerance.
- A parallel lumped resistor dissipates more than 10% of the corresponding reference field energy over the smoke interval; a capacitor-only element remains finite.

## v0.8.0 RF/SPICE gates

Validated on 2026-09-16 in the development container:

- GCC 14/OpenMP clean configuration passes **43/43 CTest targets**.
- Focused Clang 17 and serial-GCC RF/SPICE builds pass the RF/network, PEEC/MoM, circuit-analysis, hierarchy/file-parser and OSDI seams.
- SPICE regressions cover linear/nonlinear DC, AC, BE/trapezoidal/BDF2 transient, controlled sources, switches, mutual inductance, diode/MOS/BJT/JFET models, sampled N-port stamping, source/parameter/temperature sweeps, resistor noise, sensitivity, Monte Carlo, Fourier/THD, `.PARAM`, expressions, nested `.SUBCKT`, `.FUNC`, `.INCLUDE`, `.LIB` and scoped models.
- RF regressions cover Z/S/ABCD and generic N-port round-trips, Touchstone, de-embedding, mixed-mode conversion, passivity/reciprocity, return loss/VSWR, K/mu, stability circles, transducer gain/load-pull, transmission-line/waveguide baselines, array/polarization utilities, PEEC reciprocity/skin-effect behavior and center-fed thin-wire MoM finite input impedance/current distribution.
- OSDI regression loads a synthetic shared library, validates ABI version/count metadata and move semantics; descriptor evaluation is intentionally still open.
- Focused ASan+UBSan with leak detection passes SPICE file/hierarchy parsing, OSDI loading and a PEEC + thin-wire-MoM + nonlinear circuit smoke.
- The broad sanitizer build is intentionally not used as the RF/SPICE gate because the monolithic project build produces excessive debug/instrumentation compile cost; focused instrumentation covers the newly changed paths directly.

## v0.8.1 adaptive-circuit and coupled-wire gates

Validated on 2026-09-17 in the development container. The release CPU configuration passes **47/47 CTest targets**, including the three new CLI smoke cases and the focused `cfd-next-analysis-tests` target.

- Adaptive RC transient regression begins from a deliberately coarse step, exercises at least one LTE rejection, accepts a genuinely nonuniform step sequence, lands exactly on the requested stop time and agrees with the analytical first-order step response.
- Periodic steady-state RC regression converges by phase-aligned cycle residual and its extracted fundamental amplitude agrees with the small-signal RC transfer magnitude.
- Low-pass RC pole fitting recovers the analytical pole at `s=-1/(RC)` with no finite zero; high-pass RC fitting additionally recovers the zero at the origin.
- Two identical separated parallel wires produce nonzero induced current and the full current solution is reciprocal under swapping the driven wire.
- A 1 Mohm center-segment series load reduces the driven-wire feed current in the coupled-wire regression.
- Adaptive BDF2 is validated with a uniform-step refinement gate: halving the accepted step reduces RC endpoint error by more than the second-order regression threshold, while the implementation uses unequal-step coefficients in variable-step operation.
- The arbitrary-orientation MoM kernel is validated by rigidly rotating the reciprocal two-wire geometry onto a non-axis-aligned 3-D direction; the entire current vector remains invariant to numerical precision.
- Crossing/touching oriented wires are rejected explicitly until connected-junction basis functions are available.
- The ideal transformer regression verifies a 2:1 voltage ratio and the expected 4x impedance reflection using a 100 ohm secondary load; the tolerance includes the configured MNA `gmin`.
- A matched 50 ohm sampled TEM line with `alpha=0.1 Np/m` and `L=1 m` reproduces `|S21|=exp(-0.1)` with negligible `S11`.
- The focused target is `cfd-next-analysis-tests`; these checks intentionally do not claim connected-wire junction basis functions, NEC-grade canonical convergence, harmonic balance, exact descriptor pole-zero extraction or BDF3+ automatic order selection.
