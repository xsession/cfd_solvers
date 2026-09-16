# Unreleased CPU continuation

This work follows v0.6.1 and advances all existing solver families. It is not a sealed release.
The tracker reports **165/402 capabilities (41.0%)**, up from the workspace's 149/398 baseline
(which already included work beyond the historical README count). Four explicit baseline entries
were added rather than treating narrow implementations as complete production capabilities.

## APIs and limits

| Area | Public API | Scope |
|---|---|---|
| LBM | `OneStepPullSolver(config, odd_tau)` | Periodic CPU TRT for D2Q9/D3Q19/D3Q27. `config.tau` controls viscosity/even relaxation; `odd_tau` controls odd relaxation. Equal times preserve the BGK path. Guo forcing is split consistently. In-place, boundary and GPU TRT remain open. |
| FVM | `ScalarSchemeWorkspace`, scalar `*_into` gradient/interpolation/divergence | Caller-owned output and reusable scratch; no allocation after sizing on a fixed mesh. Output/input spans must not overlap. Vector interpolation and collocated pressure/momentum still allocate. |
| FVM reconstruction | `minmod`, `van_leer`, least-squares gradient, corrected normal gradient/Laplacian | Linear-field and flux-balance regressions. TVD reconstruction follows the upwind/downwind centre line and clamps to endpoint bounds. Arbitrary skew-mesh multidimensional TVD time integration is not claimed. |
| FVM transport | `ScalarTransport::operator_assemblies()` | Retains matrix, ILU and RHS buffers until boundary/face-flux changes. Source and initial-value changes retain the operator. GMRES still owns temporary Krylov storage. |
| Reaction rates | `ElementaryReaction::equilibrium`, `EquilibriumConstant`, buffer overloads | Forward-minus-reverse mass action, with `kr=kf/Kc`; constant-enthalpy van't Hoff dependence. Kc must use the same concentration units/powers as the reaction, not an unconverted activity constant. |
| Stiff chemistry | `integrate_isothermal` | Small constant-volume/isothermal networks; backward Euler with two half steps, step-doubling error control, finite-difference Newton Jacobian and a nonnegative line search. Caller state is unchanged if integration fails. No energy equation or large-mechanism sparse Jacobian yet. |
| Aqueous chemistry | `equilibrate_acids` | Ideal dilute acid/base families, water and spectator charge. Log-H charge balance and normalized mass-action distributions. This API uses mol/L; transport uses mol/m^3, so concentration transfer requires multiplication by 1000. No activity corrections, precipitation or automatic electrode/transport coupling. |
| FDTD | `Boundary1D::pml`, Debye/Drude/Lorentz setters | Polynomial matched electric/magnetic loss layer in 1-D, local CFL and conservative explicit ADE-pole limits. Changed/dispersive background material cannot overlap the layer. Multidimensional UPML/CPML remains open. |
| Optics | `SurfaceType::conic`, `even_asphere`, `surface_sag` | Rotational conic sag and r^4/r^6/r^8/r^10 corrections with Newton ray intersection. Misses/nonconvergence return invalid traces. No local-coordinate transforms or non-sequential tracing yet. |
| Gaussian optics | `GaussianBeam` | SI lengths, q=z+i*z_R and 1/e field radius. Propagation, thin lens and real ABCD matrices with determinant n_in/n_out. No sampled diffraction/PSF/MTF. |
| FEM/coupling | `Elasticity2D::solve_thermal`, `JouleHeatingCoupler2D` | Shared Tri3 DC conduction -> Joule heat -> transient heat -> small-strain expansion. Plane strain includes the out-of-plane thermal constraint. No electrical temperature feedback or mesh deformation. |

The 1-D FDTD magnetic field now follows the physical convention `Hy=-Ez/Z` for +x propagation,
matching the port and 3-D convention. Before wave decomposition, average E over the previous/current
times and H over adjacent spatial edges to collocate the staggered fields. `energy()` reports only
electromagnetic field energy, excluding ADE oscillator energy. Gaussian initialization resets time,
H and polarization. The PML target reflection controls the continuum profile; it is not a guarantee
of discrete reflected amplitude.

The FEM heat mass term now uses old-time boundary values; the new-time values only enter the
implicit boundary contribution. The Joule coupler rejects unsolved/stale electrical fields, and
the fixed-point driver checks the residual after its final allowed update.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFD_ENABLE_OPENMP=ON
cmake --build build --parallel
OMP_NUM_THREADS=2 ctest --test-dir build --output-on-failure
OMP_NUM_THREADS=2 ./build/cfd-solve lbm-trt
OMP_NUM_THREADS=2 ./build/cfd-solve fvm-workspace
./build/cfd-solve fdtd-pml1d
./build/cfd-solve chemistry-reactor
./build/cfd-solve chemistry-equilibrium
./build/cfd-solve optics-gaussian
./build/cfd-solve multiphysics-thermoelastic
python3 scripts/integration_status.py
```

`cfd-extension-tests` checks analytical ADE response and timestep refinement, pulse absorption,
physical wave direction, sheared-mesh FVM operators, TVD bounds, conservation, workspace parity,
cache invalidation, reversible/nonlinear stiff kinetics, reaction mass balance, aqueous charge
balance, Gaussian focusing, conic intersections and electro-thermoelastic patch solutions.
It is included in sanitizer CI.

## Local measurements

2026-09-16: Ubuntu 22.04 under WSL, Intel Core i7-8750H (12 logical CPUs visible), GCC 11.4,
Release/native architecture, two OpenMP threads, double precision FVM. Timings are local samples.

| Case | Allocating/rebuilding reference | Reused workspace/operator | Maximum numerical difference |
|---|---:|---:|---:|
| Scalar bounded-linear convection, 6144 cells, 100 repetitions | 0.136372 s | 0.0575047 s | 0 |
| Scalar diffusion, 256 cells, 50 steps | 0.0062817 s | 0.0024418 s | 0 |

The scalar workspace retains 452096 bytes with unchanged capacities/storage across repeated calls;
flux-balance error was 1.11e-16. The transport cache assembles/factorizes once rather than 50 times.
These observations do not establish process-wide peak-memory reduction or portable speedup.

Other case outputs: PML residual field-energy ratio 3.24e-8 after 800 steps; reversible reactor
A/B concentrations 0.2/0.8; acid buffer pH 5.00017 and charge residual 6.25e-17 mol/L; Gaussian
focus 0.0999747 m and waist 1.59135e-5 m; thermal-expansion displacement error 3.52e-10.

## Verification and remaining work

GCC/OpenMP and GCC serial each pass 34 CTest targets. MPICH adds three four-rank tests (37 total).
The new extension tests and existing sanitizer smoke pass GCC ASan+UBSan with leak detection.
Local Clang is unavailable and package installation failed on WSL DNS resolution; the Clang CI
configuration remains. AdaptiveCpp/GPU execution is unvalidated.

Use the installed MPICH wrappers explicitly if the default OpenMPI development setup is incomplete:

```bash
cmake -S . -B build-mpi -DCMAKE_BUILD_TYPE=Release -DCFD_ENABLE_MPI=ON \
  -DMPI_CXX_COMPILER=/usr/bin/mpicxx.mpich -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec.mpich
cmake --build build-mpi --parallel
OMP_NUM_THREADS=2 ctest --test-dir build-mpi --output-on-failure
```

Remaining major areas include LBM MRT/thermal/multiphase; FVM turbulence, multiphase, compressible
flow and full workspace adoption; nonlinear/mixed/distributed FEM; 3-D CPML/TFSF/GPU FDTD;
sampled optics and optimization; nonideal aqueous chemistry and implicit electrode/transport
coupling; cross-mesh multiphysics; portable case formats and interoperability. The unchecked
integration tracker remains the completion criterion.

## Mathematical references

Implementations are independent; no upstream solver code was copied.

- [Dellar, TRT LBM](https://people.maths.ox.ac.uk/dellar/papers/MHD_magic_TRT.pdf): even/odd collision decomposition.
- [Cantera reaction rates](https://www.cantera.org/dev/reference/kinetics/reaction-rates.html): mass-action conventions.
- [USGS equilibrium constraints](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc.v1/html/phqc_11.htm) and [charge balance](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc.v1/html/phqc_12.htm): the current implementation covers only ideal acid/base families.
- [Meep materials](https://meep.readthedocs.io/en/latest/Materials/) and [matched layers](https://meep.readthedocs.io/en/latest/Perfectly_Matched_Layer/): ADE equations, stability considerations and matched electric/magnetic loss.
- [Gaussian-beam ABCD propagation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10581742/): complex-q transform convention.
