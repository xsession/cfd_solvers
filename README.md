# cfd_solvers

A clean-room, C++20-first collection of numerical solvers with one shared high-performance runtime for multicore CPUs and heterogeneous GPUs.

The long-term goal is a readable solver collection covering CFD, FEM multiphysics, FDTD electromagnetics and optics without forcing those different numerical methods into one algorithm. They share memory layouts, execution backends, decomposition, linear algebra, geometry, materials, I/O and diagnostics.

The project is inspired by the capabilities and engineering lessons of OpenFOAM, FluidX3D, Elmer FEM, openEMS FDTD and Optiland. It is **not** a source-code merge or mechanical translation. The source projects have different licenses, and FluidX3D has additional restrictions, so performance techniques are independently implemented from publications and public descriptions.

## v0.15.2 status - transported RANS, Reynolds stress and SST-DES

The machine-checkable integration program now reports **608/719 capabilities (84.6%)** complete. Phase 1 common HPC runtime remains **34/34 (100.0%)**, Phase 2 FluidX3D-class LBM remains **43/46 (93.5%)**, Phase 3 OpenFOAM-class FVM advances to **94/106 (88.7%)**, and Phase 9 workflow/UX remains **50/50 (100.0%)**. The default CPU regression matrix contains **122 CTest targets**.

This checkpoint closes the five remaining turbulence-transport tracker items. `SpalartAllmarasTransport`, `KEpsilonTransport` and `KOmegaSSTTransport` add transported turbulence state on the polyhedral FVM mesh with implicit upwind advection, variable diffusion, semi-implicit destruction, positivity bounds, fixed-value/zero-gradient patches and Euler/BDF2/Crank-Nicolson time integration. SST includes F1/F2 blending, cross diffusion, production limiting and an integrated DES dissipation-length switch with optional F1/F2 zonal shielding.

`ReynoldsStressTransport` adds six symmetric Reynolds-stress equations plus an epsilon equation, LRR-style pressure-strain redistribution, production from the supplied velocity gradient, a turbulent-diffusion closure, optional source injection and a realizability projection that keeps the transported covariance tensor positive semidefinite.

Validation passes **122/122** CPU tests. The focused v0.15.2 regression covers analytic linear shear, SA growth, k-epsilon decay/production, SST blending, SST-DES switching, Reynolds-stress anisotropy/realizability and fixed-boundary turbulent diffusion. Clean-room documentation is in `docs/UPSTREAM_DOCUMENTATION_REVIEW_0_15_2.md`.

The clean-room rule remains unchanged: upstream projects and commercial-product capability descriptions are architecture/research references only; incompatible source code is not translated or copied into the MIT core.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFD_ENABLE_OPENMP=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/cfd-solve --list
```

Examples:

```bash
./build/cfd-solve lbm-d2q9-inplace-cpu
./build/cfd-solve lbm-d3q19-cpu
./build/cfd-solve lbm-d3q27-cpu
./build/cfd-solve lbm-poiseuille
./build/cfd-solve lbm-cavity
./build/cfd-solve lbm-channel-io
./build/cfd-solve fvm-operators3d
./build/cfd-solve fvm-projection2d
./build/cfd-solve fvm-taylor-green2d
./build/cfd-solve fvm-collocated-channel
./build/cfd-solve fvm-collocated-cavity
./build/cfd-solve fvm-collocated-skew
./build/cfd-solve fvm-scalar-transport
./build/cfd-solve fem-poisson3d
./build/cfd-solve fem-nonlinear-poisson
./build/cfd-solve fem-darcy
./build/cfd-solve fem-magnetostatic
./build/cfd-solve fem-modal-bar
./build/cfd-solve fdtd-mur1d
./build/cfd-solve fdtd-dispersive1d
./build/cfd-solve fdtd-port-vtk
./build/cfd-solve fdtd-cpml-tfsf
./build/cfd-solve fdtd-lumped-rlc
./build/cfd-solve multiphysics-coupling
./build/cfd-solve multiphysics-electrothermal
./build/cfd-solve optics-lens
./build/cfd-solve electrochem-corrosion1d
./build/cfd-solve electrochem-pnp1d
./build/cfd-solve electrochem-pnp-poly
./build/cfd-solve electrochem-galvanic
./build/cfd-solve em-edge3d
./build/cfd-solve em-waveport
./build/cfd-solve particle-em-pic1d
./build/cfd-solve particle-pic2d
./build/cfd-solve particle-em-pic2d
./build/cfd-solve particle-staggered-em-pic2d
./build/cfd-solve particle-pic3d
./build/cfd-solve particle-em-pic3d
./build/cfd-solve particle-distributed-round3d
./build/cfd-solve particle-field-guards3d
./build/cfd-solve particle-distributed-step3d
./build/cfd-solve particle-geant4-transport
./build/cfd-solve particle-transport-dose-bvh
./build/cfd-solve particle-campaign-doe
./build/cfd-solve particle-campaign-execution
./build/cfd-solve particle-multiserver-deploy
./build/cfd-solve particle-multiserver-execution
./build/cfd-solve particle-multiserver-supervision
./scripts/integration_status.py
```

## Benchmark

Compare the single-grid scheme with the conventional two-grid pull baseline:

```bash
./build/cfd-bench --lattice d3q19 --nx 64 --ny 64 --nz 64 \
  --warmup 10 --steps 100 --threads 8 --streaming both
```

CSV output:

```bash
./build/cfd-bench --lattice d3q27 --nx 64 --ny 64 --nz 64 \
  --warmup 10 --steps 100 --threads 8 --streaming both --csv
```

Or run the standard CPU set:

```bash
CFD_BENCH_THREADS=8 ./scripts/benchmark_lbm.sh build
```

`estimated_DDF_GB/s` is an algorithmic minimum based on one read and one write per population per update. It is not a hardware performance-counter measurement.


## Distributed MPI build

```bash
cmake -S . -B build-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_ENABLE_OPENMP=ON \
  -DCFD_ENABLE_MPI=ON
cmake --build build-mpi --parallel
ctest --test-dir build-mpi --output-on-failure
mpirun -np 4 ./build-mpi/cfd-distributed --nx 192 --ny 128 --nz 96 --steps 100
```

The Phase-3B distributed path uses MPI between bricks and OpenMP inside each rank, with persistent selective crossing-population halos instead of full-cell traffic. Optional MPI+SYCL execution keeps the local lattice device-resident and supports pinned-host staging by default or explicit GPU-aware MPI direct buffers when the platform supports them. `docs/DISTRIBUTED_RUNTIME.md` describes the 26-neighbour topology, compact halo plan, overlap strategy, device assignment, scaling tools, and restart format.

## Portable GPU build

Current AdaptiveCpp documentation recommends its CMake integration and supports a generic compilation flow for multiple device families. A typical build is:

```bash
cmake -S . -B build-sycl \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_ENABLE_SYCL=ON \
  -DAdaptiveCpp_DIR=/path/to/adaptivecpp/lib/cmake/AdaptiveCpp \
  -DACPP_TARGETS=generic
cmake --build build-sycl --parallel
ctest --test-dir build-sycl --output-on-failure
./build-sycl/cfd-solve lbm-d3q19-sycl
```

A normal GCC/Clang CPU build has no SYCL dependency.

## Performance principles

1. Keep hot fields device-resident across time steps.
2. Prefer SoA/component-major storage for coalesced and SIMD-friendly access.
3. Fuse bandwidth-bound operations only when numerical parity remains explicit.
4. Use compile-time descriptors and simple data-oriented kernels instead of deep hot-path object graphs.
5. Keep a slow/reference implementation for every aggressive optimization where practical.
6. Make reduced-precision storage opt-in and gate it on conservation/solution error.
7. Keep `fast-math` opt-in.
8. Treat measured CPU/GPU performance as hardware-specific; memory reduction and numerical equivalence are the portable guarantees.

## Reference-project relationship

| Reference | What cfd_solvers learns from it | Import policy |
|---|---|---|
| OpenFOAM 14 | FVM architecture, pressure/velocity coupling, multiphase/reacting/CHT capability map | clean-room implementation; no GPL source copied |
| FluidX3D | bandwidth-first LBM, published in-place streaming, compact population storage concepts, heterogeneous/multi-GPU design | concepts/papers only; no FluidX3D source copied or translated |
| Elmer FEM | multiphysics FEM, sparse systems, adaptive methods, coupled physics | clean-room implementation; no GPLv2 source copied |
| openEMS FDTD | Yee/FDTD EM, dispersive materials, ports, PML/NF2FF capability map | clean-room implementation; no GPLv3 source copied |
| Optiland | optical-system/ray abstractions, analysis and backend separation | MIT-compatible reference; independent C++ implementation preferred |
| Geant4 | particle-through-matter track/step/process architecture, geometry navigation and scoring | clean-room architecture reference only; no Geant4 source copied |
| SU2 | PDE-constrained optimization, solver/iteration organization, adjoint/finite-difference workflows | clean-room workflow/numerics reference only; no LGPL source copied |
| csauto | DOE campaign automation, solver-adapter boundary, registry/status/dashboard concepts | clean-room workflow reference only; no GPLv3 source copied |

See `docs/PHASE4A_FVM.md`, `docs/PHASE4B_INCOMPRESSIBLE.md`, `docs/PHASE4C_COLLOCATED.md`, `docs/DISTRIBUTED_RUNTIME.md`, `docs/PHASE3B_RESULTS.md`, `docs/LBM_BOUNDARIES.md`, `docs/SOURCE_MAP.md`, `docs/PERFORMANCE.md` and `THIRD_PARTY.md`.

## Repository layout

```text
include/cfd/core/       data layouts, CPU parallelism and shared linear solvers
include/cfd/fvm/        polyhedral mesh and reusable finite-volume operators
include/cfd/distributed/ topology, halo exchange and rank/device runtime
include/cfd/solvers/    solver APIs grouped by numerical method
src/                    implementations
src/sycl/               accelerator-specific C++20/SYCL implementations
apps/                   solver and benchmark CLIs
scripts/                repeatable benchmark helpers
tests/                  numerical regression tests
docs/                   architecture, provenance, performance and roadmap
```

## License

New cfd_solvers code is MIT licensed. That depends on keeping the implementation clean-room. Do not paste or mechanically translate GPL or FluidX3D-restricted source into this repository.
