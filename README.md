# cfd_solvers

A clean-room, C++20-first collection of numerical solvers with one shared high-performance runtime for multicore CPUs and heterogeneous GPUs.

The long-term goal is a readable solver collection covering CFD, FEM multiphysics, FDTD electromagnetics and optics without forcing those different numerical methods into one algorithm. They share memory layouts, execution backends, decomposition, linear algebra, geometry, materials, I/O and diagnostics.

The project is inspired by the capabilities and engineering lessons of OpenFOAM, FluidX3D, Elmer FEM, openEMS FDTD and Optiland. It is **not** a source-code merge or mechanical translation. The source projects have different licenses, and FluidX3D has additional restrictions, so performance techniques are independently implemented from publications and public descriptions.

## v0.19.3 status - DFN/P2D and Docker deployment

v0.19.3 completes Phase 16A with the full Doyle-Fuller-Newman / P2D distributed
porous-electrode model (`cfd/battery/dfn.hpp`): per-node spherical solid particles,
distributed electrolyte concentration, and distributed solid/electrolyte potentials.
It also adds a Docker-based deployment system (multi-stage `Dockerfile`,
`docker-compose.yml`, and a no-bash Windows helper). See `docs/RELEASE_0_19_3.md`
for usage and validation.

Validation: `cfd-v0193-dfn-tests` plus the `battery-dfn` smoke case pass; the DFN
regression, the DFN source and the updated case dispatcher compile clean under
`-Wall -Wextra -Wshadow` (C++20). The full CTest suite is exercised via the
`cfd-test` compose profile. Phase 16A is **18/18** and the tracker is
**718/825 = 87.0%**.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFD_ENABLE_OPENMP=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/cfd-solve --list
```

Practical examples and repeatable timing baselines:

```bash
python3 examples/run_benchmarks.py --build-dir build --repeat 3 --output examples/results/my_machine.json
python3 examples/compare_benchmarks.py examples/results/v0.19.0_reference.json examples/results/my_machine.json
```

Use `--quick` for CI/smoke runs. `simulation_ms` excludes case construction/setup so it can be compared separately from `setup_ms`; compare only runs with equivalent case mode/scale and comparable compiler/backend/power settings. The runner records build flags/compiler/CPU affinity, and the comparison tool flags checksum drift. See `examples/README.md` and `examples/results/v0.19.0_reference_summary.md`.

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
./build/cfd-solve fvm-euler-weno1d
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
./build/cfd-solve battery-dfn
./scripts/integration_status.py
```

## Docker deployment

A multi-stage `Dockerfile` builds the solvers (Release, OpenMP, optional MPI/HDF5)
and ships the binaries plus the source tree in a slim non-root runtime image.
Build args: `CMAKE_BUILD_TYPE`, `CFD_ENABLE_OPENMP`, `CFD_ENABLE_MPI`,
`CFD_ENABLE_NATIVE_ARCH`, `CFD_ENABLE_HDF5` (all default to a CPU/OpenMP build).

```bash
# Build the runtime image
docker build -f Dockerfile -t cfd_solvers:dev .

# Run a case
docker run --rm -e OMP_NUM_THREADS=8 --entrypoint /cfd_solvers/build/cfd-solve \
  cfd_solvers:dev lbm-d3q19-cpu
docker run --rm --entrypoint /cfd_solvers/build/cfd-solve cfd_solvers:dev battery-dfn

# docker compose: run / full CTest suite / MPI distributed case
docker compose run cfd lbm-d3q19-cpu
docker compose run --profile test cfd-test
docker compose run --profile mpi cfd-mpi          # requires the MPI build
```

On Windows use the no-bash helper (`cfd-docker.bat` delegates to `cfd-docker.ps1`):

```bat
cfd-docker.bat build
cfd-docker.bat run --threads 8 battery-dfn
cfd-docker.bat test
cfd-docker.bat mpi --ranks 4
```

`CFD_TAG` and `CFD_MPI=1` override the image tag and enable the MPI build.

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
examples/               practical solver-family examples and comparable timing baselines
scripts/                repeatable benchmark helpers
tests/                  numerical regression tests
docs/                   architecture, provenance, performance and roadmap
```

## License

New cfd_solvers code is MIT licensed. That depends on keeping the implementation clean-room. Do not paste or mechanically translate GPL or FluidX3D-restricted source into this repository.
