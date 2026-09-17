# cfd_solvers

A clean-room, C++20-first collection of numerical solvers with one shared high-performance runtime for multicore CPUs and heterogeneous GPUs.

The long-term goal is a readable solver collection covering CFD, FEM multiphysics, FDTD electromagnetics and optics without forcing those different numerical methods into one algorithm. They share memory layouts, execution backends, decomposition, linear algebra, geometry, materials, I/O and diagnostics.

The project is inspired by the capabilities and engineering lessons of OpenFOAM, FluidX3D, Elmer FEM, openEMS FDTD and Optiland. It is **not** a source-code merge or mechanical translation. The source projects have different licenses, and FluidX3D has additional restrictions, so performance techniques are independently implemented from publications and public descriptions.

## v0.9.0 status - CST-class EM, particle and bioelectromagnetic expansion

The machine-checkable integration program now reports **473/634 capabilities (74.6%)** complete. Phase 10 RF/antenna/microwave is **46/70 (65.7%)**, Phase 11 SPICE/circuit/compact-model simulation is **80/98 (81.6%)**, and the new Phase 12 CST-class expansion is **8/47 (17.0%)**. The CPU validation matrix now contains **52 CTest targets**.

This checkpoint adds or consolidates:

- a driven complex 1-D frequency-domain Maxwell/Helmholtz solver with PEC boundaries, conductive loss and staggered magnetic-field reconstruction;
- dielectric-loaded 1-D PEC cavity electromagnetic eigenmodes with analytical frequency validation;
- reusable lowest-order first-kind Nedelec Tri3 basis functions and a driven 2-D curl-conforming edge-element Maxwell FEM baseline with global edge orientation and PEC tangential constraints;
- non-relativistic 3-D charged-particle tracking through a Boris pusher plus a periodic 1-D electrostatic PIC baseline with cloud-in-cell deposition and spectral Poisson fields;
- RMS electric-field to SAR conversion and an implicit 2-D Pennes bioheat solver with conduction, perfusion, metabolic heat and spatial SAR loading;
- CLI smoke cases for the new EM, PIC and bioheat families;
- `docs/CST_COVERAGE_AUDIT.md`, which maps the actual implementation against the wider EM/thermal/structural/particle/plasma/bioelectromagnetic domains and explicitly identifies remaining gaps;
- tracker reconciliation for RF/SPICE capabilities that were already implemented and tested but still marked incomplete, including complex power-wave renormalization, passivity/causality tools, PEEC capacitance, connected-wire MoM, sparse MNA, DAE devices, B-sources, RAW I/O, noise/distortion and self-heating.

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
