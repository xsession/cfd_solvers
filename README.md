# cfd_solvers

A clean-room, C++20-first collection of numerical solvers with one shared high-performance runtime for multicore CPUs and heterogeneous GPUs.

The long-term goal is a readable solver collection covering CFD, FEM multiphysics, FDTD electromagnetics and optics without forcing those different numerical methods into one algorithm. They share memory layouts, execution backends, decomposition, linear algebra, geometry, materials, I/O and diagnostics.

The project is inspired by the capabilities and engineering lessons of OpenFOAM, FluidX3D, Elmer FEM, openEMS FDTD and Optiland. It is **not** a source-code merge or mechanical translation. The source projects have different licenses, and FluidX3D has additional restrictions, so performance techniques are independently implemented from publications and public descriptions.

## v0.6.1 status - nonlinear/adaptive FEM and additional Elmer-class physics

The machine-checkable integration program now reports **137/398 capabilities (34.4%)** complete. The tracker remains conservative: a capability is checked only when implementation and validation are both present. Phase 0 traceability/release infrastructure is complete; the Elmer-class FEM family is now 27/46 (58.7%), openEMS-class FDTD is 8/26, and Optiland-class optics is 12/43.

This checkpoint builds on v0.6.0 and adds:

- generic Newton nonlinear solves with backtracking and ILU(0)-GMRES linearization;
- nonlinear Tri3 Poisson/reaction validation;
- residual/jump error estimation, Dorfler marking and conforming longest-edge adaptive refinement;
- saturated Darcy porous flow;
- 2-D magnetostatics using an out-of-plane vector potential;
- generalized sparse FEM eigenanalysis and a fixed-free bar modal benchmark.

The v0.6.0 breadth work remains part of this tree, including:

- reference-element topology, shape functions, gradients and Gaussian quadrature for Line2/Tri3/Quad4/Tet4/Hex8/Prism6/Pyramid5;
- isoparametric mapping/Jacobians and physical shape gradients;
- reusable Tri3 CSR and matrix-free Laplace assembly;
- mixed Dirichlet/Neumann/Robin scalar FEM boundaries;
- 2-D Poisson/heat/linear-elasticity, electrostatics and DC conduction;
- axisymmetric linear elasticity with the cylindrical `2*pi*r` weak form and hoop strain;
- convergent 3-D Tet4 Poisson assembly;
- 3-D Cartesian Maxwell FDTD baseline;
- heterogeneous dielectric/conductive 1-D FDTD coefficients, hard/soft sources, first-order Mur absorption, time probes and DFT monitors;
- sequential real-ray tracing plus paraxial first-order optics;
- Sellmeier glass dispersion, Jones/Stokes polarization, Fresnel interfaces and normal-incidence multilayer thin films;
- benchmark NDJSON output plus history summarization;
- tag-driven release artifacts: source archives, delta patch, Git bundle and SHA-256 manifest.

The v0.5.0 sparse/FVM/electrochemistry work remains part of this tree, including ILU(0)-GMRES collocated momentum solves, Nernst-Planck/PNP transport, Scharfetter-Gummel fluxes, Butler-Volmer/Faradaic boundaries and galvanic mixed-potential calculations.

See `docs/INTEGRATION_TRACKER.md`, `docs/PHASE4E_FEM_CORE.md`, `docs/PHASE4F_FEM_ADVANCED.md`, `docs/PHASE5_FDTD_BASELINE.md`, `docs/PHASE6_OPTICS_BASELINE.md`, and `docs/VALIDATION.md`.

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
