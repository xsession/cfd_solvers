# cfd_solvers v0.19.0

## Phase 16A battery baseline

v0.19.0 starts lithium-ion battery cell physics with:

- conservative finite-volume spherical solid diffusion;
- SPM and SPMe reduced-order models;
- graphite/NMC OCV helpers and symmetric Butler-Volmer kinetics;
- through-cell electrolyte concentration and potential-loss baseline;
- solid-electrode/contact ohmic losses;
- lumped and 1-D through-cell thermal models;
- SEI, plating, active-material-loss and particle-cracking state baselines;
- battery EIS and piecewise-current drive-cycle workflows.

Phase 16A is **15/18** in this release. Full DFN/P2D, 3-D battery thermal coupling and multi-cell balancing remain open.

## Practical examples and repeatable simulation timing

The release also introduces a dedicated `examples/` tree with practical cases for 15 large feature families:

1. HPC multigrid Poisson;
2. LBM Taylor-Green vortex;
3. FVM lid-driven cavity;
4. FEM transient bar vibration;
5. FDTD dielectric pulse + CPML;
6. geometric optics lens ray bundle;
7. reactive chemistry + corrosion;
8. Joule electrothermal multiphysics;
9. RF/transmission-line frequency sweep;
10. SPICE BDF2 RC transient;
11. electromagnetic PIC beam;
12. Hertz-Mindlin DEM settling;
13. photoacoustic k-space propagation;
14. PN-junction TCAD;
15. SPMe battery drive cycle.

Phase 9 is workflow/orchestration rather than a numerical kernel; the benchmark runner itself exercises that layer.

Every example writes one `CFD_BENCH` JSON record with setup time, simulation-only time, total time, problem size, integration steps, workload units, throughput and a checksum. `run_benchmarks.py` stores repetitions plus host/build metadata. `compare_benchmarks.py` compares median simulation time and flags workload/checksum changes so a timing comparison cannot silently compare different numerical cases.

A quick all-example run is part of CTest. The repository also stores `examples/results/v0.19.0_reference.json` as the first three-repeat Release-mode comparison baseline for this environment, with a human-readable summary in `examples/results/v0.19.0_reference_summary.md`.

## Validation

- focused v0.19.0 battery regression: pass;
- all-example quick timing smoke: pass;
- strict new-source/example compile: pass;
- full default HDF5/Python/examples matrix: **163/163 CTest targets passed**.

## v0.19.0 reference timing on the release host

The repository reference was measured in Release mode with native-architecture optimization and OpenMP enabled, three repeats per case. The JSON file records CPU affinity, compiler/version, build/backend flags and workload checksums. These numbers are useful as the first regression anchor, not as cross-machine performance claims.

| Example | Median simulation time | Median throughput |
|---|---:|---:|
| `cfd-example-phase01-hpc` | 2.1312 ms | 4.54e+07 unknowns_cycles/s |
| `cfd-example-phase02-lbm` | 570.8918 ms | 1.55e+07 cell_steps/s |
| `cfd-example-phase03-fvm` | 10.1777 ms | 4.53e+05 cell_steps/s |
| `cfd-example-phase04-fem` | 10.2966 ms | 3e+07 dof_steps/s |
| `cfd-example-phase05-fdtd` | 37.3542 ms | 2.74e+08 cell_steps/s |
| `cfd-example-phase06-optics` | 3.0399 ms | 1.03e+07 rays/s |
| `cfd-example-phase07-chemistry-corrosion` | 4.4384 ms | 2.31e+08 cell_steps/s |
| `cfd-example-phase08-multiphysics` | 32.8977 ms | 5.6e+06 element_steps/s |
| `cfd-example-phase10-rf` | 8.6088 ms | 5.81e+06 frequency_points/s |
| `cfd-example-phase11-spice` | 7.5633 ms | 2.64e+06 time_steps/s |
| `cfd-example-phase12-em-pic` | 496.2486 ms | 6.45e+05 particle_steps/s |
| `cfd-example-phase13-dem` | 19.1342 ms | 5.64e+06 particle_steps/s |
| `cfd-example-phase14-acoustics` | 857.3967 ms | 3.06e+06 cell_steps/s |
| `cfd-example-phase15-tcad` | 37.3258 ms | 4.05e+04 nodes_bias_points/s |
| `cfd-example-phase16a-battery` | 8.2000 ms | 2.63e+07 state_steps/s |

## Tracker

v0.19.0 expands the machine-counted tracker with Phase 16A. After the implemented baseline, progress is **715/825 = 86.7%**, with Phase 16A at **15/18 = 83.3%**.
