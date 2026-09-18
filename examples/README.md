# Practical examples and performance baselines

Each executable demonstrates one major simulation family with a recognizable engineering case and emits one machine-readable timing line:

```text
CFD_BENCH {"schema":1,...,"setup_ms":...,"simulation_ms":...,"throughput":...,"checksum":...}
```

`simulation_ms` intentionally excludes object construction, mesh/profile creation, and input setup. The checksum is not a validation tolerance; it is a convenient fingerprint to detect accidental changes to the workload when comparing timings. The suite JSON also records compiler version, build type, major backend flags, CPU affinity and benchmark mode so future comparisons can be audited.

## Covered feature families

| Executable | Practical workload |
|---|---|
| `cfd-example-phase01-hpc` | geometric-multigrid Poisson solve |
| `cfd-example-phase02-lbm` | D3Q19 Taylor-Green vortex |
| `cfd-example-phase03-fvm` | incompressible lid-driven cavity |
| `cfd-example-phase04-fem` | transient vibrating cantilever bar |
| `cfd-example-phase05-fdtd` | dielectric pulse with CPML |
| `cfd-example-phase06-optics` | biconvex-lens ray bundle |
| `cfd-example-phase07-chemistry-corrosion` | implicit reactor + phase-field corrosion |
| `cfd-example-phase08-multiphysics` | DC conduction -> Joule heat -> transient thermal FEM |
| `cfd-example-phase10-rf` | microstrip/transmission-line frequency sweep |
| `cfd-example-phase11-spice` | BDF2 RC transient |
| `cfd-example-phase12-em-pic` | periodic electromagnetic PIC beam |
| `cfd-example-phase13-dem` | granular settling with Hertz-Mindlin DEM |
| `cfd-example-phase14-acoustics` | 2-D photoacoustic k-space propagation |
| `cfd-example-phase15-tcad` | PN-junction equilibrium and I-V sweep |
| `cfd-example-phase16a-battery` | SPMe current/rest/charge drive cycle |

Phase 9 is workflow/orchestration rather than a simulation kernel, so the benchmark runner itself is the practical example for that layer.

## Reproducible timing

Build Release mode with the same compiler flags and backend, then run:

```bash
python3 examples/run_benchmarks.py --build-dir build --output examples/results/my_machine.json --repeat 3
```

Use `--quick` for CI/smoke execution. Use the normal mode for performance comparisons. `--scale N` scales the larger workloads where supported.

Compare two runs:

```bash
python3 examples/compare_benchmarks.py examples/results/v0.19.0_reference.json examples/results/new.json
```

The comparison uses median `simulation_ms` for repeated cases and reports speedup plus percentage runtime change. It warns when mode/build metadata differs and reports `DRIFT` when the numerical checksum changes beyond the selected tolerance. Meaningful comparisons require the same case mode/scale and similar CPU/GPU power-management settings.

The repository baseline is summarized in `results/v0.19.0_reference_summary.md`; the full machine-readable records are in `results/v0.19.0_reference.json`.

## Stored reference baseline

`examples/results/v0.19.0_reference.json` is the first repository reference run. It records repetitions together with CPU affinity, compiler identity, build type and enabled backends. Treat it as a baseline for regression comparison on similar execution conditions, not as a universal performance number.

The JSON format is described by `examples/benchmark_suite.schema.json`. `compare_benchmarks.py` also checks the case problem size, step count and checksum; use `--fail-on-workload-change` in automated performance gates.
