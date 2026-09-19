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
| `cfd-example-phase03-vof` | bounded VOF column advection with CSF diagnostics |
| `cfd-example-phase03-thin-film` | pressure-coupled capillary thin-film relaxation |
| `cfd-example-phase03-euler-euler` | coupled two-fluid volume fraction and drag transport |
| `cfd-example-phase03-compressible-vof` | barotropic compressible VOF mass/fraction transport |
| `cfd-example-phase03-spray` | coupled spray injection, drag, evaporation, and carrier sources |
| `cfd-example-phase03-particle-rheology` | many-particle collisions, breakup, dense drag, and granular stress |
| `cfd-example-phase03-reacting-hooks` | premixed/non-premixed flame and soot/radiation source hooks |
| `cfd-example-phase04-fem` | transient vibrating cantilever bar |
| `cfd-example-phase05-fdtd` | dielectric pulse with CPML |
| `cfd-example-phase06-optics` | biconvex-lens ray bundle |
| `cfd-example-phase07-chemistry-corrosion` | implicit reactor + phase-field corrosion |
| `cfd-example-phase08-multiphysics` | DC conduction -> Joule heat -> transient thermal FEM |
| `cfd-example-phase10-rf` | microstrip/transmission-line frequency sweep |
| `cfd-example-phase11-spice` | BDF2 RC transient |
| `cfd-example-phase12-em-pic` | periodic electromagnetic PIC beam |
| `cfd-example-phase12-surface-mom` | small PEC plate with RWG surface-current MoM |
| `cfd-example-phase13-dem` | granular settling with Hertz-Mindlin DEM |
| `cfd-example-phase14-acoustics` | 2-D photoacoustic k-space propagation |
| `cfd-example-phase15-tcad` | PN-junction equilibrium and I-V sweep |
| `cfd-example-phase16a-battery` | SPMe current/rest/charge drive cycle |
| `cfd-example-battery-dfn` | distributed DFN/P2D discharge profile |
| `cfd-example-battery-pack` | series-pack balancing profile |
| `cfd-example-battery-thermal` | electrothermal pack and 3-D voxel heat balance |
| `cfd-example-real-openfoam-channel` | ventilated channel / wind-tunnel FVM case with minimal OpenFOAM case and ParaView VTK export |
| `cfd-example-real-fdtd-radome-vtk` | dielectric-window/radome FDTD pulse with ParaView VTK snapshots |
| `cfd-example-real-lbm3d-model` | 3-D voxelized Ahmed-style bluff body with live JSON frames, VTK snapshots, and Q-criterion output |

Phase 9 is workflow/orchestration rather than a simulation kernel, so the benchmark runner itself is the practical example for that layer.

## Guided walkthroughs

Every example accepts `--quick` for a small deterministic case and `--scale N` for a larger
case. The finite-volume and battery walkthroughs also accept `--output path.csv`; they write a CSV field
profile while keeping the `CFD_BENCH` record on standard output for automated timing.

Start with the finite-volume examples:

```bash
build/cfd-example-phase03-fvm --quick
build/cfd-example-phase03-vof --quick --output build/vof_final.csv
build/cfd-example-phase03-thin-film --quick --output build/thin_film_final.csv
build/cfd-example-phase03-euler-euler --quick --output build/euler_euler_final.csv
build/cfd-example-phase03-compressible-vof --quick --output build/compressible_vof_final.csv
build/cfd-example-phase03-spray --quick --output build/spray_final.csv
build/cfd-example-phase03-particle-rheology --quick --output build/particle_rheology_final.csv
build/cfd-example-phase03-reacting-hooks --quick --output build/reacting_hooks_final.csv
```

The cavity example reports the final continuity residual. The VOF example transports a liquid
column through a closed channel, reports the initial/final liquid volume and boundedness, and
writes `x_m,y_m,volume_fraction` plus CSF force components. Plot `volume_fraction` in any CSV
viewer to see the interface; the force columns are the capillary diagnostic used by the existing
CSF regression.

The thin-film example couples a lubrication flux to capillary pressure reconstructed from the
height curvature. It uses impermeable boundaries, so the inventory check exposes conservation
while the CSV shows the capillary pressure and relaxed height field. The explicit baseline is
intended for stable small time steps and unit-depth 2-D coating or film cases.

The Euler–Euler walkthrough starts with a dispersed-phase volume-fraction blob and two phase
velocities. It advances conservative internal phase fluxes, applies equal-and-opposite drag,
and writes the volume fraction, phase velocities, and local drag force. Boundary faces are
impermeable in this compact baseline, so the dispersed-phase inventory is a useful conservation
check while the slip velocity relaxes.

The compressible-VOF walkthrough transports liquid fraction and mixture mass with the same
prescribed face flux, then recovers pressure from a barotropic liquid bulk-modulus and gas
power-law EOS. It is an isothermal baseline: energy, phase change, and a pressure-velocity
projection are intentionally outside this increment. Closed internal fluxes make both liquid
volume and total mixture mass straightforward conservation checks.

The spray walkthrough injects parcels through a time-windowed mass-flow schedule, applies Stokes
drag and thermal relaxation against a carrier field, and evaporates parcels with the D² law. It
exports carrier vapor-mass, momentum-reaction, and latent-heat source fields. The parcel model is
intentionally compact: collision statistics, wall films, and a coupled carrier solve remain
separate future layers.

The particle/rheology walkthrough uses a cell-linked broad phase to resolve many simultaneous
collisions, Weber-gated coalescence and breakup, then evaluates hindered dense-phase drag and a
regularized granular stress closure on the same mesh. The carrier and particle momentum sources
are equal and opposite, making the CSV useful for checking coupling signs.

The reacting-hooks walkthrough evaluates both progress-variable flame regimes. Premixed source
strength follows the bounded reaction progress, while the non-premixed source peaks around the
configured stoichiometric mixture fraction. A soot model turns heat release into a transported
soot source and optically-thin radiative sink that can be passed directly to thermal transport.

The battery examples build from the simplest reduced model toward the distributed model:

```bash
build/cfd-example-phase16a-battery --quick
build/cfd-example-battery-pack --quick --output build/pack.csv
build/cfd-example-battery-thermal --quick --output build/thermal.csv
build/cfd-example-battery-dfn --quick --output build/dfn.csv
```

The SPMe example is the fast drive-cycle baseline. The pack profile makes cell-to-cell balancing
visible. The thermal profile exposes supplied heat, boundary cooling and the energy residual. The
DFN profile adds one row per distributed through-cell time step, including surface stoichiometry,
electrolyte bounds, reaction overpotential and ohmic drop. A useful first check is that the DFN
negative surface stoichiometry decreases during discharge while the positive surface
stoichiometry increases.

The optics library also provides clean-room sequential-lens interchange for the common surface
subset shared by Zemax ZMX, CODE V SEQ and OSLO LEN text files. Vendor-specific operands are
left untouched by the core; import/export is intentionally limited to surface type, radius or
curvature, spacing, aperture, refractive index, conic constant and four even-asphere terms.

The surface-MoM walkthrough solves a two-triangle PEC plate with a delta-gap feed. It writes one
row per RWG/half-RWG edge when `--output` is supplied, and reports the feed impedance, residual
and a far-field power checksum. The assembler uses an explicit regularized three-point reference
rule, so it is suitable for learning and small regression cases rather than large production
meshes:

```bash
build/cfd-example-phase12-surface-mom --quick --output build/surface_mom.csv
```

## Visualization examples

These examples write files that can be opened directly in ParaView:

```bash
build/cfd-example-real-openfoam-channel --quick
build/cfd-example-real-fdtd-radome-vtk --quick
build/cfd-example-real-lbm3d-model --quick
```

Default outputs are written under `examples/output/openfoam_channel`, `examples/output/fdtd_radome`, and `examples/output/lbm3d_model`. Use `--output <dir>` to redirect them. The channel example writes `paraview/channel_final.vtk` plus a minimal OpenFOAM case under `openfoam` (`constant/polyMesh`, `constant/transportProperties`, `system/controlDict`, and `0/U`/`0/p`). The FDTD example writes `radome_mid.vtk`, `radome_late.vtk`, and `radome_final.vtk`. The 3-D LBM example writes full-volume `vtk/frame_*.vtk`, `vtk/q_criterion_final.vtk`, and a polling state at `live/state.json`.

### Live 3-D viewer

Start the viewer in one terminal, then run the model example in another:

```bash
python3 examples/live_3d_viewer.py --dir build/example_outputs/cfd-example-real-lbm3d-model/repeat_0/live
build/cfd-example-real-lbm3d-model --quick --output build/example_outputs/cfd-example-real-lbm3d-model/repeat_0
```

Open `http://127.0.0.1:8765/`. The browser polls the latest state, so it works while the simulation is still running. With `--stl path/to/model.stl`, the example loads and normalizes an ASCII STL model before voxelization. The example applies stationary mid-link bounce-back to the voxelized body and a small streamwise body force; the outer domain remains periodic, making this a compact obstacle-flow regression rather than a calibrated Ahmed-body drag benchmark.

## Reproducible timing

Build Release mode with the same compiler flags and backend, then run:

```bash
python3 examples/run_benchmarks.py --build-dir build --output examples/results/my_machine.json --repeat 3
```

Use `--quick` for CI/smoke execution. Use the normal mode for performance comparisons. `--scale N` scales the larger workloads where supported. The runner includes the VOF, thin-film, Euler–Euler, compressible-VOF, spray, and all three battery walkthroughs, so a passing smoke run covers every practical example target built by CMake.

Compare two runs:

```bash
python3 examples/compare_benchmarks.py examples/results/v0.19.0_reference.json examples/results/new.json
```

The comparison uses median `simulation_ms` for repeated cases and reports speedup plus percentage runtime change. It warns when mode/build metadata differs and reports `DRIFT` when the numerical checksum changes beyond the selected tolerance. Meaningful comparisons require the same case mode/scale and similar CPU/GPU power-management settings.

The repository baseline is summarized in `results/v0.19.0_reference_summary.md`; the full machine-readable records are in `results/v0.19.0_reference.json`.

## Stored reference baseline

`examples/results/v0.19.0_reference.json` is the first repository reference run. It records repetitions together with CPU affinity, compiler identity, build type and enabled backends. Treat it as a baseline for regression comparison on similar execution conditions, not as a universal performance number.

The JSON format is described by `examples/benchmark_suite.schema.json`. `compare_benchmarks.py` also checks the case problem size, step count and checksum; use `--fail-on-workload-change` in automated performance gates.
