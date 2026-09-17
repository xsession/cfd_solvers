# cfd_solvers v0.8.0

This release freezes the first broad RF/antenna and SPICE-class circuit checkpoint while retaining the previously validated CFD, LBM, FEM, FDTD, optics, chemistry/electrochemistry and multiphysics work.

## Completion status

The authoritative integration tracker reports **429/576 capabilities (74.5%)** validated. A checked capability requires implementation plus a deterministic regression or physical/manufactured validation path; narrow baselines remain distinct from production-complete features.

- Common HPC runtime: 30/33
- FluidX3D-class LBM: 29/46
- OpenFOAM-class FVM: 83/106
- Elmer-class FEM: 40/49
- openEMS-class FDTD: 22/26
- Optiland-class optics: 40/44
- Chemistry/electrochemistry/corrosion: 57/72
- Coupled multiphysics: 12/15
- Interoperability/workflow: 12/14
- RF/antenna/microwave: 30/62
- SPICE/circuit/compact models: 61/96

## RF and microwave

The RF layer includes generic square complex N-port matrices, S/Z/ABCD conversion, Touchstone SnP v1/v2 handling, real-reference renormalization, mixed-mode 4-port conversion, reference-plane shifting, passivity/reciprocity checks, return-loss/VSWR, K/mu stability, stability circles, transducer gain and load-pull sampling.

Analytical/reference components include TEM lines, microstrip, stripline, CPW, coax and rectangular-waveguide TE10. Antenna utilities include Hertzian/sinusoidal dipoles, phased arrays and polarization metrics. A clean-room center-fed thin-wire MoM baseline provides current distribution and feed impedance, and a PEEC filament model provides resistance/partial-inductance matrices, skin-effect resistance, driven-current solves and SPICE coupled-inductor export.

## SPICE and compact-model baseline

The circuit engine uses modified nodal analysis and supports R/C/L, independent sources, VCVS/VCCS/CCVS/CCCS, switches, mutual inductance, diode, MOS Level-1, BJT, JFET, generic static residual/Jacobian devices and sampled RF N-port elements.

Analyses include DC operating point with homotopy/PN limiting, small-signal AC, backward-Euler/trapezoidal/BDF2 transient, DC/AC/parameter/temperature sweeps, resistor thermal noise integration, numerical sensitivity, Monte Carlo and Fourier/THD. Circuit ports can be directly converted to S-parameters.

Netlist support includes engineering suffixes, `.MODEL`, `.PARAM`, expressions, recursive hierarchical `.SUBCKT` flattening, instance/default/local parameters, scoped local models, `.FUNC`, `.INCLUDE` and `.LIB` section selection.

The OSDI/OpenVAF integration is deliberately split: v0.8.0 provides safe dynamic-library discovery, ABI version/count metadata and descriptor-pointer access. Full OSDI descriptor evaluation/state/charge/noise integration remains open rather than being claimed prematurely.

## Validation

- GCC 14 + OpenMP: **43/43 CTest targets pass**.
- Focused Clang 17 and serial GCC RF/SPICE paths pass.
- Focused ASan+UBSan + leak detection passes SPICE hierarchy/file parsing, OSDI loader, PEEC/MoM and nonlinear circuit smoke paths.
- `git diff --check` is required clean before archive generation.

## Research/licensing

Primary RF references are Apache-2.0 Palace, MIT OpenSEMBA FDTD and MIT OpenNEC, plus published/behavioral references for PEEC and microwave network workflows. Primary circuit behavior reference is modified-BSD ngspice; Xyce/QucsatorRF are clean-room capability references. OpenVAF/OSDI is the planned external compact-model route. Exact pins and license policy are recorded in `docs/upstreams.json`, `THIRD_PARTY.md` and `docs/RF_CIRCUIT_RESEARCH.md`.

## Deliberate remaining gaps

High-value remaining gaps include production thin-wire MoM/NEC accuracy and geometry, 3-D RF FEM wave ports/eigenmodes, vector fitting/model-order reduction, PEEC capacitance/proximity effects, sparse MNA, adaptive LTE transient stepping, pole-zero/PSS/harmonic balance, full OSDI descriptor evaluation, broader vendor-model compatibility, RF result-file interoperability and hardware-backed GPU/distributed RF.
