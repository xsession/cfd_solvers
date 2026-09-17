# RF, antenna and SPICE-class solver research

Research refresh: 2026-09-16

This document extends the unified solver program with two additional families:

- **Phase 10 — RF/antenna and microwave electromagnetic analysis**
- **Phase 11 — SPICE-class circuit and RF network simulation**

The project remains clean-room. Copyleft projects are behavioral/algorithmic references unless an explicitly separate adapter is used. No source code from GPL projects is copied or mechanically translated into the MIT core.

## RF / antenna solver landscape

| Project | License | Numerical focus | High-value capabilities for this project | Integration policy |
|---|---|---|---|---|
| Palace | Apache-2.0 | 3-D high-order FEM for Maxwell | electrostatic/magnetostatic extraction, driven frequency domain, eigenmodes/Q, transient EM, lumped/wave ports, adaptive frequency sweeps, AMR, CPU/GPU/MPI | permissive reference; independent lightweight C++ implementations may reuse published equations/architecture concepts |
| OpenSEMBA FDTD | MIT | EMC-oriented 3-D FDTD | graded mesh, CPML/Mur, anisotropic/dispersive media, wires/multiwires, MTL embedded in FDTD, SPICE coupling, plane waves, Hertzian dipoles, Huygens surfaces, thin slots/composites, transfer probes, NF2FF | strongest permissive reference for future EMC/RF coupling |
| OpenNEC | MIT | NEC-2 thin-wire MoM | center-fed wires, antenna currents, impedance, far-field/radiation patterns, wire networks | preferred permissive thin-wire antenna reference |
| openEMS | GPL-3.0 | EC-FDTD | Cartesian/cylindrical grids, PML, TFSF, lumped RLC, multipole Drude/Lorentz/Debye, NF2FF, SAR, ports, MPI/SIMD | clean-room capability/reference only |
| Meep | GPL-2.0 | FDTD | 1/2/3-D/cylindrical, MPI, CW frequency-domain, eigenmodes, HDF5, mode decomposition/S-parameters, near-to-far, LDOS, force/stress, adjoint inverse design | clean-room capability/reference only |
| scuff-em | GPL-2.0/3.0 | surface-current BEM | RF antennas, coax/MRI/passive RF structures, scattering, electrostatics, capacitance | clean-room BEM/RF reference only |

Pinned permissive references used in the tracker:

- Palace: `awslabs/palace` commit `dff1ef9d1794ff8733ac8b939fb185b181786c1d`.
- OpenSEMBA FDTD: `OpenSEMBA/fdtd` commit `82e310c9ac43d729256b23a42d8c20799ad42b8a`.
- OpenNEC: `maurymarkowitz/OpenNEC` commit `acafb3796066562cf3331fc9106485f33b6dbcc6`.

### RF features worth integrating

The RF phase should not be a second, disconnected EM stack. Existing FDTD/FEM/material/IO infrastructure should be reused.

1. **Network layer**
   - complex Z/Y/S/ABCD conversion;
   - cascading and de-embedding;
   - Touchstone S1P/S2P/SnP import/export;
   - renormalization to arbitrary real/complex port impedances;
   - passivity/reciprocity/stability metrics;
   - impedance/admittance matching utilities;
   - frequency sweeps.
2. **Transmission structures**
   - ideal/lossy TEM transmission line;
   - microstrip/stripline/coplanar quasi-static models;
   - coax/waveguide analytic modes;
   - multiconductor transmission-line matrices;
   - discontinuity/equivalent circuit models.
3. **Antenna analysis**
   - sinusoidal thin-wire dipole validation baseline;
   - NEC-like thin-wire MoM;
   - wire junctions/loads/ground models;
   - input impedance, VSWR, return loss;
   - current distribution;
   - 2-D/3-D radiation pattern, gain/directivity, efficiency;
   - polarization/axial ratio;
   - array factors and phased arrays;
   - mutual impedance/coupling;
   - near-field and NF2FF bridge from FDTD.
4. **Full-wave frequency-domain FEM**
   - Nedelec/edge elements;
   - curl-curl complex system;
   - PEC/PMC/impedance/radiation boundaries;
   - driven current/voltage excitation;
   - lumped and wave ports;
   - eigenmode/resonance/Q;
   - S-parameter extraction;
   - adaptive frequency sweep/model-order reduction;
   - AMR error indicators.
5. **EMC/RF co-simulation**
   - thin wires/slots/sheets;
   - cable and multiconductor TL models embedded in field solver;
   - lumped circuit coupling;
   - plane-wave/Hertzian/Huygens excitations;
   - transfer functions and shielding effectiveness;
   - SAR and field exposure;
   - shared geometry/material database.

## SPICE / circuit solver landscape

| Project | License | Focus | High-value capabilities | Integration policy |
|---|---|---|---|---|
| ngspice 47 | modified BSD (project FAQ) | SPICE3/XSPICE/CIDER mixed-level simulator | DCOP, DC sweeps, AC, transient, noise, pole-zero, distortion, sensitivity, S-parameters, KLU, XSPICE, OSDI/OpenVAF Verilog-A, broad compact-model library | primary permissive behavioral/numerical reference |
| Xyce 7.10 | GPL-3.0 | scalable DAE-based SPICE | DC/TRAN/AC/noise, harmonic balance, multi-time PDE, model-order reduction, MPI, modern compact models, ADMS Verilog-A | clean-room architecture/validation reference only |
| QucsatorRF 1.0.7 development line | GPL-2.0 | microwave/RF circuit analysis | S-parameters, microstrip/waveguide components, harmonic-balance-related RF workflows, matching/stability/noise formulas | clean-room RF network/device reference only |
| Qucs-S | GPL-2.0 | GUI/front-end to ngspice/Xyce/QucsatorRF | backend abstraction, SPICE model import, sweeps, RF workflow, Touchstone viewer | UX/workflow reference only |
| GNUCAP | GPL family | extensible circuit simulator | plugin/model architecture, SPICE mode, mixed-signal direction | clean-room architecture reference only |

Current official releases found during this refresh:

- ngspice **47**, released 2026-08-11.
- Xyce **7.10**, released 2025-08.
- QucsatorRF development NEWS lists **1.0.7**.

### SPICE features worth integrating

1. **Circuit kernel**
   - sparse modified nodal analysis (MNA);
   - DAE residual/Jacobian abstraction;
   - Newton with voltage limiting, source stepping and gmin stepping;
   - direct/Krylov sparse linear solvers;
   - adaptive timestep with LTE control;
   - Gear/BDF and trapezoidal integration.
2. **Analyses**
   - operating point;
   - DC sweep;
   - AC small signal;
   - transient;
   - noise;
   - pole-zero;
   - sensitivity;
   - distortion;
   - S-parameter/network analysis;
   - FFT/Fourier;
   - periodic steady state / harmonic balance;
   - parameter/temperature/Monte-Carlo sweeps.
3. **Primitive and controlled devices**
   - R, L, C;
   - independent voltage/current sources;
   - VCVS/VCCS/CCVS/CCCS;
   - switches;
   - ideal/real transmission lines;
   - mutual inductors/transformers;
   - behavioral sources.
4. **Semiconductor compact models**
   - diode;
   - BJT;
   - JFET;
   - MOS level 1/2/3;
   - BSIM3/4/CMG-class model adapters;
   - VBIC/HICUM/MEXTRAM;
   - VDMOS/power devices;
   - MESFET/HEMT;
   - thermal/self-heating terminals.
5. **Model interchange**
   - `.MODEL`, `.SUBCKT`, `.PARAM`, `.FUNC`, `.LIB`, `.INCLUDE`;
   - SPICE suffix/expression parser;
   - vendor model compatibility modes;
   - Verilog-A parser/compiler interface;
   - OpenVAF/OSDI dynamic model adapter;
   - optional external ngspice/Xyce adapters for cross-validation.
6. **RF circuit co-simulation**
   - S/Y/Z/ABCD N-port blocks;
   - Touchstone devices;
   - microstrip/waveguide components;
   - harmonic balance;
   - field-solver port extraction -> circuit N-port;
   - circuit loads -> FDTD/FEM lumped ports;
   - EM/circuit iterative and monolithic coupling.

## Implemented in the first RF/SPICE checkpoint

The current local continuation implements the following clean-room baselines:

### RF network and antenna

- two-port complex Z <-> S conversion;
- S <-> ABCD conversion and two-port cascade;
- ideal/lossy transmission-line ABCD/S model;
- Hammerstad/Jensen-style quasi-static microstrip impedance/effective-permittivity model;
- Touchstone v1 S2P RI/MA/DB reader and RI writer;
- sinusoidal thin-wire dipole far-field solver;
- numerical radiated-power/radiation-resistance integration;
- dipole directivity and normalized radiation pattern.

Validation includes the half-wave dipole values `Rrad ~= 73.13 ohm` and `D ~= 1.64`, exact S/Z round-trip, matched lossless-line `|S21|=1`, and Touchstone parsing.

### SPICE-class circuit kernel

- clean-room dense-MNA correctness baseline;
- resistor, capacitor and inductor stamping;
- independent voltage/current sources;
- nonlinear diode model with Newton linearization;
- MOSFET Level-1/Shichman-Hodges baseline with channel-length modulation;
- DC operating point;
- complex small-signal AC analysis;
- backward-Euler transient analysis;
- SPICE engineering suffix parser;
- basic netlist parser for R/C/L/V/I/D/M and `.MODEL D/NMOS/PMOS`.

Validation includes an exact resistive divider, RC `-3 dB / -45 degree` corner, nonlinear diode bias point, Level-1 MOS operating point, netlist/model parsing, and RC charging transient.

## Next RF/SPICE implementation order

1. move MNA storage from dense correctness baseline to shared sparse CSR/KLU-style abstraction;
2. source/gmin stepping, junction voltage limiting and adaptive transient LTE;
3. controlled sources, switches, mutual inductance and transmission-line elements;
4. DC/parameter sweeps, noise, pole-zero and sensitivity;
5. BJT/JFET/MOS2/3 compact models and `.SUBCKT/.PARAM/.FUNC` parser;
6. generic N-port/Touchstone component integrated into MNA;
7. NEC-like thin-wire MoM current/impedance solver;
8. RF FDTD/FEM ports -> N-port extraction and circuit co-simulation;
9. edge-element driven/eigenmode FEM RF solver;
10. harmonic balance and Verilog-A/OSDI adapter layer.


## 2026-09-16 continuation: hierarchy, OSDI seam and thin-wire MoM

The circuit front end now resolves `.FUNC`, recursive `.INCLUDE`, `.LIB <file> <section>`, nested `.SUBCKT` parameter passing, and per-instance local `.MODEL` namespacing before MNA stamping. File paths are resolved relative to the containing netlist. The string-only parser rejects file directives intentionally so path semantics remain explicit.

The OSDI integration has a first dynamic-library discovery seam. It loads the four standard exported entry symbols (`OSDI_VERSION_MAJOR`, `OSDI_VERSION_MINOR`, `OSDI_DESCRIPTORS`, `OSDI_NUM_DESCRIPTORS`) and validates descriptor presence without depending on the still-evolving `OsdiDescriptor` layout. This is sufficient to validate OpenVAF-produced library discovery, but descriptor instantiation/evaluation/stamping is still tracked separately.

The RF antenna path now includes a readable thin-wire EFIE Method-of-Moments reference solver. Its present regression checks matrix-driven finite/passive input impedance and the mirror symmetry of a center-fed straight wire. It is not yet claimed to reproduce NEC-grade feed reactance/resistance across canonical benchmarks; that remains a separate validation item.
