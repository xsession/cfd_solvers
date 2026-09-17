# CST-style coverage audit and expansion map

This document maps the current `cfd_solvers` implementation against the CST-style domain coverage supplied for the continuation work. It is a clean-room capability map, not an assertion of product parity. The goal is to identify which numerical families already exist, which have reference baselines, and which still need production-grade implementation.

Status meanings:

- **Implemented**: a public API plus deterministic regression/validation exists in this repository.
- **Baseline**: a useful numerical reference exists, but the dimensionality, model breadth, or validation set is not yet comparable to a mature specialist solver.
- **Open**: no solver of the requested class exists yet.

## Coverage matrix

| Domain | Current repository state | Status | Principal implementation |
| --- | --- | --- | --- |
| Electrostatics | 2-D FEM electrostatic potential/field solve | Implemented baseline | `src/fem/electromagnetics2d.cpp` |
| Magnetostatics | 2-D out-of-plane vector-potential FEM with reconstructed B | Implemented baseline | `src/fem/magnetostatics2d.cpp` |
| Low-frequency harmonic EM | 2-D harmonic eddy-current formulation | Implemented baseline | `src/fem/eddy_current2d.cpp` |
| Low-frequency transient EM / machines | No moving-band machine formulation, nonlinear B-H transient solver, or coil/circuit machine workflow | Open | Phase 12 |
| Full-wave time-domain EM | 1-D/3-D Cartesian Yee FDTD, cylindrical TM, materials, PML/Mur, ports, probes, NF2FF, SAR | Implemented baseline | `src/fdtd/` |
| Full-wave frequency-domain EM | 1-D complex Maxwell/Helmholtz PEC solver plus 2-D curl-conforming Nedelec Tri3 driven FEM; 3-D edge-element FEM still open | Baseline | `src/em/frequency_domain.cpp`, `src/em/edge_fem2d.cpp` |
| Electromagnetic eigenmodes | New generalized 1-D dielectric-loaded PEC cavity eigensolver; 3-D cavity/waveguide modes still open | Baseline | `src/em/frequency_domain.cpp` |
| Integral-equation EM | Thin-wire EFIE MoM with arbitrary orientation, feeds, loads, junctions, loss and ground image model | Baseline | `src/rf/thin_wire_mom.cpp` |
| MLFMM | No fast multipole acceleration | Open | Phase 12 |
| Asymptotic SBR / physical optics | Optical ray/physical-optics modules exist, but no electromagnetic SBR/PO/RCS solver | Open | Phase 12 |
| Planar/multilayer EM | Microstrip, stripline and CPW quasi-static models; no multilayer full-wave MoM | Baseline | `src/rf/network.cpp` |
| RF/microwave/antennas | N-port networks, Touchstone, transmission lines, waveguide baseline, thin-wire MoM, arrays, matching, efficiency and optimization | Implemented baseline | `src/rf/` |
| Circuit/system EM | SPICE-class MNA, nonlinear devices, AC/transient/noise/PSS/PZ, RF N-port bridges | Implemented baseline | `src/circuit/` |
| Signal integrity | Transmission-line/N-port/Touchstone primitives exist; eye/BER and channel-compliance workflow open | Baseline | `src/rf/`, `src/circuit/` |
| Power integrity | PEEC R/L/C extraction exists; PDN IR-drop/decoupling workflow and surface/volume PEEC open | Baseline | `src/rf/peec.cpp` |
| EMC/EMI | FDTD, antennas, shielding magnitude helper and RF networks exist; dedicated cable/ESD/BCI/lightning/co-site workflows open | Baseline | `src/fdtd/`, `src/rf/` |
| Cable/harness | No frequency-dependent multiconductor harness solver yet | Open | Phase 12 |
| Thermal steady/transient | FEM transient heat, FVM thermal transport, phase change | Implemented baseline | `src/fem/heat2d.cpp`, `src/fvm/thermal_transport.cpp` |
| Conjugate heat transfer | Shared CHT coupling layer over fluid/thermal fields | Implemented baseline | `src/multiphysics/cht.cpp` |
| Structural mechanics | Linear/axisymmetric elasticity, material laws, contact, dynamics and modal baseline | Implemented baseline | `src/fem/` |
| EM -> thermal -> mechanical | Joule/electrothermal and thermoelastic coupling paths exist | Implemented baseline | `src/multiphysics/electro_thermal.cpp`, FEM elasticity/heat |
| Charged-particle tracking | New 3-D non-relativistic Boris Lorentz pusher with time-dependent field callback | Implemented baseline | `src/particle/electromagnetic.cpp` |
| Electrostatic PIC | New periodic 1-D CIC deposition + spectral Poisson + particle push | Implemented baseline | `src/particle/electromagnetic.cpp` |
| Electromagnetic PIC | No self-consistent Yee-grid EM PIC yet | Open | Phase 12 |
| Wakefields | No beam wake-potential/impedance solver | Open | Phase 12 |
| Plasma/discharge | No collision, ionization, secondary-emission, multipactor or corona solver | Open | Phase 12 |
| Bioelectromagnetics | Existing SAR post-processing plus new SAR material helper and implicit 2-D Pennes bioheat coupling | Implemented baseline | `src/fdtd/postprocess.cpp`, `src/multiphysics/bioheat.cpp` |
| Photonics | Sequential/non-sequential optics, diffraction/POP, Gaussian beam, polarization and dispersive glasses; dedicated full-wave photonic eigenmode/band solver open | Baseline | `src/optics/`, `src/fdtd/` |
| Optimization / UQ | Optical optimization/tolerancing, antenna optimizer, circuit Monte Carlo/sensitivity primitives | Implemented baseline | `src/optics/optimization.cpp`, `src/rf/antenna.cpp`, `src/circuit/analysis.cpp` |

## New numerical baselines in this expansion

### Driven frequency-domain Maxwell

`cfd::em::solve_pec_driven_maxwell_1d` solves the transverse 1-D phasor equation with PEC end walls, material permittivity/permeability, conductivity and impressed current density. The implementation uses a complex finite-difference Helmholtz system and reconstructs staggered magnetic field samples. Validation uses a sinusoidal manufactured current whose continuum field amplitude is known analytically.

### Electromagnetic cavity eigenmodes

`cfd::em::pec_cavity_eigenmodes_1d` assembles the generalized dielectric-loaded 1-D curl-curl analogue and solves the symmetric transformed eigenproblem. Regression checks the first three modes and the analytical PEC-cavity fundamental frequency.

### 2-D Nedelec edge-element Maxwell

`cfd::em::solve_driven_edge_maxwell_2d` adds a true curl-conforming vector-FEM baseline on triangular meshes. It uses lowest-order first-kind Nedelec basis functions, deterministic global edge orientation, complex curl-curl/mass/conduction assembly, impressed-current loading and PEC tangential constraints on boundary edges. A manufactured in-plane field is solved on two mesh resolutions and the regression requires the refinement error to decrease.

This remains a reference stepping stone. It closes the 2-D vector-element baseline, but it does **not** close the production 3-D RF FEM requirements for tetrahedral edge elements, wave ports, impedance boundaries, adaptive meshing or large sparse complex solves.

### Charged-particle tracking and electrostatic PIC

`cfd::particle::boris_push` implements the non-relativistic Boris rotation. A uniform magnetic-field regression checks long-run speed conservation. `cfd::particle::ElectrostaticPic1D` adds periodic cloud-in-cell deposition, a spectral periodic Poisson field solve, field interpolation and particle advance. The Poisson kernel is validated against a sinusoidal charge-density field and the PIC loop against a neutralized uniform distribution.

### SAR to Pennes bioheat

`cfd::multiphysics::sar_from_rms_electric_field` converts conductivity, tissue density and RMS electric field to W/kg SAR. `PennesBioheat2D` advances the Pennes equation implicitly with conduction, blood perfusion, metabolic heat and spatial SAR loading. Periodic and fixed-temperature thermal boundaries are supported. A uniform analytical implicit step and fixed-boundary regression validate the implementation.

## Highest-value next implementation order

1. **3-D frequency-domain edge-element Maxwell**: Nedelec basis, complex sparse curl-curl assembly, PEC/PMC/impedance boundaries, lumped and wave ports.
2. **Waveguide/cavity eigenmodes and Q**: reuse the edge-element matrices in a generalized complex eigenproblem.
3. **Cable/harness MTL + EMC coupling**: frequency-dependent RLCG, shield transfer impedance, field-to-wire excitation, SPICE/network export.
4. **Surface-current MoM + MLFMM**: RWG basis first, then acceleration and hybridization with volumetric solvers.
5. **Electromagnetic PIC**: reuse the Yee engine and new Boris particle core for self-consistent charge/current deposition.
6. **Plasma/breakdown**: collisions, secondary emission and ionization on top of PIC.
7. **Wakefields**: bunch excitation, longitudinal/transverse wakes and impedance transforms.
8. **Full bioelectromagnetic feedback**: heterogeneous tissues, SAR averaging transfer, Pennes thermal feedback to EM material properties.

## Explicit non-parity statement

The repository now covers a wide set of the same *physics categories* as a CST-style multiphysics environment, but it is not yet CST-equivalent in solver maturity. The largest remaining gaps are 3-D frequency-domain vector FEM, MLFMM/asymptotic EM, cable/harness EMC, electromagnetic PIC/plasma/wakefields, and production-grade geometry/meshing/adaptivity for those families.
