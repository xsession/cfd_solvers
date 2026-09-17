# Third-party references and provenance

cfd_solvers is an independent implementation. The projects below are capability/architecture references; their source code is not vendored into this repository.

## OpenFOAM 14

- Repository: https://github.com/OpenFOAM/OpenFOAM-14
- License: GNU GPL v3 or later.
- Used as a feature map for finite-volume CFD, solver modularity and multiphysics workflows.

## FluidX3D

- Repository: https://github.com/ProjectPhysX/FluidX3D
- License: custom source license with non-commercial, non-military and other restrictions.
- No source is copied, translated, generated from, or vendored here.
- Performance ideas are implemented independently from public LBM literature and published concepts such as in-place streaming, compact distribution storage, domain decomposition and bandwidth-oriented kernel fusion.
- In-place streaming reference: Moritz Lehmann, “Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs”, Computation 2022, 10(6), 92, DOI 10.3390/computation10060092.
- In-place streaming reference: Moritz Lehmann, “Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs,” Computation 10(6), 92 (2022), DOI: 10.3390/computation10060092.

## Elmer FEM

- Repository: https://github.com/ElmerCSC/elmerfem
- License file identifies GPL 2.0.
- Used as a capability map for FEM, sparse linear systems, adaptivity and coupled multiphysics.

## openEMS

The URL initially supplied for this project, https://github.com/OpenEMS/openems, is the OpenEMS energy-management platform and is not the electromagnetic FDTD solver.

For the electromagnetic capability map this repository uses the established FDTD project:

- Repository: https://github.com/thliebig/openEMS
- License: GNU GPL v3.
- Used only as a feature/validation reference; no source copied.

## Optiland

- Repository: https://github.com/optiland/optiland
- License: MIT.
- Used as a capability map for geometric optics, ray tracing, optical systems, analyses and backend separation.

## AdaptiveCpp

- Project: https://github.com/AdaptiveCpp/AdaptiveCpp
- Used optionally as the SYCL implementation for portable C++ accelerator kernels.

## Contribution rule

Every implementation derived from a restrictive/copyleft reference must be based on mathematical descriptions, standards, publications, public API behaviour, or independently written test cases. Do not use line-by-line translation or close paraphrase of source code.

## Electrochemistry and reactive chemistry references

### EchemFEM

- Repository: https://github.com/LLNL/echemfem
- License: MIT.
- Used as a capability/validation reference for Nernst-Planck transport, electromigration, electroneutral/Poisson potential formulations, porous electrodes and finite-size electrochemical models.

### echemAMR

- Repository: https://github.com/NatLabRockies/echemAMR
- License: BSD 3-Clause.
- Used as a capability/validation reference for microstructure-resolved electrochemistry, immersed electrode interfaces, Butler-Volmer fluxes, AMR and heterogeneous HPC execution.

### Cantera

- Repository: https://github.com/Cantera/cantera
- License: permissive BSD-style license.
- Used as a capability reference for thermodynamic phases, chemical kinetics, transport properties, mechanism handling and reacting-flow coupling.

### Reaktoro

- Repository: https://github.com/reaktoro/reaktoro
- License: LGPL 2.1 or later.
- Used as a capability/reference map for multiphase equilibrium, kinetics, thermochemical databases and automatic-differentiation-based chemistry. Any future adapter must preserve a clean library boundary.

### PHREEQC / PhreeqcRM

- Project: USGS PHREEQC Version 3 and PhreeqcRM.
- Distribution: U.S. Government public-domain software.
- Used as an aqueous-geochemistry/speciation/reactive-transport reference and potential optional external reaction-module integration.

### Corrosion-specific scientific references

- `mrshariati/FEMCorrosionSimulation` is used only as a published Poisson-Nernst-Planck corrosion/FCT validation reference because no clear permissive repository license was identified during the research pass.
- MOOSE phase-field corrosion concepts and OpenPNM porous reactive-transport concepts are tracked as later capability references rather than source donors.

Exact research pins are recorded in `docs/upstreams.json`.

## RF and microwave references

- Palace: https://github.com/awslabs/palace — Apache-2.0; RF FEM/eigenmode/driven/wave-port/HPC capability reference.
- OpenSEMBA FDTD: https://github.com/OpenSEMBA/fdtd — MIT; EMC/FDTD, wire/Huygens/NF2FF capability reference.
- OpenNEC: https://github.com/maurymarkowitz/OpenNEC — MIT; thin-wire antenna/MoM capability reference.
- PEEC and other RF extraction projects with copyleft/unclear licensing are used only as mathematical/behavioral references unless explicitly documented otherwise.

## SPICE and compact-model references

- ngspice — modified-BSD simulator; primary behavioral reference for SPICE analyses and model interoperability.
- Xyce — GPL-3.0; clean-room capability reference for scalable DAE, harmonic balance and Verilog-A workflows.
- QucsatorRF — GPL-2.0; clean-room RF-network/circuit capability reference.
- OpenVAF/OSDI — external Verilog-A/compact-model toolchain/interface target. The core only implements an independent loader/adapter boundary; OpenVAF source is not vendored.
