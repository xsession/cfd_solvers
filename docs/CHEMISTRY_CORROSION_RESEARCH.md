# Open-source chemistry, electrochemistry and corrosion research

Research snapshot: 2026-09-16

The goal is not to concatenate upstream source trees. `cfd_solvers` keeps one C++20 architecture and reimplements well-documented numerical/physical methods, or optionally links to clearly separated third-party libraries when that provides a better long-term result.

## Recommended references

### EchemFEM - LLNL/echemfem

Repository: https://github.com/LLNL/echemfem
License: MIT.

High-value capabilities:

- finite-element electrochemical transport;
- continuous Galerkin and discontinuous Galerkin formulations;
- electroneutral Nernst-Planck;
- diffusion, advection and electromigration;
- electroneutral or Poisson ionic-potential formulations;
- porous and non-porous domains;
- electronic potential in porous electrodes;
- SUPG for advection/migration in its CG formulation;
- generalized modified Poisson-Nernst-Planck finite-size models;
- microkinetic coupling examples.

Use in this project: primary reference for equations, validation cases and the future FEM electrochemistry branch. Its permissive license also makes it legally much simpler than GPL solver code, although clean independent implementation remains preferable for a coherent architecture.

### echemAMR - NatLabRockies/echemAMR

Repository: https://github.com/NatLabRockies/echemAMR
License: BSD-3-Clause.

High-value capabilities:

- 3-D microstructure-resolved electrochemical transport;
- electrode/electrolyte interfacial chemistry;
- immersed-interface representation of complex electrodes;
- nonlinear Butler-Volmer flux conditions;
- AMReX mesh adaptivity;
- MPI parallel execution;
- CUDA GPU execution;
- HYPRE for stiff linear systems.

Use in this project: primary architecture/performance reference for microstructure-scale electrochemistry, immersed interfaces, AMR and GPU/distributed electrochemical transport.

### PHREEQC / PhreeqcRM

Official software page: https://www.usgs.gov/software/phreeqc-version-3
Source mirror: https://github.com/phreeqc-dev/phreeqc3
Status: U.S. Geological Survey public-domain software.

High-value capabilities:

- aqueous speciation;
- equilibrium reactions;
- kinetic reactions;
- mineral precipitation/dissolution;
- surface complexation;
- ion exchange;
- Pitzer and SIT activity models;
- transport/reaction workflows;
- PhreeqcRM reaction module designed to be embedded in transport solvers;
- MPI and OpenMP execution in PhreeqcRM.

Use in this project: aqueous corrosion chemistry/speciation target and optional reaction-module interoperability. This is especially useful for pH, precipitation/corrosion products and complex electrolyte chemistry that should not be reduced to a few hard-coded reactions.

### Reaktoro

Repository: https://github.com/reaktoro/reaktoro
License: LGPL-2.1-or-later.

High-value capabilities:

- chemical equilibrium and kinetics with general constraints;
- arbitrary numbers of phases and species;
- PHREEQC, SUPCRT/SUPCRTBL, NASA and ThermoFun database support;
- automatic differentiation for derivatives;
- modern C++ implementation.

Use in this project: design reference for a generic chemical-system/state API, constrained equilibrium, multiphase chemistry and derivative-aware solvers. Prefer an optional dynamic adapter rather than copying source.

### Cantera

Repository: https://github.com/Cantera/cantera
License: permissive BSD-style 3-clause conditions.

Current repository documentation describes:

- thermodynamic and transport properties;
- chemical equilibrium;
- species production rates;
- large kinetic mechanisms;
- one-dimensional flames;
- reaction-path analysis;
- stirred-reactor networks;
- non-ideal fluids;
- interfaces in C++, Python, C, Fortran and Matlab.

Use in this project: reference for reaction mechanism organization, thermo/transport property interfaces, reactor networks and reacting-flow coupling. A future Cantera-YAML importer would give immediate access to existing mechanisms while keeping our numerical solvers independent.

## Corrosion-specific references

### FEMCorrosionSimulation / FEniC(orr)S

Repository: https://github.com/mrshariati/FEMCorrosionSimulation

This published research code implements a parallel Poisson-Nernst-Planck corrosion model with an algebraic flux-correction/FCT approach, using FEniCS/PETSc/SUNDIALS in the associated workflow.

Important licensing decision: the repository does not expose a normal permissive OSS license in the way EchemFEM/echemAMR do. Treat it as a scientific/paper reference. Do not translate or copy its source into `cfd_solvers` without explicit licensing clarification.

### MOOSE / phase-field ecosystem

Project: https://mooseframework.inl.gov/

Useful concepts include:

- phase-field free-energy formulations;
- multiphase models;
- anisotropy;
- CALPHAD coupling;
- mechanics coupling;
- nucleation and grain evolution;
- automatic differentiation of free-energy functions.

The wider MOOSE ecosystem has included corrosion/oxidation-oriented applications. This is a valuable reference for the later moving-interface/phase-field corrosion phase, not the first electrochemical transport implementation.

### BioDeg

Repository: https://github.com/mbarzegary/BioDeg-UI

BioDeg targets metallic biomaterial biodegradation/corrosion with FreeFEM/C++/Python and supports arbitrary geometry and parallel models. It is a useful validation/problem-selection reference for geometry recession and biomedical corrosion scenarios.

### OpenPNM

Repository: https://github.com/PMEAL/OpenPNM

Useful for future porous-electrode/corrosion-product layers because it provides pore-network transport, reactive transport infrastructure and electrochemical/Butler-Volmer-related models in its ecosystem.

## Selected architecture for cfd_solvers

Rather than picking one chemical codebase, use layers:

1. **chemistry core** - species, phases, reactions, thermo/activity models, kinetics and equilibrium;
2. **electrochemistry core** - electrode reactions, Nernst, Butler-Volmer/Marcus-type kinetics, Faraday laws;
3. **transport** - advection/diffusion/electromigration on FVM and FEM meshes;
4. **potential** - electroneutral, Poisson, solid/electrolyte conduction;
5. **interface physics** - immersed or conformal electrode boundaries, current transfer, double layer;
6. **corrosion state** - dissolution, films, products, pitting/localized attack and geometry recession;
7. **aqueous chemistry adapter** - native equilibrium first, optional PHREEQC/Reaktoro integration for advanced chemistry;
8. **HPC** - use the same OpenMP/SYCL/MPI/AMR runtime as CFD/FDTD rather than a separate execution model.

## Implemented in the first chemistry checkpoint

- species metadata;
- elementary mass-action reaction network;
- Arrhenius forward rates;
- Nernst equation helper;
- Butler-Volmer current-density helper;
- Faradaic molar flux and penetration-rate conversion;
- conservative 1-D Nernst-Planck diffusion/advection/electromigration solver;
- nonlinear 1-D corrosion cell coupling electrolyte ohmic loss and Butler-Volmer interface kinetics;
- analytical/regression tests for reaction rate, Butler-Volmer symmetry, linearized corrosion current and Fourier-mode diffusion decay.

## Next electrochemistry/corrosion validation sequence

1. Binary symmetric electrolyte: Poisson-Nernst-Planck equilibrium.
2. Electroneutral Nernst-Planck benchmark.
3. 1-D electrode with concentration-dependent Butler-Volmer.
4. Tafel polarization curve.
5. galvanic couple with two dissimilar metal boundary patches.
6. diffusion-limited cathodic reaction.
7. pH/speciation coupling through PHREEQC-compatible chemistry.
8. 2-D/3-D metal dissolution with geometry recession.
9. PNP corrosion benchmark derived from published FEniC(orr)S equations/results without copying source.
10. phase-field corrosion and stress-assisted coupling.
