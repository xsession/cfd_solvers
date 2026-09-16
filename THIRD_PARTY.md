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
