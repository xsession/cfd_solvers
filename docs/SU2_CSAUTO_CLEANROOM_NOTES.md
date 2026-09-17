# SU2 and csauto clean-room notes for v0.10.9

This note records what was learned from public repository inspection and what was implemented independently in `cfd_solvers`.

## Sources inspected

- `su2code/SU2` default branch `master`.
- `simvia-tech/csauto` default branch `main`.

No source code was copied or mechanically translated. SU2 is LGPL-2.1 licensed and csauto is GPL-3.0 licensed, so both are treated as clean-room architecture and workflow references for this MIT repository.

## SU2 lessons applied

SU2 is organized as a C++ PDE/CFD and optimization code with separate solver, iteration, numerics, output and Python workflow layers. The public tree shows dedicated C++ directories such as `drivers`, `iteration`, `numerics`, `solvers`, `variables` and Python-side workflows including continuous/discrete adjoint, direct differentiation, finite differences and uncertainty scripts.

Transferable clean-room ideas implemented here:

- keep primal solver kernels separate from campaign/optimization orchestration;
- represent a solver's externally visible capabilities with a small descriptor rather than hard-wiring workflow code to solver names;
- add simple finite-difference gradient and gradient-descent helpers as validation-facing building blocks before attempting full adjoints;
- keep gradient utilities regression-testable on tiny analytic objectives.

Implemented files:

- `include/cfd/workflow/campaign.hpp`
- `src/workflow/campaign.cpp`
- `tests/test_v0109_transport_campaign.cpp`

## csauto lessons applied

csauto demonstrates a campaign layer around CFD solvers: one DOE row generates one case, template files use placeholders and conditional blocks, runtime-specific logic is separated from solver-specific logic, and solver capabilities are provided through an adapter boundary.

Transferable clean-room ideas implemented here:

- factorial and Latin-hypercube DOE generation with deterministic seeds;
- `{placeholder}` rendering plus `<!-- IF key -->`, `<!-- IF key=value -->` and `<!-- IF key!=value -->` conditional blocks;
- `SolverAdapterDescriptor` and `SolverCapability` as a small boundary for future external solver adapters;
- in-memory campaign registry updates and summaries to support future dashboard/export work.

## Current limitations

- No external solver launching yet.
- No persistent campaign directory writer yet.
- No Slurm/native/container runtime integration yet.
- No web dashboard yet.
- No true continuous/discrete adjoint implementation yet.

These limitations are deliberate: the current release validates the reusable in-process primitives first.
