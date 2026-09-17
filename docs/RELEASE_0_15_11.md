# cfd_solvers v0.15.11 - shared resident FVM fields, RANS and thermo/species sources

v0.15.11 advances GPU-R3.3 from a flow-plus-generic-scalar baseline to a shared resident equation environment for turbulence and reacting thermal/species transport.

## Added

- `ResidentFvmFieldRegistrySycl`: one named device-field registry layered on a single `ResidentPolyMeshSycl` mesh/queue/context.
- `ResidentScalarEquationSycl`: reusable implicit-Euler, first-order-upwind scalar equation with device cell diffusivity, explicit volumetric source and semi-implicit sink fields.
- Variable-diffusivity resident scalar matrix assembly with face interpolation and fixed-value/zero-gradient boundaries.
- Device scalar clamping and velocity strain-rate evaluation built from resident Gauss gradients.
- Resident Spalart-Allmaras transport with wall-distance, production/destruction and eddy-viscosity fields.
- Resident k-epsilon transport with segregated k/epsilon solves, production limiting, positivity floors and bounded eddy viscosity.
- Resident k-omega SST transport with F1/F2 blending, cross-diffusion, DES length-scale switching and resident eddy viscosity.
- `ResidentThermoSpeciesSourceSycl`: one-step Arrhenius reaction-rate, fuel-consumption and heat-release source generation on-device.
- Direct device-to-device source application into shared temperature/species equations.

## Residency architecture

The main change is ownership rather than just another set of kernels. A flow solver can expose its existing `ResidentPolyMeshSycl`; turbulence, temperature and species equations then attach to that same mesh/queue through `ResidentFvmFieldRegistrySycl`.

The covered execution chain is now:

`resident U/phi -> strain/gradients -> RANS closures -> RANS BiCGStab solves -> reacting source kernel -> T/Y BiCGStab solves`

No extra mesh mirror is required for these new shared equations and no complete volume field needs to cross the host boundary in the covered hot loop.

## Numerical scope and limitations

- Resident RANS currently uses implicit Euler, first-order upwind convection and orthogonal diffusion.
- The GPU SA/k-epsilon/SST implementations are clean-room ports of the project's existing CPU model equations; they are not bitwise-equivalent discretizations because the resident path uses Gauss gradients while the CPU reference uses its mature least-squares path in several places.
- Wall-function treatment and the full CPU boundary-condition catalogue are not yet resident.
- The SST path includes the existing DES baseline but does not add new DDES/IDDES physics.
- The reacting source is a compact one-step Arrhenius coupling seam, not a replacement for the project's full CPU chemistry mechanism/integrator stack.
- Krylov convergence still synchronizes scalar reductions with the host.
- Physical NVIDIA/AMD/Intel runtime and performance qualification remain open in this environment.

## Validation

- Full Python-enabled/default regression matrix: **131/131 CTest targets passed**.
- Python C ABI smoke reports **0.15.11**.
- The new default regression keeps the CPU k-epsilon reference convergence/positivity path exercised.
- The SYCL-only regression covers one shared resident mesh/registry with SA, k-epsilon, SST, temperature, fuel and Arrhenius heat/species source coupling when built on a real SYCL target.
- Modified/new accelerator sources and the new test compile with the fake-SYCL C++20 harness under `-Wall -Wextra -Wpedantic -Werror`.

## Tracker

The machine-counted capability total remains **613/719 (85.3%)**. GPU-R3.3 increases accelerator residency/parity for already-tracked FVM turbulence/thermal/species capabilities. The combined `SYCL FVM/FEM sparse algebra` item remains open until resident FEM assembly/operators and physical accelerator qualification are complete.
