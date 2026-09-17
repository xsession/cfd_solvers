# cfd_solvers v0.15.9 - resident incompressible FVM control loop

v0.15.9 advances GPU-R3 from isolated resident FVM operators and a constant-mobility projection into a complete first-order resident incompressible coupling baseline.

## Added

- On-device momentum diagonal/RHS assembly for Euler time integration, orthogonal viscous diffusion and first-order upwind convection.
- Matrix-free resident Jacobi momentum sweeps, avoiding nonlinear CSR rebuilds on the host.
- Per-cell pressure mobility (`V/aP`) refresh on-device after every momentum predictor.
- Device regeneration of numeric pressure-CSR values against a fixed pinned-Neumann sparsity pattern.
- Device-to-device `SyclCsrLinearAlgebra::update_values_device()` plus host compatibility update API.
- Device-native velocity boundary descriptors for fixed-value, zero-gradient and slip faces.
- Device Gauss pressure gradients and explicit non-orthogonal pressure-flux correction on general `PolyMesh` geometry.
- Resident SIMPLE, PISO and PIMPLE control-loop baseline with persistent velocity, pressure, `H/A`, mobility, face flux, pressure RHS and correction fields.
- Device reductions for continuity, velocity-change and kinetic-energy diagnostics.

## Residency contract

After boundary/setup transfers and a transfer-counter reset, the covered SIMPLE/PISO/PIMPLE hot loop does not stage a cell or face field through host memory. Pressure CG still synchronizes scalar residual reductions with the host for convergence decisions.

## Current scope and limitations

This is intentionally a baseline rather than a production-parity replacement for `CollocatedIncompressible`:

- first-order Euler temporal discretization only;
- first-order upwind convection in the resident momentum path;
- matrix-free Jacobi momentum solution rather than the CPU ILU0/GMRES path;
- pressure uses an all-Neumann domain with a pinned reference cell;
- velocity boundaries support fixed-value, zero-gradient and slip handling;
- fixed-pressure/outlet pressure conditions, turbulence/species/energy equations, moving meshes and distributed halo overlap are not yet wired into this resident loop;
- physical NVIDIA/AMD/Intel SYCL execution is still unqualified in the current environment.

These constraints are documented to avoid claiming full CPU/SYCL feature parity prematurely.

## Validation

- Default build: **129/129 CTest targets passed**.
- Python C ABI smoke reports **0.15.9**.
- New variable-mobility pinned pressure matrix regression passes on a sheared non-orthogonal mesh.
- New SYCL-only regression exercises quiescent SIMPLE, PISO and PIMPLE and asserts zero bulk host-field transfer when built/run on a SYCL device.
- The complete modified SYCL source and v0.15.9 test compile under the strict fake-SYCL C++20 harness with `-Wall -Wextra -Wpedantic -Werror`.

## Tracker

The machine-checkable capability count remains **613/719 (85.3%)**. This release strengthens GPU execution architecture for already-implemented FVM pressure-velocity physics; the combined `SYCL FVM/FEM sparse algebra` gate remains open until resident FEM assembly/operators are implemented and physical accelerator qualification succeeds.
