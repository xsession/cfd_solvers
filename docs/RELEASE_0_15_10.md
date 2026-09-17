# cfd_solvers v0.15.10 - resident FVM Krylov, outlet pressure and scalar transport

v0.15.10 advances GPU-R3.2 from a Jacobi-only incompressible baseline toward a reusable resident finite-volume equation stack.

## Added

- Device-resident BiCGStab in `SyclCsrLinearAlgebra` for nonsymmetric CSR systems.
- Persistent BiCGStab work vectors; direct device-USM entry points keep volume vectors on the accelerator.
- Device assembly of the first-order Euler/upwind/orthogonal-diffusion momentum CSR numeric values against fixed mesh sparsity.
- Resident incompressible momentum prediction now uses BiCGStab as the primary solver, retaining Jacobi only as a startup/pathology fallback.
- Per-patch fixed-value and zero-gradient pressure descriptors.
- Fixed-pressure/outlet pressure matrix, RHS, face-flux correction and pressure-gradient handling without rebuilding host connectivity.
- Dynamic choice between a pinned all-Neumann pressure reference and an unpinned Dirichlet pressure system.
- `ResidentScalarTransportSycl`: implicit-Euler, first-order-upwind, orthogonal-diffusion scalar transport with fixed/zero-gradient boundaries and optional device volumetric source.
- The scalar solver can consume `ResidentIncompressibleSycl::face_flux_device()` directly when both objects share a SYCL queue/context.

## Residency contract

After setup/boundary uploads and transfer-counter reset, the covered flow + scalar hot path can execute:

`momentum CSR -> BiCGStab(U) -> V/aP -> pressure CG -> U/flux correction -> scalar CSR -> BiCGStab(phi)`

without staging a complete cell or face field through host memory. Krylov convergence still synchronizes scalar reductions with the CPU.

## Numerical scope

The GPU path remains deliberately narrower than the mature CPU FVM implementation:

- momentum: implicit Euler, first-order upwind convection, orthogonal diffusion;
- pressure: variable mobility, explicit non-orthogonal correction, fixed-pressure or pinned all-Neumann baseline;
- scalar: implicit Euler, first-order upwind, orthogonal diffusion, fixed/zero-gradient boundaries;
- BiCGStab currently has no ILU/Jacobi preconditioner;
- the standalone scalar object mirrors mesh geometry in its own resident mesh object even when sharing the flow queue/context;
- turbulence transport, detailed thermophysical coupling, reacting species source integration, moving mesh and distributed FVM halo overlap are still future work;
- physical NVIDIA/AMD/Intel SYCL runtime qualification is not available in the current environment.

## Validation

- Full Python-enabled/default regression matrix: **130/130 CTest targets passed**.
- Python C ABI smoke reports **0.15.10**.
- Host regression verifies the fixed-pressure orthogonal matrix remains symmetric with positive diagonal entries.
- SYCL-only regression covers a nonsymmetric BiCGStab solve, a quiescent fixed-pressure SIMPLE path and resident scalar transport with zero measured hot-loop host field transfer when executed on a real SYCL build.
- Modified SYCL sources and the new v0.15.10 test compile under the strict fake-SYCL C++20 harness with `-Wall -Wextra -Wpedantic -Werror`.

## Tracker

The machine-counted capability total remains **613/719 (85.3%)**. GPU-R3.2 strengthens execution residency for already-tracked FVM capabilities; `SYCL FVM/FEM sparse algebra` stays open until resident FEM assembly/operators and physical accelerator qualification are complete.
