# cfd_solvers v0.15.8 - GPU-R3 resident FVM infrastructure

v0.15.8 implements the first device-owned finite-volume mesh/operator/pressure-correction path. The goal is not to claim that the production SIMPLE/PISO/PIMPLE solver is fully GPU-resident yet; instead, this release establishes the mesh, field-operator and sparse-solver seams required to keep a complete pressure-correction loop on the accelerator once resident momentum assembly is added.

## Added

- `cfd::fvm::ResidentPolyMeshSycl`
  - immutable `PolyMesh` owner/neighbour topology mirrored to device USM;
  - flattened cell-to-face adjacency;
  - resident cell geometry and face geometry;
  - precomputed interpolation, face-normal and orthogonal pressure metrics;
  - explicit setup/output upload/download helpers with transfer accounting.
- Device FVM operators using cell-owned or face-owned kernels:
  - scalar face interpolation;
  - Gauss scalar gradient;
  - vector divergence;
  - orthogonal scalar Laplacian;
  - predictor face flux;
  - face-flux divergence;
  - pressure RHS formation;
  - pressure face-flux correction;
  - cell velocity correction from a resident pressure gradient.
- `build_orthogonal_pressure_matrix()` for constant or cell-varying positive mobility.
  - fixed-value boundary mode produces an SPD Dirichlet operator;
  - all-Neumann mode can pin a symmetric zero-pressure reference cell without destroying CG compatibility.
- `SyclCsrLinearAlgebra(const CsrMatrix&, sycl::queue)` so a resident FVM context and the persistent CSR/Krylov solver can share the exact same SYCL queue/context and USM allocations.
- `cfd::fvm::ResidentPressureProjectionSycl` as a constant-mobility closed-domain projection baseline. `H/A`, pressure, gradient, corrected velocity, face flux, RHS, mobility and continuity work storage all remain device-owned through `project()`.

## Residency contract

After setup and `reset_transfer_stats()`, the covered pressure-correction sequence is:

1. interpolate resident `H/A` to faces and form predictor flux;
2. assemble the resident continuity RHS;
3. solve the pinned pressure system with direct-device CG;
4. compute the pressure gradient on-device;
5. correct velocity and face flux on-device;
6. optionally reduce continuity to one host-visible scalar.

No volume field is copied through host memory inside that block. The present CG implementation still synchronizes scalar dot-product/convergence values with the host; those synchronizations are tracked separately from bulk transfer bytes.

## Validation

- Default CPU/HDF5/Python build: **128/128 CTest targets passed**.
- New CPU regression validates:
  - fixed-pressure orthogonal matrix symmetry/SPD;
  - symmetric pinned-Neumann pressure reference handling.
- Strict fake-SYCL C++20 syntax validation passes with `-Wall -Wextra -Wpedantic -Werror` for:
  - resident FVM mesh/operators;
  - queue-sharing sparse linear algebra;
  - resident pressure-projection API and regression source.
- Python ABI smoke reports **0.15.8**.

A physical SYCL accelerator is not available in the development environment. Therefore v0.15.8 does **not** claim NVIDIA/AMD/Intel runtime parity, performance, queue-overlap quality, or full production SIMPLE/PISO/PIMPLE GPU qualification.

## Remaining GPU-R3 work

- resident momentum-equation assembly and coefficient update;
- non-orthogonal/corrected gradient and Laplacian kernels;
- boundary-condition descriptors consumed entirely on-device;
- variable-mobility pressure-matrix coefficient refresh without host rebuild;
- resident turbulence/species/energy equation assembly;
- distributed FVM halo exchange and interior/boundary overlap;
- real-device numerical parity and performance qualification.

The combined tracker gate `SYCL FVM/FEM sparse algebra` remains open because v0.15.8 completes the FVM pressure side but not resident FEM assembly/operator integration.
