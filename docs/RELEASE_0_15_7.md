# cfd_solvers v0.15.7 - resident FDTD GPU-R2

## Scope

v0.15.7 applies the device-resident execution policy to the 3-D structured Yee FDTD family. This release is a clean-room implementation and does not translate openEMS or FluidX3D source code.

## Added

- `cfd::fdtd::ResidentMaxwell3DSycl`;
- persistent SYCL E/H fields;
- resident anisotropic permittivity, permeability and conductivity coefficients;
- 3-D CPML with per-axis graded profiles and twelve convolution-memory fields;
- device-side Gaussian initialization, soft source injection and material-box updates;
- device scalar energy reduction and direct device probe capture;
- six-face device halo pack/unpack;
- physical-face masks for decomposition-ready internal faces;
- `cfd-v0157-gpu-fdtd-residency-tests`.

## Residency contract

After setup, stepping, source injection, material use, probe capture and halo packing do not require whole-field host transfers. Output/checkpoint extraction remains an explicit boundary operation. Scalar convergence/diagnostic observations may synchronize.

## Validation

- default project matrix: **127/127 CTest targets passed**;
- strict fake-SYCL C++20 syntax build: PASS with `-Wall -Wextra -Wpedantic -Werror`;
- Python C ABI smoke: v0.15.7 PASS.

No physical SYCL accelerator is available in the development environment. Real-device numerical parity, CPML reflection, bandwidth and GPU-aware MPI qualification therefore remain open and are not implied by this release.

## Tracker impact

The Phase-5 `SYCL FDTD kernels` item is now complete. Overall progress becomes **613/719 (85.3%)** and Phase 5 becomes **23/26 (88.5%)**. `MPI + SYCL domain decomposition` remains open.

## Next

GPU-R3 should introduce a resident finite-volume mesh/field/operator layer so FVM gradients, fluxes, pressure correction and sparse solves can share device-owned state instead of moving vectors through host memory.
