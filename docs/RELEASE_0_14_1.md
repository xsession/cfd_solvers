# Release 0.14.1 - LBM physics, geometry, precision and diagnostics wave

## Highlights

v0.14.1 advances the FluidX3D-class LBM phase from 29/46 to 40/46 tracked capabilities and moves the full integration tracker from 583/719 to 594/719 capabilities.

Closed capabilities:

- D2Q9 cumulant collision baseline;
- Smagorinsky-Lilly LES in the in-place/pull solver;
- conservative VOF/free-surface baseline with phase classification, interface normals/curvature and surface-tension forcing;
- link-wise interpolated curved-wall bounce-back;
- immersed-boundary point particles with passive and two-way coupling;
- optional local LBM body-acceleration field for distributed particle/capillary reaction forces;
- optional SYCL triangle-mesh voxelization into the common solid mask;
- compressed 16-bit population storage with explicit conservation/macroscopic error gates;
- IEEE FP16, BF16, shifted-half and independent fixed-deviation packed storage experiments;
- 3-D Q-criterion derived field;
- dependency-free interactive HTML scalar viewer kept outside the solver core.

## Collision and free-surface models

`CumulantD2Q9Solver` provides a readable low-Mach cumulant-space baseline. Second-order cumulants retain viscosity control while higher-order cumulants relax toward factorized equilibrium. It is validated with a periodic Taylor-Green decay case.

`EsotericPullSolver` can enable a local Smagorinsky relaxation-time correction based on the non-equilibrium stress tensor. Existing configuration remains source-compatible because the option is appended with a disabled default.

`free_surface.*` adds a conservative periodic VOF transport baseline, fill-fraction phase classification, interface normals, curvature and a continuum-surface-force acceleration field. This release does not claim analytic PLIC reconstruction or full FluidX3D free-surface parity.

## Geometry and particles

D2Q9 supports per-link wall fractions for interpolated curved boundaries. The half-way distance reduces to the existing bounce-back behavior in regression coverage.

The particle layer interpolates fluid velocity to Lagrangian particles through a compact regularized kernel. Passive mode moves particles without reacting on the fluid; two-way mode spreads the equal-and-opposite drag into the local acceleration field. The focused test verifies force conservation explicitly.

An optional SYCL voxelization path performs triangle parity tests per voxel and emits the same mask as the CPU geometry path. Accelerator runtime parity remains a conditional test because this packaging host has no SYCL device runtime.

## Population-memory experiments

`compressed_pull.hpp` keeps collision arithmetic in float while population memory uses 16-bit codecs. The focused regression compares compressed results against the float in-place solver and enforces explicit mass/macroscopic-field tolerances for each codec rather than treating reduced memory as automatically acceptable.

The encodings are clean-room experiments and do not claim bit compatibility with another solver's custom compression format.

## Diagnostics and viewer

`q_criterion_3d(...)` computes the incompressible Q-criterion from the velocity-gradient tensor. Analytic solid-rotation and simple-shear regressions validate sign/magnitude behavior.

`write_html_scalar_viewer(...)` emits a standalone canvas-based HTML viewer with pan, zoom and hover values. It has no runtime dependency and stays in the visualization layer rather than the numerical kernel.

## Validation

Executed together on the release source tree:

- `cfd-tests` - passed;
- `cfd-completion-tests` - passed;
- `cfd-v0140-hpc-runtime-tests` - passed;
- `cfd-v0141-lbm-les-qcriterion-tests` - passed.

The optional SYCL voxelizer also passed a C++ API/syntax check against the SYCL interface used. No real accelerator was available here, so GPU execution is not claimed.

## Remaining Phase-2 items

Six intentionally unchecked items remain:

- hardware CI on NVIDIA;
- hardware CI on AMD;
- hardware CI on Intel GPU;
- kernel-fusion/autotuning database by device;
- multiple GPUs per rank;
- adaptive domain repartitioning.
