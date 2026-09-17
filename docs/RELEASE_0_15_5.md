# cfd_solvers v0.15.5 - GPU-resident execution groundwork

v0.15.5 is an optimization/architecture release. It does not close a new solver-family capability checkbox; instead it makes the existing SYCL paths follow a stricter device-residency contract derived independently from public accelerator/HPC literature and FluidX3D's publicly documented design patterns.

## Changes

- Added `cfd::core::DeviceTransferStats` for explicit host/device/device transfer bytes and synchronization accounting.
- Changed common SYCL LBM uniform and Taylor-Green initialization to generate populations directly on-device.
- Added persistent device macroscopic scratch to the in-place Esoteric-Pull SYCL solver.
- Added device-side mass reductions so benchmark/CLI mass checks do not download population lattices.
- Changed macroscopic host extraction to transfer only `rho, ux, uy, uz` instead of all q populations.
- Added resident-byte and transfer-stat inspection APIs to the SYCL LBM paths.
- Reused persistent device Krylov vectors in the SYCL CSR solver instead of allocating/freeing work vectors for every multiply/CG solve.
- Added direct device-USM SpMV and CG entry points so a future resident FVM/FEM field layer can bypass host vector staging entirely.
- Added a focused no-host-transfer LBM residency regression and fake-SYCL syntax validation path.
- Added `FLUIDX3D_GPU_RESIDENCY_RESEARCH.md`, including the clean-room/license boundary and a solver-family roadmap.

## What this release does not claim

- It does not make the complete FVM, FEM, chemistry or workflow stack GPU-resident.
- It does not claim physical GPU qualification from the development container; physical hardware CI remains required.
- It does not copy or translate FluidX3D source code. FluidX3D is source-available non-commercial rather than open source, so this project keeps an explicit clean-room boundary.
- It does not add a new tracker checkbox; integration remains at the v0.15.4 capability count.

## Device-residency regression

The focused regression resets transfer accounting after initialization, performs multiple D3Q19 time steps and device-side mass reductions, and requires zero explicit host-transfer bytes during that interval. An explicit macroscopic download is then checked against the expected `4*N*sizeof(float)` transfer budget.

The test still counts synchronization points because a scalar reduction/result observation can serialize the queue even when it does not copy a volume field.

The default CPU regression configuration with Python ABI smoke coverage passes **125/125 CTest targets**. The SYCL-specific sources additionally pass the strict local fake-SYCL syntax harness; this checks source compatibility only and is not a substitute for physical accelerator execution.

## Next accelerator work

The next highest-value residency work is to keep the LBM thermal/free-surface/particle extensions and derived visualization fields on-device. In parallel, a GPU `PolyMesh`/field view is required before OpenFOAM-class FVM can truthfully claim a fully resident SIMPLE/PISO/PIMPLE loop.
