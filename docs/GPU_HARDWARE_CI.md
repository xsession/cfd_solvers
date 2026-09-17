# GPU hardware CI runner contract

The repository contains `.github/workflows/gpu-hardware-self-hosted.yml` for real accelerator validation. The workflow intentionally targets self-hosted runners because the project needs direct SYCL access to NVIDIA, AMD and Intel GPUs rather than CPU-only emulation.

## Required runner labels

Register Linux x64 runners with the normal `self-hosted`, `linux`, `x64` labels plus exactly one accelerator label:

- `gpu-nvidia`
- `gpu-amd`
- `gpu-intel`

GitHub Actions matches all labels in a `runs-on` array, so these jobs will not accidentally land on generic CPU runners.

## Provisioning contract

Each runner must provide:

- a working vendor driver/runtime visible to the service account;
- AdaptiveCpp and `acpp-info` on `PATH`;
- CMake and a C++20 compiler compatible with the installed AdaptiveCpp toolchain;
- enough local storage for a Release build;
- no container/device filtering that hides the intended accelerator.

The workflow prints `acpp-info`, configures `CFD_ENABLE_SYCL=ON`, builds the common HPC, LBM and portability regressions, and executes them on hardware.

## Qualification rule

The three Phase-2 hardware-CI tracker items remain unchecked until the corresponding vendor job has completed successfully on a real runner. Merely having the workflow file is not treated as hardware qualification.

## Security note

Self-hosted runners execute repository workflow code with access to the host. Restrict runner access appropriately, keep drivers/toolchains patched, and do not expose privileged credentials to untrusted pull-request code.
