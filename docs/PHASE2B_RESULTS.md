# Phase 2B results

Phase 2B turns the LBM performance baseline into a physically useful validation baseline.

## Implemented

- Guo-style body acceleration with half-step macroscopic correction.
- Halfway bounce-back walls on D2Q9 CPU.
- Moving-wall bounce-back correction on D2Q9 CPU.
- Left velocity inlet and right pressure/density outlet using Zou-He reconstruction.
- Poiseuille/channel analytical regression.
- Lid-driven cavity regression.
- Velocity-inlet/pressure-outlet channel regression.
- Taylor-Green analytical decay/convergence regression.
- Generic FP32/FP64 periodic pull reference for precision validation.
- Body-force parity between single-grid and two-grid CPU implementations for D2Q9, D3Q19 and D3Q27.
- Equivalent force term added to the optional in-place SYCL implementation, with forced CPU/SYCL parity tests for all three lattices when `CFD_ENABLE_SYCL=ON`; accelerator compilation/execution remains pending because AdaptiveCpp/GPU hardware are unavailable in this environment.
- `set_solid()` now resets any previous moving-wall velocity, making stationary-wall semantics deterministic.

## Numerical snapshot

Development runner: 5 exposed CPU threads, GCC Release build.

### Force-driven Poiseuille flow

CLI case: `cfd-solve lbm-poiseuille`

- grid: 96 x 32
- tau: 0.8
- acceleration: 5e-6 lattice units
- relative L2 velocity-profile error versus the analytical parabola: **2.22e-4**
- maximum speed: **5.62e-3**

The smaller regression case in `cfd-tests` uses 24 x 14 and requires less than 1% relative L2 error; observed error during development was about 1.5e-3.

### Lid-driven cavity

CLI case: `cfd-solve lbm-cavity`

- grid: 96 x 96
- lid speed: 0.05
- tau: 0.8
- center `ux`: **-1.032e-2**
- center `uy`: **1.266e-3**
- maximum fluid speed: **4.844e-2**

The negative center horizontal velocity together with positive near-lid flow is the expected clockwise primary vortex signature.

### Velocity inlet / pressure outlet channel

CLI case: `cfd-solve lbm-channel-io`

- inlet mean `ux`: **0.020000**
- mid-channel mean `ux`: **0.020154**
- outlet mean `ux`: **0.020232**
- maximum speed: **0.030325**

The outlet density regression checks the prescribed value of 1.0 on every non-solid outlet node.

### Taylor-Green convergence

Double-precision D2Q9 periodic pull, 20 time steps, low-Mach initial amplitude 0.001:

| Resolution | relative velocity L2 error |
|---:|---:|
| 16 x 16 | 1.948e-2 |
| 32 x 32 | 4.340e-3 |
| 64 x 64 | 1.122e-3 |

The regression requires the error to decrease by at least roughly 2.86x per refinement step and the 64 x 64 error to remain below 0.2%.

### FP32 versus FP64

For a 40 x 40 Taylor-Green case, tau 0.73, amplitude 0.02, 120 time steps, the maximum difference across density and velocity fields was **4.28e-7** in the development run. The regression gate is 1e-6.

## Scope boundary

The physical wall and Zou-He boundary implementations currently live in the conventional two-grid D2Q9 reference solver. The single-grid CPU and SYCL kernels now support forcing, but boundary-aware in-place streaming is intentionally deferred until the reference cases above can be reproduced exactly enough to avoid hiding boundary errors inside an optimized memory layout.

## Toolchain validation

Validated on 2026-09-16:

- GCC 14 Release + OpenMP + complete CTest numerical suite: pass.
- Clang 17 Release, serial backend + complete CTest numerical suite: pass.
- GCC Release, explicit serial backend + complete CTest numerical suite: pass.
- GCC Debug + AddressSanitizer + UndefinedBehaviorSanitizer: Phase-2B API sanitizer smoke pass.
- AdaptiveCpp/SYCL: source and conditional parity tests updated for forcing, but no AdaptiveCpp compiler/GPU runtime is installed in this environment, so accelerator execution remains pending.
