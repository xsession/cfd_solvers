# cfd_solvers v0.15.6 - resident LBM multiphysics GPU-R1

v0.15.6 completes the first multiphysics extension of the device-residency architecture introduced in v0.15.5. It is an execution-quality release: the integration tracker remains 612/719 because the underlying thermal/free-surface/particle/Q-criterion capabilities were already present on the CPU path.

## Implemented

- Persistent three-component device-local LBM acceleration field.
- Direct device-force access for coupled physics; optional host upload remains a compatibility boundary.
- Resident D3Q7 thermal/passive-scalar lattice with D3Q19 velocity coupling.
- Resident scalar reduction and Boussinesq acceleration coupling.
- Resident conservative periodic VOF/free-surface transport.
- Resident interface normals, curvature and continuum-surface-force acceleration.
- Resident immersed-boundary particle state and push.
- Trilinear velocity interpolation directly from resident macroscopic fields.
- Atomic two-way equal-and-opposite particle reaction spreading to the LBM acceleration field.
- Resident 3-D Q-criterion generation.

## Intended hot-loop sequence

1. clear the persistent device acceleration field;
2. reconstruct D3Q19 macroscopic values into device scratch;
3. advance D3Q7 thermal transport and accumulate buoyancy;
4. advect free-surface fill and accumulate capillary acceleration;
5. advance particles and atomically accumulate reaction acceleration;
6. execute D3Q19 stream/collide using the combined local force;
7. repeat without downloading complete fields.

Scalar reductions may synchronize. Setup uploads and explicit output/checkpoint/download calls are host boundaries by design.

## Validation

- Default Release/OpenMP/HDF5/Python matrix: **126/126 CTest targets passed**.
- New resident SYCL headers and v0.15.6 test source pass strict C++20 fake-SYCL syntax compilation with `-Wall -Wextra -Wpedantic -Werror`.
- No physical SYCL accelerator is available in the development environment; NVIDIA/AMD/Intel numerical parity and throughput qualification remain open and are not claimed by this release.

## Next GPU work

- device-native visualization/export of thermal/fill/Q fields;
- particle binning and local force accumulation after real-hardware contention profiling;
- physical multi-vendor bandwidth/transfer qualification;
- GPU-R2 resident FDTD fields and probes;
- GPU-R3 FVM mesh/field/operator residency.
