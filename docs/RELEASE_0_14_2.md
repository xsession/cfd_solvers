# Release 0.14.2 - device tuning, multi-device ranks and adaptive repartition planning

## Highlights

v0.14.2 advances the FluidX3D-class portability baseline from 40/46 to 43/46 tracked capabilities and the full tracker from 594/719 to 597/719.

Closed capabilities:

- kernel-fusion/autotuning database keyed by portable SYCL device identity;
- multiple accelerator devices assignable to one local MPI rank with balanced work slices and queue groups;
- adaptive weighted contiguous domain repartition planning with migration-segment generation and churn control.

The three remaining Phase-2 items are real hardware CI executions on NVIDIA, AMD and Intel GPUs.

## Device-keyed tuning

`cfd/core/autotune.hpp` provides a small persistent database keyed by vendor, device name and driver version. A tuning record stores work-group size, vector width, fused-stage count and measured score. Generic benchmark-driven selection keeps policy outside individual kernels while still allowing solver code to retrieve a per-device choice.

The storage format is a human-readable TSV file with explicit validation. Invalid candidate sizes and non-finite/non-positive benchmark scores fail closed.

## Multi-device rank execution

`cfd/distributed/device_assignment.hpp` now supports groups rather than only one device ordinal per rank. When devices outnumber local ranks, non-overlapping device groups are assigned deterministically; when ranks outnumber devices, the existing deterministic oversubscription behavior remains. Optional per-rank caps allow operators to reserve accelerators for other work.

`split_rank_work(...)` divides a contiguous range across the assigned devices. With SYCL enabled, `make_rank_device_group(...)` constructs one in-order queue per selected accelerator plus the matching slices. This establishes a real multi-device submission boundary without violating SYCL context/USM ownership rules.

## Adaptive repartition planning

`cfd/distributed/repartition.hpp` adds weighted contiguous load balancing for structured/slab ownership. It returns partition boundaries, per-part loads and imbalance, derives the exact old-to-new overlap segments that must migrate, and provides a benefit threshold so small fluctuations do not cause repeated repartitions.

This is a storage-agnostic ownership planner. Solver-specific population/field migration uses the returned segments and remains separate from the load model.

## GPU hardware CI definition

`.github/workflows/gpu-hardware-self-hosted.yml` defines NVIDIA, AMD and Intel jobs using cumulative self-hosted labels. `docs/GPU_HARDWARE_CI.md` documents the runner contract. These jobs are not counted as complete until they execute successfully on physical accelerators.

## Validation

Performed in this environment:

- `cfd-v0142-portability-tests`: passed;
- full default CPU/OpenMP CTest matrix after building all registered targets: passed;
- SYCL API/syntax checks for device-key extraction and multi-device queue assignment: passed against the interfaces used.

This host does not provide real NVIDIA/AMD/Intel SYCL accelerator runners, so hardware qualification is explicitly not claimed.
