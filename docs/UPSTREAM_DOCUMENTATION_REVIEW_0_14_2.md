# Upstream documentation review - v0.14.2 portability and load balance

This is a clean-room documentation review. It extracts portable interface/design lessons only; no upstream implementation source was copied.

## SYCL device identity and multi-device execution

Reviewed:

- Khronos SYCL device reference: https://github.khronos.org/SYCL_Reference/iface/device.html
- Khronos SYCL platform reference: https://github.khronos.org/SYCL_Reference/iface/platform.html
- Khronos SYCL queue reference: https://github.khronos.org/SYCL_Reference/iface/queue.html
- Khronos SYCL USM reference: https://github.khronos.org/SYCL_Reference/iface/usm_basic_concept.html

Applied lessons:

1. Portable device identity can be built from standard device `vendor`, `name` and `driver_version` queries without depending on a vendor UUID extension.
2. SYCL exposes all root devices through platform/device enumeration, so one host process or MPI rank can own more than one queue/device when hardware is available.
3. USM allocations belong to a context; cross-device algorithms must not silently assume pointers allocated for one context are valid in every queue. v0.14.2 therefore adds a multi-device scheduling boundary rather than pretending existing single-device LBM allocations are transparently shareable.
4. Queue ownership stays explicit. Work slices are associated with concrete device ordinals and each queue can progress independently before an explicit wait.

Implemented:

- device-keyed kernel tuning database with work-group, vector-width and fused-stage choices;
- generic benchmark-driven candidate selection and persistent TSV storage;
- standard-SYCL device key extraction;
- deterministic multi-device groups per local MPI rank;
- balanced contiguous work slicing across each rank's assigned devices;
- SYCL queue-group construction for actual multi-device submission by higher solver layers.

## MPI topology and adaptive repartition boundaries

Reviewed:

- MPI 4.1 standard: https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report.pdf
- MPI 4.1 neighborhood collective/topology section: https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node232.htm

Applied lessons:

1. Communication topology is metadata separate from the numerical kernel, so a load-balancing layer should first produce new ownership boundaries and migration segments.
2. Sparse neighborhood communication benefits from stable topology information; repartition should therefore be triggered deliberately rather than every iteration.
3. Rebalancing decisions need hysteresis/benefit thresholds to avoid partition churn when the expected load improvement is small.

Implemented:

- weighted contiguous partition planner with a configurable minimum cells per partition;
- deterministic zero-weight fallback;
- measured maximum/average load imbalance;
- old/new ownership overlap analysis producing migration segments;
- trigger and minimum-improvement gating for adaptive repartition decisions.

The current baseline plans ownership changes and the exact segments that must migrate. Solver-specific MPI field migration remains an integration concern for each distributed storage layout; the planner is intentionally storage-agnostic.

## GitHub self-hosted hardware CI

Reviewed GitHub Actions documentation for self-hosted runner labels and runner selection. The new workflow targets cumulative labels such as `[self-hosted, linux, x64, gpu-nvidia]` so vendor jobs cannot silently run on generic CPU hosts.

The workflow definitions are committed, but the NVIDIA/AMD/Intel tracker boxes remain open until real hardware jobs pass.
