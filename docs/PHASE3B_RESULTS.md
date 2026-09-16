# Phase 3B results - device-resident distributed transport

Phase 3B moves the distributed LBM reference path toward heterogeneous multi-device execution while keeping the CPU/MPI reference intact.

## Implemented

- Selective D3Q19/D3Q27 halo plans that send only populations crossing each process boundary.
- Shared host/device `constexpr` mapping for face/edge/corner payload linearization.
- Persistent MPI requests for compact CPU halo traffic.
- Device-resident distributed two-grid D3Q19/D3Q27 pull block using SYCL device USM.
- On-device selective halo pack/unpack kernels.
- Pinned-host staging using SYCL host USM as the default portable MPI transport.
- Explicit `--gpu-aware-mpi` direct-device transport mode for accelerator-aware MPI stacks.
- One device per node-local MPI rank assignment.
- Strict-interior GPU work scheduled while MPI halo traffic is in flight.
- GPU checkpoint/restart using the Phase-3 rank-local checkpoint format.
- `cfd-halo-plan` exact communication-volume calculator.
- `scripts/benchmark_distributed.sh` strong/weak rank-scaling benchmark driver with CSV output.
- `cfd-distributed --csv` reporting MLUPS plus aggregate halo messages/bytes per step for reproducible scaling analysis.
- D3Q19/D3Q27 MPI/SYCL staged parity tests compiled whenever both MPI and SYCL are enabled.

## Compact halo traffic

For a representative 64x64x64 local brick using FP32 populations:

| Lattice | Full-cell halo bytes/step | Selective bytes/step | Reduction |
|---|---:|---:|---:|
| D3Q19 | 1,926,752 | 494,592 | 74.33% |
| D3Q27 | 2,738,016 | 893,984 | 67.35% |

These are exact layout counts from `cfd-halo-plan`; they are not performance-counter estimates.

D3Q19 sends five populations per face cell and one per edge cell. Corner messages are omitted because D3Q19 has no corner lattice directions. D3Q27 sends nine populations per face cell, three per edge cell, and one per corner cell.

## Validation completed in the development environment

- GCC/OpenMP CPU build and full regression suite.
- Selective virtual-rank exchange is bitwise-identical to the Phase-3A full-cell halo exchange for D3Q19 and D3Q27.
- The compact layout is exercised on the uneven 17x15x13, 2x2x2 virtual decomposition.
- MPI-only translation units compile in a local MPI API syntax harness.
- Combined MPI+SYCL Phase-3B headers, tests and CLI compile in a local API syntax harness.
- The distributed benchmark driver passes a mocked strong/weak-scaling CSV smoke test, including rank-dependent weak-scaling dimensions.

## Release validation

The CPU-side Phase-3B implementation passes GCC/OpenMP, Clang serial fallback, serial GCC and ASan+UBSan smoke validation. The Phase-3 -> Phase-3B patch is also applied to a fresh Phase-3 source tree and rebuilt/tested as part of release packaging.

## Hardware validation still required

This development environment does not contain a real MPI installation, AdaptiveCpp, or accelerator runtime. Therefore the following remain hardware/CI gates rather than claimed measurements:

- actual persistent MPI execution;
- pinned-host staged SYCL+MPI execution;
- direct device-buffer GPU-aware MPI execution;
- multi-GPU strong/weak scaling;
- NVIDIA, AMD and Intel accelerator numerical parity/performance.

The normal OpenMPI CI job executes the compact persistent CPU path. If CI is later provisioned with AdaptiveCpp and GPUs, `cfd-mpi-tests` automatically adds staged SYCL distributed parity coverage.
