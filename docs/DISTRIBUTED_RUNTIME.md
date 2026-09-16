# Distributed runtime

Phase 3 provides a distributed-memory runtime shared by CPU/OpenMP and optional SYCL accelerator execution.

## Topology

`cfd::distributed::choose_process_grid()` searches factorisations of the rank count and chooses an aspect-aware 2-D/3-D Cartesian process grid. Uneven global dimensions are supported: each axis distributes the remainder over the lowest process coordinates, so every global cell belongs to exactly one non-empty brick.

`MpiCartesianRuntime` attaches an MPI Cartesian communicator to the chosen grid. A shared-memory communicator discovers `local_rank` and `local_size` for node-local device assignment.

## Correct 3-D halo topology

The structured LBM runtime models all 26 potential 3-D neighbour regions:

- 6 faces;
- 12 edges;
- 8 corners.

D3Q19 needs face and edge neighbours. D3Q27 additionally needs corner neighbours. A one-stage six-face-only exchange would therefore be numerically incomplete for diagonal process-boundary streaming.

## Selective population exchange

Phase 3A exchanged every population stored in every boundary cell. Phase 3B keeps the same audited ghost-cell geometry but sends only populations whose lattice velocity actually crosses that neighbour boundary.

For offset `(ox,oy,oz)`, a population is sent when every non-zero offset component matches that population's lattice-velocity component. Tangential velocity components remain unrestricted. Consequently:

- D3Q19: 5 populations/face cell, 1/edge cell, 0/corner cell;
- D3Q27: 9 populations/face cell, 3/edge cell, 1/corner cell.

`selective_send_entry()` and `selective_receive_entry()` provide one shared `constexpr` linearization used by host and device pack/unpack implementations. This avoids maintaining two subtly different edge/corner layouts.

`MpiSelectiveHaloExchange` allocates compact buffers once and uses persistent `MPI_Send_init` / `MPI_Recv_init` requests. Each time step only repacks the stable buffers, starts the request set, waits, and unpacks.

Use `cfd-halo-plan` to inspect exact traffic:

```bash
./build/cfd-halo-plan --lattice d3q19 --nx 64 --ny 64 --nz 64
./build/cfd-halo-plan --lattice d3q27 --nx 64 --ny 64 --nz 64
```

## CPU overlap

A CPU time step is:

1. pack compact current-boundary populations and start persistent MPI requests;
2. compute the strict local interior with OpenMP while communication is in flight;
3. wait for halo completion and unpack ghost populations;
4. compute the boundary shell;
5. swap current/next population grids.

MPI is called by the main thread only, so `MPI_THREAD_FUNNELED` is sufficient while OpenMP workers execute numerical kernels.

## Device-resident SYCL path

`DistributedSyclPullBlock` stores both distributed population grids in device USM. Collision/pull kernels read halo cells directly from the device-resident current grid and write the next grid. Only checkpointing/macroscopic downloads intentionally copy the full local state to the host.

`MpiSyclSelectiveHaloExchange` has two transports:

### Pinned-host staged mode (default)

1. pack compact halos on the device;
2. copy device send buffers into `sycl::malloc_host` pinned host buffers;
3. start persistent MPI requests using those host pointers;
4. execute strict-interior collision on the device while MPI is in flight;
5. wait for MPI;
6. copy host receive buffers to compact device receive buffers;
7. unpack on-device, then execute the boundary shell.

This mode requires no accelerator-aware MPI implementation.

### Direct device-buffer mode

`--gpu-aware-mpi` creates persistent MPI requests directly against compact device-USM buffers. The application still waits for device packing before starting MPI and uses separate communication buffers so the strict-interior kernel can overlap safely with network activity.

This is an explicit opt-in because MPI accelerator-awareness is implementation/backend dependent. Users must validate their MPI/UCX/libfabric stack for the SYCL backend in use.

## Device assignment

`device_assignment.hpp` maps `local_rank` deterministically onto visible accelerator devices. One accelerator per local rank is preferred. If more ranks than devices are visible, assignment wraps and reports oversubscription.

## Checkpoint/restart

Both CPU and distributed SYCL reference blocks use the same rank-local Phase-3 checkpoint payload:

- magic/version/endian marker;
- descriptor population count;
- global dimensions;
- brick global begin and local extent;
- time step;
- relaxation and acceleration parameters;
- interior populations in component-major order.

Restart currently requires the same decomposition. A decomposition-independent collective format is still planned.

## Build and run

CPU MPI:

```bash
cmake -S . -B build-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_ENABLE_OPENMP=ON \
  -DCFD_ENABLE_MPI=ON
cmake --build build-mpi --parallel
mpirun -np 4 ./build-mpi/cfd-distributed \
  --backend cpu --lattice d3q27 --nx 192 --ny 128 --nz 96 --steps 100
```

MPI + SYCL, staged transport:

```bash
cmake -S . -B build-mpi-sycl \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_ENABLE_OPENMP=ON \
  -DCFD_ENABLE_MPI=ON \
  -DCFD_ENABLE_SYCL=ON \
  -DAdaptiveCpp_DIR=/path/to/AdaptiveCpp/lib/cmake/AdaptiveCpp \
  -DACPP_TARGETS=generic
cmake --build build-mpi-sycl --parallel
mpirun -np 4 ./build-mpi-sycl/cfd-distributed \
  --backend sycl --lattice d3q19 --nx 256 --ny 256 --nz 256 --steps 200
```

Direct accelerator-aware MPI is the same command plus `--gpu-aware-mpi`.

## Scaling benchmark

Strong scaling keeps the global mesh fixed:

```bash
MODE=strong RANKS_LIST="1 2 4 8" REPEATS=3 \
BACKEND=cpu LATTICE=d3q19 \
NX=256 NY=256 NZ=256 \
EXE=./build-mpi/cfd-distributed \
./scripts/benchmark_distributed.sh
```

Weak scaling keeps total cells per rank constant by multiplying one global axis by the rank count. The default is the x axis and can be changed with `WEAK_AXIS=y` or `WEAK_AXIS=z`:

```bash
MODE=weak WEAK_AXIS=x RANKS_LIST="1 2 4 8" REPEATS=3 \
BACKEND=cpu LATTICE=d3q19 \
NX=64 NY=128 NZ=128 \
EXE=./build-mpi/cfd-distributed \
./scripts/benchmark_distributed.sh
```

For a GPU-aware MPI stack, set `BACKEND=sycl GPU_AWARE_MPI=1` and point `EXE` at the MPI+SYCL build. The benchmark script emits stable CSV columns including process-grid shape, total messages per step, total halo bytes per step, elapsed seconds, MLUPS, transferred GiB and mass.
