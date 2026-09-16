# GPU backend strategy

The accelerator language remains C++20 and SYCL is the portability API. AdaptiveCpp is the first implementation target.

## Memory model

Hot solver state uses device USM (`sycl::malloc_device`) so population arrays remain resident on the selected accelerator across time steps. Host USM (`sycl::malloc_host`) is used only for the default MPI staging path and behaves like pinned/device-accessible host memory.

The distributed runtime never copies the entire local lattice merely to exchange a halo. Phase 3B first packs only boundary-crossing populations into compact device buffers.

## Distributed transports

### Staged host - portable default

Compact device send buffers are copied into pinned host buffers before MPI. Received pinned-host buffers are copied back into compact device buffers and unpacked by a device kernel. This works with an ordinary host-buffer MPI implementation.

### Direct device MPI - explicit opt-in

When the selected MPI stack is accelerator-aware, `--gpu-aware-mpi` passes the compact device-USM buffers directly to persistent MPI point-to-point requests. This removes the host staging copies, but support depends on the MPI transport and the accelerator backend behind SYCL.

Do not infer support only because the code compiled. Validate the MPI stack and backend on the target machine.

## Overlap

The distributed device time step uses separate compact send/receive buffers, allowing this order:

1. device halo pack;
2. start MPI;
3. device strict-interior collision while MPI is active;
4. finish MPI and device halo unpack;
5. device boundary-shell collision;
6. population-grid swap.

The current implementation uses an in-order SYCL queue so staging copies, unpack and boundary kernels have deterministic dependencies without accidental implementation-specific scheduling.

## Device mapping

Node-local MPI rank selects a visible accelerator ordinal deterministically. External visibility controls can still be used to map ranks to NUMA-local GPUs. Oversubscription is reported rather than silently hidden.

## Numerical policy

- D3Q19/D3Q27 distributed GPU kernels retain the conventional two-grid pull algorithm as the distributed correctness reference.
- The single-grid in-place scheme remains the local high-performance path and will be introduced into distributed execution only after the two-grid GPU/MPI parity matrix is hardware-validated.
- CPU reference solvers remain available for field-level comparisons.
- Fast-math and reduced-precision storage stay opt-in.

## Current validation gap

The development container has GCC/Clang but no AdaptiveCpp or physical accelerator runtime. The MPI+SYCL code is source/syntax validated locally, and staged distributed parity tests are compiled whenever both features are enabled, but execution still requires GPU-capable CI or target hardware.
