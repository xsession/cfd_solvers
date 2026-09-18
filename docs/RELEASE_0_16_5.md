# cfd_solvers v0.16.5 - GPU-aware distributed DEM transport

v0.16.5 connects the resident SYCL DEM state introduced in v0.16.2 with the distributed slab/contact architecture completed on the CPU path in v0.16.4. The release focuses on preserving particle/contact identity and moving only compact migration/halo records rather than copying the full DEM state through host memory.

## Stable device identity and history

`ResidentDemParticle` now carries a stable 64-bit `global_id` and an `owner_rank`. Existing single-device callers that leave the ID unset receive deterministic local IDs on upload, while distributed callers should provide globally unique IDs.

The resident Mindlin history table is now keyed by the neighbor's **global ID** rather than its transient local array index. Pair ownership also follows global-ID ordering. Therefore particle reordering, compaction and rank migration do not invalidate contact identity.

New packed device PODs are provided:

- `ResidentDemPackedParticle`;
- `ResidentDemPackedHistory`.

`ResidentDemSycl::pack_indices_device()` packs selected particles and all fixed history slots on-device. `append_packed_device()` restores them directly into resident storage, and `retain_x_slab_device()` compacts owned state on-device after outgoing migrants have been packed.

## MPI/SYCL exchange modes

`MpiResidentDemSyclExchange` adds a slab-transport layer with two explicit modes.

### staged_host

Only selected compact migrant/ghost payloads cross the PCIe/host boundary:

1. classify migration/halo indices on-device;
2. pack compact records on-device;
3. copy only packed records to host;
4. exchange with `MPI_Alltoallv`;
5. copy only received compact records back to device;
6. compact resident ownership and append incoming migrants on-device.

This is the portable baseline for MPI stacks that do not accept accelerator pointers.

### direct_device

Packed SYCL USM pointers are passed directly to `MPI_Alltoallv`. This mode requires `assume_gpu_aware_mpi=true`; it is intentionally opt-in because accepting SYCL device USM is a capability of the MPI implementation/runtime and cannot be inferred from the MPI API alone.

Small count exchanges and scalar synchronization remain host-visible. Full particle fields do not need to be downloaded.

## Ghost buffer and overlap seam

Incoming ghost particles remain in a dedicated device buffer instead of being appended to the owned DEM state. The exchange layer also produces a one-byte boundary mask for the current owned particle array. Interior particles can therefore be scheduled independently from boundary particles while halo traffic is in flight in a future asynchronous implementation.

This release does **not** yet claim resident local<->ghost Hertz-Mindlin evaluation or reverse force exchange on the GPU. Those are the next kernels required before the distributed-memory DEM tracker item can close.

## Validation

The v0.16.5 checks include:

- default CPU regression and Python ABI smoke;
- strict C++20 `-Wall -Wextra -Wpedantic -Werror` compilation of the modified resident DEM source;
- strict combined MPI+SYCL compilation of the new transport layer against the project compile harnesses;
- a CPU-executing SYCL emulation that validates stable IDs through device pack -> slab compaction -> device append;
- re-execution of the v0.16.2 resident Hertz/Mindlin regression against the new global-ID history representation;
- a one-rank MPI/SYCL transport harness for the staged transport path.

The complete default project matrix passes **137/137 CTest targets**, including Python ABI **0.16.5**. Real GPU hardware and real multi-rank GPU-aware MPI are not present in this environment, so direct-device runtime qualification remains open.

## Tracker

No capability checkbox is closed in this release. Phase 13 remains **30/37 (81.1%)** and overall completion remains **643/756 (85.1%)**. This is transport/residency infrastructure for the still-open distributed-memory DEM item rather than a claim that the item is fully qualified.
