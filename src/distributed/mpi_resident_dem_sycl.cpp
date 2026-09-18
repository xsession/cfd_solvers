#include "cfd/multibody/mpi_resident_dem_sycl.hpp"

#if defined(CFD_HAS_MPI) && defined(CFD_HAS_SYCL)

#include "cfd/distributed/mpi_runtime.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace cfd::multibody {
namespace {
using GlobalAtomicU32 = sycl::atomic_ref<std::uint32_t,sycl::memory_order::relaxed,
    sycl::memory_scope::device,sycl::access::address_space::global_space>;

[[nodiscard]] int checked_int(std::size_t value,const char* what) {
    if(value>static_cast<std::size_t>(std::numeric_limits<int>::max())) throw std::overflow_error(what);
    return static_cast<int>(value);
}

template<class T>
[[nodiscard]] std::vector<int> exchange_counts(
    MPI_Comm communicator,
    const std::uint32_t* send_counts,
    int ranks,
    const char* operation) {
    static_assert(std::is_trivially_copyable_v<T>);
    std::vector<int> send(static_cast<std::size_t>(ranks),0),receive(static_cast<std::size_t>(ranks),0);
    for(int r=0;r<ranks;++r) send[static_cast<std::size_t>(r)]=checked_int(send_counts[static_cast<std::size_t>(r)],"resident DEM MPI count exceeds int range");
    cfd::distributed::check_mpi(MPI_Alltoall(send.data(),1,MPI_INT,receive.data(),1,MPI_INT,communicator),operation);
    return receive;
}

[[nodiscard]] std::size_t sum_counts(const std::uint32_t* counts,int ranks) {
    std::size_t total=0U;
    for(int r=0;r<ranks;++r) total+=counts[static_cast<std::size_t>(r)];
    return total;
}
[[nodiscard]] std::size_t sum_counts(const std::vector<int>& counts) {
    std::size_t total=0U;
    for(int count:counts){if(count<0)throw std::runtime_error("negative resident DEM MPI count");total+=static_cast<std::size_t>(count);}
    return total;
}

template<class T>
void exchange_fixed_segments(
    sycl::queue& queue,
    MPI_Comm communicator,
    int ranks,
    std::size_t stride,
    const std::uint32_t* send_counts_u32,
    const std::vector<int>& receive_counts,
    T* device_send,
    T* device_receive,
    ResidentDemMpiMode mode,
    ResidentDemMpiExchangeStats& stats,
    const char* operation) {
    static_assert(std::is_trivially_copyable_v<T>);
    std::vector<int> send_bytes(static_cast<std::size_t>(ranks),0),receive_bytes(static_cast<std::size_t>(ranks),0);
    std::vector<int> send_displacements(static_cast<std::size_t>(ranks),0),receive_displacements(static_cast<std::size_t>(ranks),0);
    std::size_t send_total=0U,receive_total=0U;
    for(int r=0;r<ranks;++r){
        const auto rr=static_cast<std::size_t>(r);
        const auto send_count=static_cast<std::size_t>(send_counts_u32[rr]);
        const auto receive_count=static_cast<std::size_t>(receive_counts[rr]);
        send_bytes[rr]=checked_int(send_count*sizeof(T),"resident DEM MPI send byte count exceeds int range");
        receive_bytes[rr]=checked_int(receive_count*sizeof(T),"resident DEM MPI receive byte count exceeds int range");
        send_displacements[rr]=checked_int(rr*stride*sizeof(T),"resident DEM MPI send displacement exceeds int range");
        receive_displacements[rr]=checked_int(rr*stride*sizeof(T),"resident DEM MPI receive displacement exceeds int range");
        send_total+=send_count*sizeof(T);receive_total+=receive_count*sizeof(T);
    }

    if(mode==ResidentDemMpiMode::direct_device){
        cfd::distributed::check_mpi(
            MPI_Alltoallv(device_send,send_bytes.data(),send_displacements.data(),MPI_BYTE,
                          device_receive,receive_bytes.data(),receive_displacements.data(),MPI_BYTE,
                          communicator),operation);
        stats.direct_device_bytes+=send_total+receive_total;
        return;
    }

    std::vector<T> host_send(static_cast<std::size_t>(ranks)*stride);
    std::vector<T> host_receive(static_cast<std::size_t>(ranks)*stride);
    for(int r=0;r<ranks;++r){
        const auto rr=static_cast<std::size_t>(r),count=static_cast<std::size_t>(send_counts_u32[rr]);
        if(count) queue.memcpy(host_send.data()+rr*stride,device_send+rr*stride,count*sizeof(T));
    }
    queue.wait_and_throw();++stats.scalar_synchronizations;stats.staged_host_bytes+=send_total;
    cfd::distributed::check_mpi(
        MPI_Alltoallv(host_send.data(),send_bytes.data(),send_displacements.data(),MPI_BYTE,
                      host_receive.data(),receive_bytes.data(),receive_displacements.data(),MPI_BYTE,
                      communicator),operation);
    for(int r=0;r<ranks;++r){
        const auto rr=static_cast<std::size_t>(r),count=static_cast<std::size_t>(receive_counts[rr]);
        if(count) queue.memcpy(device_receive+rr*stride,host_receive.data()+rr*stride,count*sizeof(T));
    }
    queue.wait_and_throw();++stats.scalar_synchronizations;stats.staged_host_bytes+=receive_total;
}

template<class T>
struct PendingFixedExchange {
    ResidentDemMpiMode mode{ResidentDemMpiMode::staged_host};
    std::size_t stride{};
    std::vector<int> receive_counts{};
    std::vector<int> send_bytes{},receive_bytes{},send_displacements{},receive_displacements{};
    std::vector<T> host_send{},host_receive{};
    MPI_Request request{MPI_REQUEST_NULL};
    std::size_t receive_total_bytes{};
};

template<class T>
[[nodiscard]] PendingFixedExchange<T> begin_fixed_segments_nonblocking(
    sycl::queue& queue,
    MPI_Comm communicator,
    int ranks,
    std::size_t stride,
    const std::uint32_t* send_counts_u32,
    std::vector<int> receive_counts,
    T* device_send,
    T* device_receive,
    ResidentDemMpiMode mode,
    ResidentDemMpiExchangeStats& stats,
    const char* operation) {
    static_assert(std::is_trivially_copyable_v<T>);
    PendingFixedExchange<T> pending;
    pending.mode=mode;pending.stride=stride;pending.receive_counts=std::move(receive_counts);
    const auto nr=static_cast<std::size_t>(ranks);
    pending.send_bytes.assign(nr,0);pending.receive_bytes.assign(nr,0);pending.send_displacements.assign(nr,0);pending.receive_displacements.assign(nr,0);
    std::size_t send_total=0U,receive_total=0U;
    for(int r=0;r<ranks;++r){
        const auto rr=static_cast<std::size_t>(r),sc=static_cast<std::size_t>(send_counts_u32[rr]),rc=static_cast<std::size_t>(pending.receive_counts[rr]);
        pending.send_bytes[rr]=checked_int(sc*sizeof(T),"resident DEM nonblocking send bytes exceed int range");
        pending.receive_bytes[rr]=checked_int(rc*sizeof(T),"resident DEM nonblocking receive bytes exceed int range");
        pending.send_displacements[rr]=checked_int(rr*stride*sizeof(T),"resident DEM nonblocking send displacement exceeds int range");
        pending.receive_displacements[rr]=checked_int(rr*stride*sizeof(T),"resident DEM nonblocking receive displacement exceeds int range");
        send_total+=sc*sizeof(T);receive_total+=rc*sizeof(T);
    }
    pending.receive_total_bytes=receive_total;
    if(mode==ResidentDemMpiMode::direct_device){
        cfd::distributed::check_mpi(
            MPI_Ialltoallv(device_send,pending.send_bytes.data(),pending.send_displacements.data(),MPI_BYTE,
                           device_receive,pending.receive_bytes.data(),pending.receive_displacements.data(),MPI_BYTE,
                           communicator,&pending.request),operation);
        stats.direct_device_bytes+=send_total+receive_total;
    }else{
        pending.host_send.resize(nr*stride);pending.host_receive.resize(nr*stride);
        for(int r=0;r<ranks;++r){const auto rr=static_cast<std::size_t>(r),count=static_cast<std::size_t>(send_counts_u32[rr]);if(count)queue.memcpy(pending.host_send.data()+rr*stride,device_send+rr*stride,count*sizeof(T));}
        queue.wait_and_throw();++stats.scalar_synchronizations;stats.staged_host_bytes+=send_total;
        cfd::distributed::check_mpi(
            MPI_Ialltoallv(pending.host_send.data(),pending.send_bytes.data(),pending.send_displacements.data(),MPI_BYTE,
                           pending.host_receive.data(),pending.receive_bytes.data(),pending.receive_displacements.data(),MPI_BYTE,
                           communicator,&pending.request),operation);
    }
    ++stats.nonblocking_exchanges;
    return pending;
}

template<class T>
void finish_fixed_segments_nonblocking(
    sycl::queue& queue,
    PendingFixedExchange<T>& pending,
    int ranks,
    T* device_receive,
    ResidentDemMpiExchangeStats& stats,
    const char* operation) {
    cfd::distributed::check_mpi(MPI_Waitall(1,&pending.request,MPI_STATUSES_IGNORE),operation);
    if(pending.mode==ResidentDemMpiMode::staged_host){
        for(int r=0;r<ranks;++r){const auto rr=static_cast<std::size_t>(r),count=static_cast<std::size_t>(pending.receive_counts[rr]);if(count)queue.memcpy(device_receive+rr*pending.stride,pending.host_receive.data()+rr*pending.stride,count*sizeof(T));}
        queue.wait_and_throw();++stats.scalar_synchronizations;stats.staged_host_bytes+=pending.receive_total_bytes;
    }
}

struct PendingCompactForceExchange {
    ResidentDemMpiMode mode{ResidentDemMpiMode::staged_host};
    std::vector<int> receive_counts{},send_bytes{},receive_bytes{},send_displacements{},receive_displacements{};
    std::vector<ResidentDemPackedForce> host_send{},host_receive{};
    MPI_Request request{MPI_REQUEST_NULL};
    std::size_t receive_total{};
    std::size_t receive_total_bytes{};
};

[[nodiscard]] PendingCompactForceExchange begin_compact_force_exchange(
    sycl::queue& queue,MPI_Comm communicator,int ranks,
    const std::uint32_t* send_counts,const std::uint32_t* send_offsets,
    std::vector<int> receive_counts,
    ResidentDemPackedForce* device_send,ResidentDemPackedForce* device_receive,
    ResidentDemMpiMode mode,ResidentDemMpiExchangeStats& stats) {
    PendingCompactForceExchange pending;pending.mode=mode;pending.receive_counts=std::move(receive_counts);
    const auto nr=static_cast<std::size_t>(ranks);pending.send_bytes.assign(nr,0);pending.receive_bytes.assign(nr,0);pending.send_displacements.assign(nr,0);pending.receive_displacements.assign(nr,0);
    std::size_t recv_offset=0U,send_total=0U;
    for(int r=0;r<ranks;++r){const auto rr=static_cast<std::size_t>(r),sc=static_cast<std::size_t>(send_counts[rr]),rc=static_cast<std::size_t>(pending.receive_counts[rr]);pending.send_bytes[rr]=checked_int(sc*sizeof(ResidentDemPackedForce),"reverse DEM force send bytes exceed int range");pending.receive_bytes[rr]=checked_int(rc*sizeof(ResidentDemPackedForce),"reverse DEM force receive bytes exceed int range");pending.send_displacements[rr]=checked_int(static_cast<std::size_t>(send_offsets[rr])*sizeof(ResidentDemPackedForce),"reverse DEM force send displacement exceeds int range");pending.receive_displacements[rr]=checked_int(recv_offset*sizeof(ResidentDemPackedForce),"reverse DEM force receive displacement exceeds int range");send_total+=sc*sizeof(ResidentDemPackedForce);recv_offset+=rc;}
    pending.receive_total=recv_offset;pending.receive_total_bytes=recv_offset*sizeof(ResidentDemPackedForce);
    if(mode==ResidentDemMpiMode::direct_device){
        cfd::distributed::check_mpi(MPI_Ialltoallv(device_send,pending.send_bytes.data(),pending.send_displacements.data(),MPI_BYTE,device_receive,pending.receive_bytes.data(),pending.receive_displacements.data(),MPI_BYTE,communicator,&pending.request),"MPI_Ialltoallv(resident DEM reverse forces)");
        stats.direct_device_bytes+=send_total+pending.receive_total_bytes;
    }else{
        const std::size_t send_records=send_total/sizeof(ResidentDemPackedForce);pending.host_send.resize(send_records);pending.host_receive.resize(recv_offset);
        if(send_records)queue.memcpy(pending.host_send.data(),device_send,send_total);
        queue.wait_and_throw();++stats.scalar_synchronizations;stats.staged_host_bytes+=send_total;
        cfd::distributed::check_mpi(MPI_Ialltoallv(pending.host_send.data(),pending.send_bytes.data(),pending.send_displacements.data(),MPI_BYTE,pending.host_receive.data(),pending.receive_bytes.data(),pending.receive_displacements.data(),MPI_BYTE,communicator,&pending.request),"MPI_Ialltoallv(resident DEM reverse forces)");
    }
    ++stats.nonblocking_exchanges;return pending;
}

void finish_compact_force_exchange(sycl::queue& queue,PendingCompactForceExchange& pending,ResidentDemPackedForce* device_receive,ResidentDemMpiExchangeStats& stats){
    cfd::distributed::check_mpi(MPI_Waitall(1,&pending.request,MPI_STATUSES_IGNORE),"MPI_Waitall(resident DEM reverse forces)");
    if(pending.mode==ResidentDemMpiMode::staged_host){if(pending.receive_total)queue.memcpy(device_receive,pending.host_receive.data(),pending.receive_total_bytes);queue.wait_and_throw();++stats.scalar_synchronizations;stats.staged_host_bytes+=pending.receive_total_bytes;}
}

[[nodiscard]] bool valid_bond_host(const ResidentDemPackedBond& b){
    return b.particle_a!=invalid_dem_global_id && b.particle_a<b.particle_b && b.rest_length>0.0F &&
        b.normal_stiffness>=0.0F && b.shear_stiffness>=0.0F && b.normal_damping>=0.0F && b.shear_damping>=0.0F &&
        b.tensile_failure_force>0.0F && b.shear_failure_force>0.0F && b.damage_onset_ratio>=0.0F && b.damage_onset_ratio<1.0F &&
        b.damage>=0.0F && b.damage<=1.0F && std::isfinite(b.rest_length) && std::isfinite(b.damage);
}
}

MpiResidentDemSyclExchange::MpiResidentDemSyclExchange(
    sycl::queue& communication_queue,
    ResidentDemMpiExchangeConfig config,
    MPI_Comm communicator)
    : queue_(communication_queue),config_(config) {
    if(config_.mode==ResidentDemMpiMode::direct_device && !config_.assume_gpu_aware_mpi) throw std::invalid_argument("direct-device DEM MPI requires explicit GPU-aware MPI opt-in");
    cfd::distributed::check_mpi(MPI_Comm_dup(communicator,&communicator_),"MPI_Comm_dup(resident DEM)");
    try {
        cfd::distributed::check_mpi(MPI_Comm_rank(communicator_,&rank_),"MPI_Comm_rank(resident DEM)");
        cfd::distributed::check_mpi(MPI_Comm_size(communicator_,&size_),"MPI_Comm_size(resident DEM)");
        config_.decomposition.ranks=size_;config_.decomposition.validate();
        if(config_.bond_capacity) reserve_bonds(config_.bond_capacity);
    } catch(...) {
        if(communicator_!=MPI_COMM_NULL) MPI_Comm_free(&communicator_);
        communicator_=MPI_COMM_NULL;
        throw;
    }
}

MpiResidentDemSyclExchange::~MpiResidentDemSyclExchange(){ release(); }

void MpiResidentDemSyclExchange::release() noexcept {
    try{queue_.wait_and_throw();}catch(...){}
    auto freep=[&](auto*& p){if(p){sycl::free(p,queue_);p=nullptr;}};
    freep(migrate_counts_);freep(ghost_counts_);freep(receive_counts_);freep(migrate_indices_);freep(ghost_indices_);freep(boundary_mask_);
    freep(remote_force_count_);freep(reverse_counts_);freep(reverse_cursor_);freep(reverse_offsets_);
    freep(migrate_send_);freep(migrate_receive_);freep(migrate_history_send_);freep(migrate_history_receive_);freep(ghost_send_);freep(ghost_receive_);
    freep(remote_force_raw_);freep(reverse_force_send_);freep(reverse_force_receive_);
    freep(bonds_);freep(bond_send_);freep(bond_receive_);freep(bond_counts_);freep(bond_receive_counts_);freep(bond_compaction_count_);
    reserved_capacity_=0U;reserved_slots_=0U;ghost_count_=0U;boundary_mask_count_=0U;remote_force_capacity_=0U;bond_capacity_=0U;bond_count_=0U;
    if(communicator_!=MPI_COMM_NULL){int finalized=0;MPI_Finalized(&finalized);if(!finalized)MPI_Comm_free(&communicator_);communicator_=MPI_COMM_NULL;}
}

void MpiResidentDemSyclExchange::reserve_bonds(std::size_t capacity){
    if(capacity==0U) capacity=1U;
    if(capacity<=bond_capacity_ && bonds_ && bond_send_ && bond_receive_) return;
    const auto old_count=bond_count_;std::vector<ResidentDemPackedBond> old;
    if(bonds_&&old_count){old.resize(old_count);queue_.memcpy(old.data(),bonds_,old_count*sizeof(ResidentDemPackedBond)).wait_and_throw();}
    auto freep=[&](auto*& p){if(p){sycl::free(p,queue_);p=nullptr;}};freep(bonds_);freep(bond_send_);freep(bond_receive_);freep(bond_counts_);freep(bond_receive_counts_);freep(bond_compaction_count_);
    bond_capacity_=capacity;const auto ranks=static_cast<std::size_t>(size_);
    bonds_=sycl::malloc_device<ResidentDemPackedBond>(bond_capacity_,queue_);bond_send_=sycl::malloc_device<ResidentDemPackedBond>(ranks*bond_capacity_,queue_);bond_receive_=sycl::malloc_device<ResidentDemPackedBond>(ranks*bond_capacity_,queue_);
    bond_counts_=sycl::malloc_shared<std::uint32_t>(ranks,queue_);bond_receive_counts_=sycl::malloc_shared<std::uint32_t>(ranks,queue_);bond_compaction_count_=sycl::malloc_shared<std::uint32_t>(1U,queue_);
    if(!bonds_||!bond_send_||!bond_receive_||!bond_counts_||!bond_receive_counts_||!bond_compaction_count_)throw std::bad_alloc{};
    if(!old.empty()) queue_.memcpy(bonds_,old.data(),old.size()*sizeof(ResidentDemPackedBond)).wait_and_throw();
    bond_count_=old_count;
}

void MpiResidentDemSyclExchange::reserve(const ResidentDemSycl& dem) {
    if(reserved_capacity_==dem.capacity() && reserved_slots_==dem.history_slots_per_particle()) return;
    if(reserved_capacity_!=0U) {
        auto freep=[&](auto*& p){if(p){sycl::free(p,queue_);p=nullptr;}};
        freep(migrate_counts_);freep(ghost_counts_);freep(receive_counts_);freep(migrate_indices_);freep(ghost_indices_);freep(boundary_mask_);freep(remote_force_count_);freep(reverse_counts_);freep(reverse_cursor_);freep(reverse_offsets_);
        freep(migrate_send_);freep(migrate_receive_);freep(migrate_history_send_);freep(migrate_history_receive_);freep(ghost_send_);freep(ghost_receive_);freep(remote_force_raw_);freep(reverse_force_send_);freep(reverse_force_receive_);
    }
    reserved_capacity_=dem.capacity();reserved_slots_=dem.history_slots_per_particle();
    if(bond_capacity_==0U)reserve_bonds(config_.bond_capacity?config_.bond_capacity:std::max<std::size_t>(1U,4U*reserved_capacity_));
    remote_force_capacity_=reserved_capacity_*reserved_slots_+bond_capacity_;if(remote_force_capacity_<reserved_capacity_)throw std::overflow_error("resident DEM remote force capacity overflow");
    const auto ranks=static_cast<std::size_t>(size_),segments=ranks*reserved_capacity_,history_segments=segments*reserved_slots_;
    migrate_counts_=sycl::malloc_shared<std::uint32_t>(ranks,queue_);ghost_counts_=sycl::malloc_shared<std::uint32_t>(ranks,queue_);receive_counts_=sycl::malloc_shared<std::uint32_t>(ranks,queue_);
    reverse_counts_=sycl::malloc_shared<std::uint32_t>(ranks,queue_);reverse_cursor_=sycl::malloc_shared<std::uint32_t>(ranks,queue_);reverse_offsets_=sycl::malloc_shared<std::uint32_t>(ranks+1U,queue_);remote_force_count_=sycl::malloc_shared<std::uint32_t>(1U,queue_);
    migrate_indices_=sycl::malloc_device<std::uint32_t>(segments,queue_);ghost_indices_=sycl::malloc_device<std::uint32_t>(segments,queue_);boundary_mask_=sycl::malloc_device<std::uint8_t>(reserved_capacity_,queue_);
    migrate_send_=sycl::malloc_device<ResidentDemPackedParticle>(segments,queue_);migrate_receive_=sycl::malloc_device<ResidentDemPackedParticle>(segments,queue_);
    migrate_history_send_=sycl::malloc_device<ResidentDemPackedHistory>(history_segments,queue_);migrate_history_receive_=sycl::malloc_device<ResidentDemPackedHistory>(history_segments,queue_);
    ghost_send_=sycl::malloc_device<ResidentDemPackedParticle>(segments,queue_);ghost_receive_=sycl::malloc_device<ResidentDemPackedParticle>(segments,queue_);
    remote_force_raw_=sycl::malloc_device<ResidentDemPackedForce>(remote_force_capacity_,queue_);reverse_force_send_=sycl::malloc_device<ResidentDemPackedForce>(remote_force_capacity_,queue_);reverse_force_receive_=sycl::malloc_device<ResidentDemPackedForce>(remote_force_capacity_,queue_);
    if(!migrate_counts_||!ghost_counts_||!receive_counts_||!reverse_counts_||!reverse_cursor_||!reverse_offsets_||!remote_force_count_||!migrate_indices_||!ghost_indices_||!boundary_mask_||!migrate_send_||!migrate_receive_||!migrate_history_send_||!migrate_history_receive_||!ghost_send_||!ghost_receive_||!remote_force_raw_||!reverse_force_send_||!reverse_force_receive_)throw std::bad_alloc{};
}

void MpiResidentDemSyclExchange::upload_bonds(std::span<const ResidentDemPackedBond> bonds){
    for(const auto& bond:bonds)if(!valid_bond_host(bond))throw std::invalid_argument("invalid resident distributed DEM bond");
    if(bonds.size()>bond_capacity_)reserve_bonds(std::max<std::size_t>(bonds.size(),config_.bond_capacity));
    if(!bonds_&&bonds.empty()){bond_count_=0U;return;}
    if(!bonds.empty()) queue_.memcpy(bonds_,bonds.data(),bonds.size_bytes()).wait_and_throw();
    bond_count_=bonds.size();
}

std::vector<ResidentDemPackedBond> MpiResidentDemSyclExchange::download_bonds() const{
    std::vector<ResidentDemPackedBond> out(bond_count_);if(bond_count_)queue_.memcpy(out.data(),bonds_,bond_count_*sizeof(ResidentDemPackedBond)).wait_and_throw();return out;
}

void MpiResidentDemSyclExchange::clear_counts(){for(int r=0;r<size_;++r){const auto rr=static_cast<std::size_t>(r);migrate_counts_[rr]=0U;ghost_counts_[rr]=0U;receive_counts_[rr]=0U;reverse_counts_[rr]=0U;reverse_cursor_[rr]=0U;}if(remote_force_count_)*remote_force_count_=0U;}

void MpiResidentDemSyclExchange::classify_migration(const ResidentDemSycl& dem) {
    for(int r=0;r<size_;++r)migrate_counts_[static_cast<std::size_t>(r)]=0U;
    const auto view=const_cast<ResidentDemSycl&>(dem).device_view();const auto n=view.count,cap=dem.capacity();const double xmin=config_.decomposition.x_min,width=config_.decomposition.slab_width();const int ranks=size_,local_rank=rank_;auto* counts=migrate_counts_;auto* indices=migrate_indices_;const float* x=view.x;
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const std::size_t i=id[0];long long destination=static_cast<long long>((static_cast<double>(x[i])-xmin)/width);if(destination<0)destination=0;if(destination>=ranks)destination=ranks-1;if(destination==local_rank)return;const auto d=static_cast<std::size_t>(destination);const auto slot=GlobalAtomicU32(counts[d]).fetch_add(1U);if(slot<cap)indices[d*cap+slot]=static_cast<std::uint32_t>(i);});
    queue_.wait_and_throw();++stats_.scalar_synchronizations;
}

void MpiResidentDemSyclExchange::classify_ghosts(const ResidentDemSycl& dem,double ghost_width) {
    if(!(ghost_width>=0.0))throw std::invalid_argument("resident DEM ghost width must be non-negative");
    for(int r=0;r<size_;++r)ghost_counts_[static_cast<std::size_t>(r)]=0U;
    const auto view=const_cast<ResidentDemSycl&>(dem).device_view();const auto n=view.count,cap=dem.capacity();const float lo=static_cast<float>(config_.decomposition.slab_min(rank_)),hi=static_cast<float>(config_.decomposition.slab_max(rank_)),halo=static_cast<float>(ghost_width);const int local_rank=rank_,ranks=size_;auto* counts=ghost_counts_;auto* indices=ghost_indices_;auto* boundary=boundary_mask_;const float* x=view.x;const float* radius=view.radius;
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const std::size_t i=id[0];const bool left=local_rank>0&&x[i]-radius[i]-halo<lo,right=local_rank+1<ranks&&x[i]+radius[i]+halo>=hi;boundary[i]=static_cast<std::uint8_t>(left||right);if(left){const auto r=static_cast<std::size_t>(local_rank-1);const auto slot=GlobalAtomicU32(counts[r]).fetch_add(1U);if(slot<cap)indices[r*cap+slot]=static_cast<std::uint32_t>(i);}if(right){const auto r=static_cast<std::size_t>(local_rank+1);const auto slot=GlobalAtomicU32(counts[r]).fetch_add(1U);if(slot<cap)indices[r*cap+slot]=static_cast<std::uint32_t>(i);}});
    queue_.wait_and_throw();++stats_.scalar_synchronizations;boundary_mask_count_=n;
}

void MpiResidentDemSyclExchange::pack_particle_segments(ResidentDemSycl& dem,const std::uint32_t* indices,const std::uint32_t* counts,ResidentDemPackedParticle* particles,ResidentDemPackedHistory* history){const auto cap=dem.capacity(),slots=dem.history_slots_per_particle();for(int r=0;r<size_;++r){const auto rr=static_cast<std::size_t>(r),count=static_cast<std::size_t>(counts[rr]);if(count)dem.pack_indices_device(indices+rr*cap,count,particles+rr*cap,history+rr*cap*slots);}}

void MpiResidentDemSyclExchange::exchange_bonds_for_migration(ResidentDemSycl& dem){
    if(bond_count_==0U) return;
    for(int r=0;r<size_;++r){bond_counts_[static_cast<std::size_t>(r)]=0U;bond_receive_counts_[static_cast<std::size_t>(r)]=0U;}
    const auto view=dem.exchange_device_view();const auto n=view.count;const auto* gid=view.global_id;const auto* x=view.x;auto* bonds=bonds_;auto* send=bond_send_;auto* counts=bond_counts_;const auto bcount=bond_count_,bcap=bond_capacity_;const double xmin=config_.decomposition.x_min,width=config_.decomposition.slab_width();const int ranks=size_,local_rank=rank_;
    queue_.parallel_for(sycl::range<1>(bcount),[=](sycl::id<1> id){const auto bond=bonds[id[0]];std::size_t ia=n;for(std::size_t i=0;i<n;++i)if(gid[i]==bond.particle_a){ia=i;break;}if(ia==n)return;long long dest=static_cast<long long>((static_cast<double>(x[ia])-xmin)/width);if(dest<0)dest=0;if(dest>=ranks)dest=ranks-1;if(dest==local_rank)return;const auto d=static_cast<std::size_t>(dest);const auto slot=GlobalAtomicU32(counts[d]).fetch_add(1U);if(slot<bcap)send[d*bcap+slot]=bond;});
    queue_.wait_and_throw();++stats_.scalar_synchronizations;for(int r=0;r<size_;++r)if(bond_counts_[static_cast<std::size_t>(r)]>bond_capacity_)throw std::overflow_error("resident DEM outgoing bond buffer overflow");
    auto receive_counts=exchange_counts<ResidentDemPackedBond>(communicator_,bond_counts_,size_,"MPI_Alltoall(resident DEM bond counts)");if(sum_counts(receive_counts)>bond_capacity_)throw std::overflow_error("resident DEM incoming bond buffer overflow");
    exchange_fixed_segments(queue_,communicator_,size_,bond_capacity_,bond_counts_,receive_counts,bond_send_,bond_receive_,config_.mode,stats_,"MPI_Alltoallv(resident DEM bonds)");
    *bond_compaction_count_=0U;auto* out_count=bond_compaction_count_;const float lo=static_cast<float>(config_.decomposition.slab_min(rank_)),hi=static_cast<float>(config_.decomposition.slab_max(rank_));
    queue_.parallel_for(sycl::range<1>(1U),[=](sycl::id<1>){std::uint32_t dst=0U;for(std::size_t b=0;b<bcount;++b){std::size_t ia=n;for(std::size_t i=0;i<n;++i)if(gid[i]==bonds[b].particle_a){ia=i;break;}if(ia!=n&&x[ia]>=lo&&x[ia]<hi){if(dst!=b)bonds[dst]=bonds[b];++dst;}}*out_count=dst;}).wait_and_throw();bond_count_=static_cast<std::size_t>(*bond_compaction_count_);++stats_.scalar_synchronizations;
    for(int source=0;source<size_;++source){const auto rr=static_cast<std::size_t>(source),count=static_cast<std::size_t>(receive_counts[rr]);if(!count)continue;if(bond_count_+count>bond_capacity_)throw std::overflow_error("resident DEM bond append exceeds capacity");queue_.memcpy(bonds_+bond_count_,bond_receive_+rr*bond_capacity_,count*sizeof(ResidentDemPackedBond));bond_count_+=count;}
    queue_.wait_and_throw();++stats_.scalar_synchronizations;stats_.resident_bonds=bond_count_;
}

void MpiResidentDemSyclExchange::exchange_migrants(ResidentDemSycl& dem) {
    classify_migration(dem);stats_.migrated_out=sum_counts(migrate_counts_,size_);pack_particle_segments(dem,migrate_indices_,migrate_counts_,migrate_send_,migrate_history_send_);queue_.wait_and_throw();++stats_.scalar_synchronizations;
    auto receive_counts=exchange_counts<ResidentDemPackedParticle>(communicator_,migrate_counts_,size_,"MPI_Alltoall(resident DEM migrant counts)");exchange_fixed_segments(queue_,communicator_,size_,reserved_capacity_,migrate_counts_,receive_counts,migrate_send_,migrate_receive_,config_.mode,stats_,"MPI_Alltoallv(resident DEM migrants)");
    std::vector<int> history_receive_counts(receive_counts.size(),0);std::vector<std::uint32_t> history_send_counts(static_cast<std::size_t>(size_),0U);for(int r=0;r<size_;++r){const auto rr=static_cast<std::size_t>(r);history_send_counts[rr]=migrate_counts_[rr]*static_cast<std::uint32_t>(reserved_slots_);history_receive_counts[rr]=receive_counts[rr]*checked_int(reserved_slots_,"resident DEM history slot count exceeds int range");}
    exchange_fixed_segments(queue_,communicator_,size_,reserved_capacity_*reserved_slots_,history_send_counts.data(),history_receive_counts,migrate_history_send_,migrate_history_receive_,config_.mode,stats_,"MPI_Alltoallv(resident DEM histories)");
    exchange_bonds_for_migration(dem);
    dem.retain_x_slab_device(config_.decomposition.slab_min(rank_),config_.decomposition.slab_max(rank_),rank_);stats_.migrated_in=sum_counts(receive_counts);for(int source=0;source<size_;++source){const auto rr=static_cast<std::size_t>(source),count=static_cast<std::size_t>(receive_counts[rr]);if(count)dem.append_packed_device(migrate_receive_+rr*reserved_capacity_,migrate_history_receive_+rr*reserved_capacity_*reserved_slots_,count,rank_);}queue_.wait_and_throw();++stats_.scalar_synchronizations;
}

void MpiResidentDemSyclExchange::exchange_ghosts(ResidentDemSycl& dem,double ghost_width) {
    classify_ghosts(dem,ghost_width);const auto view=dem.exchange_device_view();const auto cap=reserved_capacity_;
    for(int r=0;r<size_;++r){const auto rr=static_cast<std::size_t>(r),count=static_cast<std::size_t>(ghost_counts_[rr]);if(!count)continue;const auto* indices=ghost_indices_+rr*cap;auto* output=ghost_send_+rr*cap;const auto* gid=view.global_id;const auto* owner=view.owner_rank;const auto *x=view.x,*y=view.y,*z=view.z,*vx=view.vx,*vy=view.vy,*vz=view.vz,*wx=view.wx,*wy=view.wy,*wz=view.wz,*radius=view.radius,*im=view.inverse_mass,*ii=view.inverse_inertia;queue_.parallel_for(sycl::range<1>(count),[=](sycl::id<1> id){const std::size_t p=id[0],i=indices[p];auto& out=output[p];out.global_id=gid[i];out.owner_rank=owner[i];out.x=x[i];out.y=y[i];out.z=z[i];out.vx=vx[i];out.vy=vy[i];out.vz=vz[i];out.wx=wx[i];out.wy=wy[i];out.wz=wz[i];out.radius=radius[i];out.inverse_mass=im[i];out.inverse_inertia=ii[i];});}
    queue_.wait_and_throw();++stats_.scalar_synchronizations;stats_.ghosts_out=sum_counts(ghost_counts_,size_);auto receive_counts=exchange_counts<ResidentDemPackedParticle>(communicator_,ghost_counts_,size_,"MPI_Alltoall(resident DEM ghost counts)");exchange_fixed_segments(queue_,communicator_,size_,reserved_capacity_,ghost_counts_,receive_counts,ghost_send_,ghost_receive_,config_.mode,stats_,"MPI_Alltoallv(resident DEM ghosts)");ghost_count_=sum_counts(receive_counts);stats_.ghosts_in=ghost_count_;for(int r=0;r<size_;++r)receive_counts_[static_cast<std::size_t>(r)]=static_cast<std::uint32_t>(receive_counts[static_cast<std::size_t>(r)]);auto* ghosts=ghost_receive_;auto* counts=receive_counts_;const auto stride=reserved_capacity_;const int ranks=size_;queue_.parallel_for(sycl::range<1>(1U),[=](sycl::id<1>){std::size_t dst=0U;for(int r=0;r<ranks;++r){const auto rr=static_cast<std::size_t>(r);for(std::uint32_t k=0;k<counts[rr];++k){const auto src=rr*stride+static_cast<std::size_t>(k);if(dst!=src)ghosts[dst]=ghosts[src];++dst;}}});queue_.wait_and_throw();++stats_.scalar_synchronizations;
}

void MpiResidentDemSyclExchange::exchange_reverse_forces(ResidentDemSycl& dem){
    const auto count=static_cast<std::size_t>(*remote_force_count_);stats_.reverse_forces_out=count;if(count==0U){stats_.reverse_forces_in=0U;return;}if(count>remote_force_capacity_)throw std::overflow_error("resident DEM remote force buffer overflow");for(int r=0;r<size_;++r){reverse_counts_[static_cast<std::size_t>(r)]=0U;reverse_cursor_[static_cast<std::size_t>(r)]=0U;}
    auto* raw=remote_force_raw_;auto* counts=reverse_counts_;const int ranks=size_;queue_.parallel_for(sycl::range<1>(count),[=](sycl::id<1> id){const auto f=raw[id[0]];if(f.owner_rank<0||f.owner_rank>=ranks)return;GlobalAtomicU32(counts[static_cast<std::size_t>(f.owner_rank)]).fetch_add(1U);});queue_.wait_and_throw();++stats_.scalar_synchronizations;
    std::uint32_t offset=0U;for(int r=0;r<size_;++r){const auto rr=static_cast<std::size_t>(r);reverse_offsets_[rr]=offset;offset+=reverse_counts_[rr];reverse_cursor_[rr]=0U;}reverse_offsets_[static_cast<std::size_t>(size_)]=offset;if(static_cast<std::size_t>(offset)>remote_force_capacity_)throw std::overflow_error("resident DEM compact reverse force overflow");
    auto* cursor=reverse_cursor_;auto* offsets=reverse_offsets_;auto* packed=reverse_force_send_;queue_.parallel_for(sycl::range<1>(count),[=](sycl::id<1> id){const auto f=raw[id[0]];if(f.owner_rank<0||f.owner_rank>=ranks)return;const auto rr=static_cast<std::size_t>(f.owner_rank);const auto slot=GlobalAtomicU32(cursor[rr]).fetch_add(1U);packed[static_cast<std::size_t>(offsets[rr])+slot]=f;});queue_.wait_and_throw();++stats_.scalar_synchronizations;
    auto receive_counts=exchange_counts<ResidentDemPackedForce>(communicator_,reverse_counts_,size_,"MPI_Alltoall(resident DEM reverse force counts)");const auto receive_total=sum_counts(receive_counts);if(receive_total>remote_force_capacity_)throw std::overflow_error("resident DEM reverse force receive overflow");auto pending=begin_compact_force_exchange(queue_,communicator_,size_,reverse_counts_,reverse_offsets_,std::move(receive_counts),reverse_force_send_,reverse_force_receive_,config_.mode,stats_);finish_compact_force_exchange(queue_,pending,reverse_force_receive_,stats_);stats_.reverse_forces_in=pending.receive_total;dem.apply_remote_forces_device(reverse_force_receive_,pending.receive_total);
}

ResidentDemMpiExchangeStats MpiResidentDemSyclExchange::exchange(ResidentDemSycl& dem,double ghost_width) {reserve(dem);clear_counts();stats_={};exchange_migrants(dem);exchange_ghosts(dem,ghost_width);stats_.resident_bonds=bond_count_;return stats_;}

ResidentDemMpiExchangeStats MpiResidentDemSyclExchange::step(ResidentDemSycl& dem,double ghost_width){
    if(&dem.queue()!=&queue_) throw std::invalid_argument("resident distributed DEM step requires the solver and exchange to share one in-order SYCL queue");
    reserve(dem);clear_counts();stats_={};exchange_migrants(dem);
    classify_ghosts(dem,ghost_width);const auto view=dem.exchange_device_view();const auto cap=reserved_capacity_;
    for(int r=0;r<size_;++r){const auto rr=static_cast<std::size_t>(r),count=static_cast<std::size_t>(ghost_counts_[rr]);if(!count)continue;const auto* indices=ghost_indices_+rr*cap;auto* output=ghost_send_+rr*cap;const auto* gid=view.global_id;const auto* owner=view.owner_rank;const auto *x=view.x,*y=view.y,*z=view.z,*vx=view.vx,*vy=view.vy,*vz=view.vz,*wx=view.wx,*wy=view.wy,*wz=view.wz,*radius=view.radius,*im=view.inverse_mass,*ii=view.inverse_inertia;queue_.parallel_for(sycl::range<1>(count),[=](sycl::id<1> id){const std::size_t p=id[0],i=indices[p];auto& out=output[p];out.global_id=gid[i];out.owner_rank=owner[i];out.x=x[i];out.y=y[i];out.z=z[i];out.vx=vx[i];out.vy=vy[i];out.vz=vz[i];out.wx=wx[i];out.wy=wy[i];out.wz=wz[i];out.radius=radius[i];out.inverse_mass=im[i];out.inverse_inertia=ii[i];});}
    queue_.wait_and_throw();++stats_.scalar_synchronizations;stats_.ghosts_out=sum_counts(ghost_counts_,size_);auto receive_counts=exchange_counts<ResidentDemPackedParticle>(communicator_,ghost_counts_,size_,"MPI_Alltoall(resident DEM overlap ghost counts)");ghost_count_=sum_counts(receive_counts);stats_.ghosts_in=ghost_count_;auto pending=begin_fixed_segments_nonblocking(queue_,communicator_,size_,reserved_capacity_,ghost_counts_,receive_counts,ghost_send_,ghost_receive_,config_.mode,stats_,"MPI_Ialltoallv(resident DEM ghosts)");
    dem.begin_distributed_contact_step();dem.apply_owned_contacts_interior_device(boundary_mask_);
    finish_fixed_segments_nonblocking(queue_,pending,size_,ghost_receive_,stats_,"MPI_Waitall(resident DEM ghosts)");for(int r=0;r<size_;++r)receive_counts_[static_cast<std::size_t>(r)]=static_cast<std::uint32_t>(pending.receive_counts[static_cast<std::size_t>(r)]);auto* ghosts=ghost_receive_;auto* counts=receive_counts_;const auto stride=reserved_capacity_;const int ranks=size_;queue_.parallel_for(sycl::range<1>(1U),[=](sycl::id<1>){std::size_t dst=0U;for(int r=0;r<ranks;++r){const auto rr=static_cast<std::size_t>(r);for(std::uint32_t k=0;k<counts[rr];++k){const auto src=rr*stride+static_cast<std::size_t>(k);if(dst!=src)ghosts[dst]=ghosts[src];++dst;}}});
    *remote_force_count_=0U;dem.apply_owned_contacts_boundary_device(boundary_mask_);dem.apply_ghost_contacts_device(ghost_receive_,ghost_count_,boundary_mask_,rank_,remote_force_raw_,remote_force_count_,remote_force_capacity_);dem.apply_bonds_device(bonds_,bond_count_,ghost_receive_,ghost_count_,rank_,remote_force_raw_,remote_force_count_,remote_force_capacity_);queue_.wait_and_throw();++stats_.scalar_synchronizations;if(static_cast<std::size_t>(*remote_force_count_)>remote_force_capacity_)throw std::overflow_error("resident DEM contact/bond reverse-force overflow");exchange_reverse_forces(dem);dem.finish_distributed_contact_step();queue_.wait_and_throw();++stats_.scalar_synchronizations;stats_.resident_bonds=bond_count_;return stats_;
}

} // namespace cfd::multibody
#endif
