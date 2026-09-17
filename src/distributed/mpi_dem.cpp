#include "cfd/multibody/distributed_dem.hpp"

#if defined(CFD_HAS_MPI)

namespace cfd::distributed { void check_mpi(int error_code,const char* operation); }

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace cfd::multibody {
namespace {

static_assert(std::is_trivially_copyable_v<DistributedDemParticle>);

[[nodiscard]] int checked_int(std::size_t value,const char* what) {
    if(value>static_cast<std::size_t>(std::numeric_limits<int>::max())) throw std::overflow_error(what);
    return static_cast<int>(value);
}

[[nodiscard]] std::vector<int> make_displacements(const std::vector<int>& counts) {
    std::vector<int> displacements(counts.size(),0);
    long long total=0;
    for(std::size_t i=0;i<counts.size();++i){
        if(counts[i]<0) throw std::runtime_error("negative DEM MPI exchange count");
        if(total>static_cast<long long>(std::numeric_limits<int>::max())) throw std::overflow_error("DEM MPI displacement exceeds int range");
        displacements[i]=static_cast<int>(total);
        total+=counts[i];
    }
    if(total>static_cast<long long>(std::numeric_limits<int>::max())) throw std::overflow_error("DEM MPI payload exceeds int range");
    return displacements;
}

[[nodiscard]] std::size_t total_count(const std::vector<int>& counts) {
    std::size_t total=0;
    for(int count:counts){ if(count<0) throw std::runtime_error("negative DEM MPI exchange count"); total+=static_cast<std::size_t>(count); }
    return total;
}

[[nodiscard]] std::vector<DistributedDemParticle> alltoall_particles(
    MPI_Comm communicator,
    const std::vector<std::vector<DistributedDemParticle>>& send_by_rank) {
    const std::size_t ranks=send_by_rank.size();
    std::vector<int> send_counts(ranks,0);
    for(std::size_t r=0;r<ranks;++r) send_counts[r]=checked_int(send_by_rank[r].size(),"DEM MPI send count exceeds int range");
    std::vector<int> receive_counts(ranks,0);
    cfd::distributed::check_mpi(MPI_Alltoall(send_counts.data(),1,MPI_INT,receive_counts.data(),1,MPI_INT,communicator),"MPI_Alltoall(DEM counts)");
    const auto send_displacements=make_displacements(send_counts);
    const auto receive_displacements=make_displacements(receive_counts);
    std::vector<DistributedDemParticle> send(total_count(send_counts));
    for(std::size_t r=0;r<ranks;++r){
        std::copy(send_by_rank[r].begin(),send_by_rank[r].end(),send.begin()+static_cast<std::ptrdiff_t>(send_displacements[r]));
    }
    std::vector<DistributedDemParticle> receive(total_count(receive_counts));
    std::vector<int> send_bytes(ranks,0),receive_bytes(ranks,0),send_byte_displacements(ranks,0),receive_byte_displacements(ranks,0);
    for(std::size_t r=0;r<ranks;++r){
        send_bytes[r]=checked_int(static_cast<std::size_t>(send_counts[r])*sizeof(DistributedDemParticle),"DEM MPI byte count exceeds int range");
        receive_bytes[r]=checked_int(static_cast<std::size_t>(receive_counts[r])*sizeof(DistributedDemParticle),"DEM MPI byte count exceeds int range");
        send_byte_displacements[r]=checked_int(static_cast<std::size_t>(send_displacements[r])*sizeof(DistributedDemParticle),"DEM MPI byte displacement exceeds int range");
        receive_byte_displacements[r]=checked_int(static_cast<std::size_t>(receive_displacements[r])*sizeof(DistributedDemParticle),"DEM MPI byte displacement exceeds int range");
    }
    cfd::distributed::check_mpi(MPI_Alltoallv(send.data(),send_bytes.data(),send_byte_displacements.data(),MPI_BYTE,
                                               receive.data(),receive_bytes.data(),receive_byte_displacements.data(),MPI_BYTE,
                                               communicator),"MPI_Alltoallv(DEM particles)");
    return receive;
}

} // namespace

MpiDemDomainExchange::MpiDemDomainExchange(DemSlabDecomposition decomposition,MPI_Comm communicator):decomposition_(decomposition) {
    cfd::distributed::check_mpi(MPI_Comm_dup(communicator,&communicator_),"MPI_Comm_dup(DEM)");
    try{
        cfd::distributed::check_mpi(MPI_Comm_rank(communicator_,&rank_),"MPI_Comm_rank(DEM)");
        cfd::distributed::check_mpi(MPI_Comm_size(communicator_,&size_),"MPI_Comm_size(DEM)");
        decomposition_.ranks=size_;
        decomposition_.validate();
    }catch(...){
        if(communicator_!=MPI_COMM_NULL) MPI_Comm_free(&communicator_);
        communicator_=MPI_COMM_NULL;
        throw;
    }
}

MpiDemDomainExchange::~MpiDemDomainExchange() {
    if(communicator_==MPI_COMM_NULL) return;
    int finalized=0;
    MPI_Finalized(&finalized);
    if(!finalized) MPI_Comm_free(&communicator_);
}

std::vector<DistributedDemParticle> MpiDemDomainExchange::exchange(std::vector<DistributedDemParticle>& owned,double ghost_width) {
    auto migration_plan=plan_dem_slab_exchange(owned,decomposition_,rank_,ghost_width);
    std::vector<std::vector<DistributedDemParticle>> migrate(static_cast<std::size_t>(size_));
    std::vector<bool> leaving(owned.size(),false);
    for(int r=0;r<size_;++r){
        for(std::size_t index:migration_plan.migrate_indices_by_rank[static_cast<std::size_t>(r)]){
            auto particle=owned[index];
            particle.owner_rank=r;
            migrate[static_cast<std::size_t>(r)].push_back(particle);
            leaving[index]=true;
        }
    }
    std::vector<DistributedDemParticle> kept;
    kept.reserve(owned.size());
    for(std::size_t i=0;i<owned.size();++i){
        if(!leaving[i]){ owned[i].owner_rank=rank_; kept.push_back(owned[i]); }
    }
    auto incoming=alltoall_particles(communicator_,migrate);
    for(auto& particle:incoming) particle.owner_rank=rank_;
    kept.insert(kept.end(),incoming.begin(),incoming.end());
    owned=std::move(kept);

    auto ghost_plan=plan_dem_slab_exchange(owned,decomposition_,rank_,ghost_width);
    std::vector<std::vector<DistributedDemParticle>> ghost_send(static_cast<std::size_t>(size_));
    for(int r=0;r<size_;++r){
        for(std::size_t index:ghost_plan.ghost_indices_by_rank[static_cast<std::size_t>(r)]) ghost_send[static_cast<std::size_t>(r)].push_back(owned[index]);
    }
    return alltoall_particles(communicator_,ghost_send);
}

} // namespace cfd::multibody
#endif
