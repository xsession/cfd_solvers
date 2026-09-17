#include "cfd/multibody/distributed_dem.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::multibody {

void DemSlabDecomposition::validate() const {
    if(!(x_max>x_min) || !std::isfinite(x_min) || !std::isfinite(x_max)) throw std::invalid_argument("DEM slab bounds must be finite and ordered");
    if(ranks<=0) throw std::invalid_argument("DEM slab decomposition requires at least one rank");
}

double DemSlabDecomposition::slab_width() const {
    validate();
    return (x_max-x_min)/static_cast<double>(ranks);
}

int DemSlabDecomposition::owner_rank(double x) const {
    validate();
    if(!std::isfinite(x)) throw std::invalid_argument("DEM particle coordinate must be finite");
    if(x<=x_min) return 0;
    if(x>=x_max) return ranks-1;
    const double relative=(x-x_min)/(x_max-x_min);
    const auto rank=static_cast<int>(std::floor(relative*static_cast<double>(ranks)));
    return std::clamp(rank,0,ranks-1);
}

double DemSlabDecomposition::slab_min(int rank) const {
    validate();
    if(rank<0 || rank>=ranks) throw std::out_of_range("DEM slab rank out of range");
    return x_min+static_cast<double>(rank)*slab_width();
}

double DemSlabDecomposition::slab_max(int rank) const {
    validate();
    if(rank<0 || rank>=ranks) throw std::out_of_range("DEM slab rank out of range");
    return rank+1==ranks ? x_max : x_min+static_cast<double>(rank+1)*slab_width();
}

DistributedDemParticle make_distributed_dem_particle(std::uint64_t global_id,const RigidBody& body,const SphereShape& sphere,int owner_rank) {
    if(sphere.body==invalid_body) throw std::invalid_argument("distributed DEM sphere must reference a body");
    if(!(sphere.radius>0.0) || !std::isfinite(sphere.radius)) throw std::invalid_argument("distributed DEM sphere radius must be finite and positive");
    body.validate();
    DistributedDemParticle particle;
    particle.global_id=global_id;
    particle.state=body.state;
    particle.mass=body.mass;
    particle.inertia_diagonal=body.inertia_diagonal;
    particle.radius=sphere.radius;
    particle.owner_rank=owner_rank;
    return particle;
}

DemExchangePlan plan_dem_slab_exchange(std::span<const DistributedDemParticle> particles,const DemSlabDecomposition& decomposition,int local_rank,double ghost_width) {
    decomposition.validate();
    if(local_rank<0 || local_rank>=decomposition.ranks) throw std::out_of_range("local DEM rank out of range");
    if(!(ghost_width>=0.0) || !std::isfinite(ghost_width) || ghost_width>decomposition.slab_width()) throw std::invalid_argument("DEM ghost width must be finite and no larger than one slab");
    DemExchangePlan plan;
    plan.migrate_indices_by_rank.resize(static_cast<std::size_t>(decomposition.ranks));
    plan.ghost_indices_by_rank.resize(static_cast<std::size_t>(decomposition.ranks));
    const double lower=decomposition.slab_min(local_rank);
    const double upper=decomposition.slab_max(local_rank);
    for(std::size_t i=0;i<particles.size();++i){
        const auto& particle=particles[i];
        if(!(particle.radius>0.0) || !std::isfinite(particle.radius)) throw std::invalid_argument("distributed DEM particle radius must be finite and positive");
        const int owner=decomposition.owner_rank(particle.state.position.x);
        if(owner!=local_rank){
            plan.migrate_indices_by_rank[static_cast<std::size_t>(owner)].push_back(i);
            continue;
        }
        const double reach=particle.radius+ghost_width;
        if(local_rank>0 && particle.state.position.x-reach<lower){
            plan.ghost_indices_by_rank[static_cast<std::size_t>(local_rank-1)].push_back(i);
        }
        if(local_rank+1<decomposition.ranks && particle.state.position.x+reach>upper){
            plan.ghost_indices_by_rank[static_cast<std::size_t>(local_rank+1)].push_back(i);
        }
    }
    return plan;
}

} // namespace cfd::multibody
