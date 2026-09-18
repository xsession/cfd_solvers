#include "cfd/multibody/distributed_dem.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
using namespace cfd::multibody;

void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
void require_close(double a,double b,double tol,const char* message){if(std::abs(a-b)>tol)throw std::runtime_error(message);}

DistributedDemParticle particle(std::uint64_t id,int owner,double x,double vx=0.0,double radius=0.1){
    DistributedDemParticle p;
    p.global_id=id;p.owner_rank=owner;p.state.position={x,0,0};p.state.linear_velocity={vx,0,0};
    p.mass=1.0;p.radius=radius;p.inertia_diagonal={0.004,0.004,0.004};
    return p;
}

void test_cross_rank_contact_force_ownership(){
    std::vector<DistributedDemParticle> rank0{particle(10,0,0.95)};
    std::vector<DistributedDemParticle> rank1{particle(20,1,1.05)};
    std::vector<DistributedDemParticle> ghost0{rank1[0]};
    std::vector<DistributedDemParticle> ghost1{rank0[0]};
    std::vector<DistributedDemContactHistory> hist0,hist1;
    HertzMindlinContactModel model;
    model.normal_stiffness=1.0e4;model.normal_damping=0.0;model.tangential_stiffness=100.0;model.tangential_damping=0.0;model.friction=0.5;model.rolling_resistance=0.0;

    auto r0=evaluate_distributed_dem_contacts(rank0,ghost0,hist0,0,1.0e-3,model);
    auto r1=evaluate_distributed_dem_contacts(rank1,ghost1,hist1,1,1.0e-3,model);
    require(r0.contacts==1U,"lower-ID owner must evaluate cross-rank contact");
    require(r1.contacts==0U,"higher-ID owner must not duplicate cross-rank contact");
    require(r0.local.size()==1U && r0.remote.size()==1U,"cross-rank contact must produce local and remote reactions");
    require(r0.remote[0].global_id==20U && r0.remote[0].owner_rank==1,"remote force owner identity");
    require_close(r0.local[0].force.x+r0.remote[0].force.x,0.0,1e-12,"cross-rank contact action/reaction");
    require(hist0.size()==1U && hist0[0].particle_a==10U && hist0[0].particle_b==20U,"contact history canonical owner");

    std::vector<DistributedDemForceContribution> f0=r0.local;
    std::vector<DistributedDemForceContribution> f1;
    accumulate_distributed_dem_forces(f1,r0.remote);
    integrate_distributed_dem_particles(rank0,f0,{0,0,0},1.0e-3);
    integrate_distributed_dem_particles(rank1,f1,{0,0,0},1.0e-3);
    require(rank0[0].state.linear_velocity.x<0.0 && rank1[0].state.linear_velocity.x>0.0,"contact must separate cross-rank particles");
    require_close(rank0[0].state.linear_velocity.x+rank1[0].state.linear_velocity.x,0.0,1e-12,"cross-rank contact momentum conservation");

    ghost0[0]=rank1[0];
    const auto second=evaluate_distributed_dem_contacts(rank0,ghost0,hist0,0,1.0e-3,model);
    require(second.contacts==1U && hist0.size()==1U && hist0[0].age==2U,"persistent cross-rank Mindlin history age");
}

void test_cross_rank_bond_force_and_damage(){
    std::vector<DistributedDemParticle> rank0{particle(1,0,0.90)};
    std::vector<DistributedDemParticle> rank1{particle(2,1,1.10)};
    BondedParticleModel model;
    model.normal_stiffness=1000.0;model.shear_stiffness=500.0;model.normal_damping=0.0;model.shear_damping=0.0;
    model.tensile_failure_force=500.0;model.shear_failure_force=500.0;model.damage_onset_ratio=0.5;
    std::vector<DistributedDemBond> bonds{make_distributed_dem_bond(rank0[0],rank1[0],model)};
    rank1[0].state.position.x=1.15;
    std::vector<DistributedDemParticle> ghosts{rank1[0]};
    auto result=evaluate_distributed_dem_bonds(rank0,ghosts,bonds,0,1.0e-3);
    require(result.active_bonds==1U && result.broken_bonds==0U,"cross-rank bond must remain active");
    require(result.local.size()==1U && result.remote.size()==1U,"cross-rank bond force exchange");
    require(result.local[0].force.x>0.0 && result.remote[0].force.x<0.0,"stretched bond must attract particles");
    require_close(result.local[0].force.x+result.remote[0].force.x,0.0,1e-12,"cross-rank bond action/reaction");
    require(bonds[0].particle_a==1U && bonds[0].particle_b==2U,"cross-rank bond canonical history owner");
}

void test_slab_migration_after_integration(){
    const DemSlabDecomposition decomposition{0.0,2.0,2};
    std::vector<DistributedDemParticle> rank0{particle(7,0,0.99,0.20,0.02)};
    std::vector<DistributedDemForceContribution> none;
    integrate_distributed_dem_particles(rank0,none,{0,0,0},0.10);
    require(rank0[0].state.position.x>1.0,"particle must cross slab during integration");
    const auto plan=plan_dem_slab_exchange(rank0,decomposition,0,0.05);
    require(plan.migrate_indices_by_rank[1].size()==1U && plan.migrate_indices_by_rank[1][0]==0U,"next exchange must migrate integrated particle");
}

void test_ghost_planning_covers_contact_halo(){
    const DemSlabDecomposition decomposition{0.0,2.0,2};
    std::vector<DistributedDemParticle> rank0{particle(3,0,0.91,0.0,0.05)};
    const auto narrow=plan_dem_slab_exchange(rank0,decomposition,0,0.01);
    require(narrow.ghost_indices_by_rank[1].empty(),"narrow halo should not ghost distant particle");
    const auto wide=plan_dem_slab_exchange(rank0,decomposition,0,0.05);
    require(wide.ghost_indices_by_rank[1].size()==1U,"expanded halo must include near-boundary contact candidate");
}

}

int main(){
    try{
        test_cross_rank_contact_force_ownership();
        test_cross_rank_bond_force_and_damage();
        test_slab_migration_after_integration();
        test_ghost_planning_covers_contact_halo();
    }catch(const std::exception& e){
        std::cerr<<"v0.16.4 distributed-contact regression failed: "<<e.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.16.4 distributed-contact regression passed\n";
    return 0;
}
