#include "cfd/multibody/distributed_dem.hpp"
#include "cfd/multibody/system.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
using namespace cfd::multibody;

void require(bool condition,const char* message){ if(!condition) throw std::runtime_error(message); }
void require_close(double a,double b,double tol,const char* message){ if(std::abs(a-b)>tol) throw std::runtime_error(message); }

RigidBody unit_body(Vec3 position,Vec3 velocity={}){
    RigidBody body;
    body.state.position=position;
    body.state.linear_velocity=velocity;
    body.mass=1.0;
    body.inertia_diagonal={1.0,1.0,1.0};
    return body;
}

void test_bond_elastic_damage_and_failure(){
    std::vector<RigidBody> bodies{unit_body({0,0,0}),unit_body({1.0,0,0})};
    BondedParticleModel model;
    model.normal_stiffness=100.0;
    model.shear_stiffness=80.0;
    model.normal_damping=0.0;
    model.shear_damping=0.0;
    model.tensile_failure_force=50.0;
    model.shear_failure_force=40.0;
    model.damage_onset_ratio=0.5;
    std::vector<ParticleBond> bonds{make_particle_bond(bodies,0,1,model)};

    bodies[1].state.position.x=1.2;
    auto stats=apply_particle_bonds(bodies,bonds,0.01);
    require(stats.active==1 && stats.damaged==0 && stats.broken_total==0,"elastic bond stats");
    require_close(bodies[0].force.x,20.0,1e-12,"elastic bond restoring force on A");
    require_close(bodies[1].force.x,-20.0,1e-12,"elastic bond action/reaction");

    bodies[0].clear_accumulators();bodies[1].clear_accumulators();
    bodies[1].state.position.x=1.4;
    stats=apply_particle_bonds(bodies,bonds,0.01);
    require(stats.active==1 && stats.damaged==1 && !bonds[0].broken,"bond progressive damage state");
    require(bonds[0].damage>0.0 && bonds[0].damage<1.0,"bond damage must be progressive before failure");
    require(std::abs(bodies[0].force.x)<40.0,"damage must reduce effective bond stiffness");

    bodies[0].clear_accumulators();bodies[1].clear_accumulators();
    bodies[1].state.position.x=1.6;
    stats=apply_particle_bonds(bodies,bonds,0.01);
    require(bonds[0].broken && bonds[0].damage==1.0,"tensile bond failure");
    require(stats.broken_this_step==1 && stats.broken_total==1 && stats.active==0,"bond failure statistics");
    require_close(norm(bodies[0].force),0.0,1e-14,"broken bond must not apply force");
}

void test_bond_shear_history(){
    std::vector<RigidBody> bodies{unit_body({0,0,0}),unit_body({1,0,0},{0,0.2,0})};
    BondedParticleModel model;
    model.normal_stiffness=100.0;
    model.shear_stiffness=100.0;
    model.normal_damping=0.0;
    model.shear_damping=0.0;
    model.tensile_failure_force=100.0;
    model.shear_failure_force=100.0;
    model.damage_onset_ratio=0.8;
    std::vector<ParticleBond> bonds{make_particle_bond(bodies,0,1,model)};
    const auto first=apply_particle_bonds(bodies,bonds,0.1);
    require(first.active==1,"bond shear first update active");
    require_close(bonds[0].tangential_displacement.y,0.02,1e-12,"bond tangential history accumulation");
    require_close(bodies[0].force.y,2.0,1e-12,"bond shear restoring force");
    bodies[0].clear_accumulators();bodies[1].clear_accumulators();
    bodies[1].state.linear_velocity={0,0,0};
    const auto second=apply_particle_bonds(bodies,bonds,0.1);
    require(second.active==1,"bond shear second update active");
    require_close(bodies[0].force.y,2.0,1e-12,"bond shear history must persist without slip");
}

void test_explicit_dem_bond_integration(){
    DemConfig cfg;
    cfg.dt=1.0e-3;
    cfg.gravity={0,0,0};
    ExplicitDemSystem dem(cfg);
    const auto a=dem.add_particle({0,0,0},0.1,1000.0);
    const auto b=dem.add_particle({0.25,0,0},0.1,1000.0);
    BondedParticleModel model;
    model.normal_stiffness=1000.0;
    model.tensile_failure_force=1.0e6;
    model.shear_failure_force=1.0e6;
    dem.add_bond(a,b,model,0.2);
    dem.step();
    require(dem.bonds().size()==1,"DEM bond registration");
    require(dem.last_bond_stats().active==1,"DEM bond update integration");
    require(dem.bodies()[a].state.linear_velocity.x>0.0 && dem.bodies()[b].state.linear_velocity.x<0.0,"DEM bond must pull stretched particles together");
}

void test_slab_exchange_planning(){
    const DemSlabDecomposition decomp{0.0,4.0,4};
    std::vector<DistributedDemParticle> particles(4);
    particles[0].global_id=10;particles[0].state.position={1.05,0,0};particles[0].radius=0.05;particles[0].owner_rank=1;
    particles[1].global_id=11;particles[1].state.position={1.95,0,0};particles[1].radius=0.05;particles[1].owner_rank=1;
    particles[2].global_id=12;particles[2].state.position={2.1,0,0};particles[2].radius=0.05;particles[2].owner_rank=1;
    particles[3].global_id=13;particles[3].state.position={1.5,0,0};particles[3].radius=0.05;particles[3].owner_rank=1;
    const auto plan=plan_dem_slab_exchange(particles,decomp,1,0.05);
    require(plan.migrate_indices_by_rank[2].size()==1 && plan.migrate_indices_by_rank[2][0]==2,"DEM migration destination planning");
    require(plan.ghost_indices_by_rank[0].size()==1 && plan.ghost_indices_by_rank[0][0]==0,"DEM left ghost planning");
    require(plan.ghost_indices_by_rank[2].size()==1 && plan.ghost_indices_by_rank[2][0]==1,"DEM right ghost planning");
    require(plan.ghost_indices_by_rank[1].empty(),"DEM must not ghost to self");
    require(decomp.owner_rank(-1.0)==0 && decomp.owner_rank(4.5)==3,"DEM out-of-domain ownership clamps to edge ranks");
}

void test_distributed_particle_conversion(){
    auto body=unit_body({0.5,0.25,-0.1},{1,2,3});
    body.mass=2.5;body.inertia_diagonal={0.4,0.5,0.6};
    SphereShape sphere{0,0.125};
    const auto p=make_distributed_dem_particle(42,body,sphere,3);
    require(p.global_id==42 && p.owner_rank==3,"distributed DEM particle identity");
    require_close(p.state.linear_velocity.z,3.0,0.0,"distributed DEM particle state copy");
    require_close(p.radius,0.125,0.0,"distributed DEM particle radius copy");
}

}

int main(){
    try{
        test_bond_elastic_damage_and_failure();
        test_bond_shear_history();
        test_explicit_dem_bond_integration();
        test_slab_exchange_planning();
        test_distributed_particle_conversion();
    }catch(const std::exception& e){
        std::cerr<<"v0.16.3 bonded/distributed DEM regression failed: "<<e.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.16.3 bonded/distributed DEM regression passed\n";
    return 0;
}
