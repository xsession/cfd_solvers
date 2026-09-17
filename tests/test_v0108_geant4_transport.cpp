#include "cfd/particle/transport.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);} 

cfd::particle::TransportWorld make_world(){
    using namespace cfd::particle;
    TransportWorld world;
    world.particles.push_back({"electron",-1.602176634e-19,9.1093837e-31});
    world.materials.push_back({"gas",1.2,5.0});
    world.materials.push_back({"absorber",2000.0,20.0});
    world.regions.push_back({"gas-slab",{{0.0,0.0,0.0},{1.0,1.0,1.0}},0U});
    world.regions.push_back({"absorber-slab",{{1.0,0.0,0.0},{2.0,1.0,1.0}},1U});
    world.processes.push_back({"ionization-loss",TransportProcessKind::continuous_energy_loss,0U,0U,100.0,0.0,0.0,0.0,0.0});
    world.processes.push_back({"hard-scatter",TransportProcessKind::discrete_interaction,0U,0U,0.0,0.25,0.2,0.25,0.0});
    world.processes.push_back({"low-energy-absorb",TransportProcessKind::absorber,0U,1U,0.0,0.0,0.0,0.0,25.0});
    return world;
}

void continuous_discrete_and_secondary(){
    using namespace cfd::particle;
    auto world=make_world();TransportConfig config;config.max_step_m=0.8;config.energy_cut_ev=1.0;
    TransportTrack track;track.position_m={0.1,0.5,0.5};track.direction={1.0,0.0,0.0};track.kinetic_energy_ev=1000.0;
    const auto step=transport_one_step(world,track,config);
    check(step.limiter==StepLimiter::discrete_process,"discrete process must define the first step");
    check(std::abs(step.step_length_m-0.25)<1.0e-14,"physical interaction length selected");
    check(std::abs(step.energy_deposit_ev-(25.0+0.2*975.0))<1.0e-12,"continuous and discrete energy deposits");
    check(step.secondaries.size()==1U,"discrete interaction creates secondary above production cut");
    check(std::abs(step.secondaries[0].track.kinetic_energy_ev-0.25*(975.0-0.2*975.0))<1.0e-12,"secondary energy fraction");
    check(track.status==TrackStatus::alive,"primary remains alive after hard scatter");
}

void geometry_boundary_and_absorber(){
    using namespace cfd::particle;
    auto world=make_world();TransportConfig config;config.max_step_m=2.0;config.energy_cut_ev=1.0;
    // Disable the discrete process so the boundary is the limiting GPIL.
    world.processes[1].physical_interaction_length_m=0.0;
    TransportTrack track;track.position_m={0.9,0.5,0.5};track.direction={1.0,0.0,0.0};track.kinetic_energy_ev=40.0;
    auto step=transport_one_step(world,track,config);
    check(step.limiter==StepLimiter::geometry_boundary,"geometry boundary must limit transport step");
    check(step.material_index==0U,"first step is in gas material");
    check(track.status==TrackStatus::alive,"track should enter adjacent slab after crossing internal boundary");
    step=transport_one_step(world,track,config);
    check(step.material_index==1U,"second step is in absorber material");
    check(track.status==TrackStatus::alive||track.status==TrackStatus::stopped||track.status==TrackStatus::escaped,"valid post-absorber status");
}

void run_scoring(){
    using namespace cfd::particle;
    auto world=make_world();TransportConfig config;config.max_step_m=0.2;config.energy_cut_ev=10.0;
    world.processes[1].physical_interaction_length_m=0.0;
    TransportTrack track;track.position_m={0.1,0.5,0.5};track.direction={1.0,0.0,0.0};track.kinetic_energy_ev=80.0;
    const auto result=transport_track(world,track,config,16U);
    check(!result.steps.empty(),"transport run records steps");
    check(result.scoring.steps==result.steps.size(),"step count scoring");
    check(result.scoring.total_track_length_m>0.0,"track length accumulates");
    check(result.scoring.total_energy_deposit_ev>0.0,"energy deposition accumulates");
    check(result.primary.status!=TrackStatus::alive||result.steps.size()==16U,"run terminates or exhausts budget");
}

void invalid_inputs(){
    using namespace cfd::particle;
    bool caught=false;try{(void)normalize_direction({0.0,0.0,0.0});}catch(const std::exception&){caught=true;}check(caught,"zero direction rejected");
    auto world=make_world();TransportTrack bad;bad.particle_index=4U;bad.direction={1.0,0.0,0.0};bad.kinetic_energy_ev=1.0;TransportConfig config;
    caught=false;try{(void)transport_track(world,bad,config,1U);}catch(const std::exception&){caught=true;}check(caught,"invalid particle index rejected");
}
}

int main(){
    continuous_discrete_and_secondary();
    geometry_boundary_and_absorber();
    run_scoring();
    invalid_inputs();
    std::cout<<"v0.10.8 Geant4-inspired transport tests passed\n";
    return 0;
}
