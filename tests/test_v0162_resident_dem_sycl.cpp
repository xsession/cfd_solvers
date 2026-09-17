#include "cfd/multibody/contact.hpp"
#include "cfd/multibody/resident_dem_sycl.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfd::multibody;

namespace {
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}

void test_cell_linked_cpu_reference(){
    std::vector<RigidBody> bodies(7);
    const std::vector<Vec3> positions{{0,0,0},{0.09,0,0},{0.31,0.01,0},{0.39,0.02,0},{0.0,0.25,0},{0.08,0.27,0},{0.65,0.65,0.65}};
    std::vector<SphereShape> spheres;
    for(std::size_t i=0;i<bodies.size();++i){bodies[i].state.position=positions[i];spheres.push_back({i,0.05});}
    auto sweep=broad_phase_sweep_and_prune(bodies,spheres);
    auto grid=broad_phase_cell_linked(bodies,spheres,0.12);
    std::sort(sweep.begin(),sweep.end());std::sort(grid.begin(),grid.end());
    require(sweep==grid,"cell-linked broad phase must match sweep-and-prune AABB candidates");
    require(grid.size()==3U,"expected three close particle pairs");
}

#if defined(CFD_HAS_SYCL)
void test_resident_dem_hot_loop(){
    sycl::queue queue{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};
    ResidentDemSyclConfig cfg;
    cfg.capacity=8U;cfg.domain_min={-1,-1,-1};cfg.domain_max={1,1,1};cfg.cell_size=0.2;cfg.dt=1.0e-4;cfg.gravity={0,0,0};
    cfg.normal_stiffness=2.0e5;cfg.normal_damping=5.0;cfg.tangential_stiffness=1.0e4;cfg.tangential_damping=2.0;cfg.friction=0.5;cfg.history_slots_per_particle=8U;
    ResidentDemSycl dem(queue,cfg);
    std::vector<ResidentDemParticle> particles(3);
    particles[0].position={-0.045,0,0};particles[0].radius=0.05;particles[0].mass=1.0;particles[0].linear_velocity={0,0.1,0};
    particles[1].position={ 0.045,0,0};particles[1].radius=0.05;particles[1].mass=1.0;particles[1].linear_velocity={0,-0.1,0};
    particles[2].position={ 0.5,0.5,0.5};particles[2].radius=0.04;particles[2].mass=0.5;
    dem.upload(particles);
    dem.reset_transfer_stats();
    dem.step(4U);
    require(dem.transfer_stats().host_transfer_bytes()==0U,"resident DEM hot loop must not transfer particle fields through host");
    require(dem.total_kinetic_energy()>0.0,"resident DEM kinetic-energy reduction");
    require(dem.transfer_stats().host_transfer_bytes()==0U,"scalar reduction must not count bulk field host transfers");
    const auto out=dem.download();
    require(out.size()==3U,"resident DEM download size");
    require(out[0].position.x<particles[0].position.x && out[1].position.x>particles[1].position.x,"Hertz contact should separate overlapping particles");
}
#endif
}

int main(){
    try{
        test_cell_linked_cpu_reference();
#if defined(CFD_HAS_SYCL)
        test_resident_dem_hot_loop();
#endif
    }catch(const std::exception& e){
        std::cerr<<"v0.16.2 resident DEM regression failed: "<<e.what()<<'\n';return 1;
    }
    std::cout<<"v0.16.2 resident DEM regression passed\n";return 0;
}
