#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);} 
bool near(double a,double b,double tol=1.0e-12){return std::abs(a-b)<=tol*std::max({1.0,std::abs(a),std::abs(b)});} 
std::size_t idx3(std::size_t x,std::size_t y,std::size_t z,std::size_t nx,std::size_t ny){return (z*ny+y)*nx+x;}

std::vector<std::vector<double>> make_component(const std::vector<cfd::particle::PicDomain3D>& domains,double value){
    std::vector<std::vector<double>> result(domains.size());
    for(const auto& domain:domains)result[domain.domain_id].assign(domain.cells_x*domain.cells_y*domain.cells_z,value);
    return result;
}

void distributed_step_pushes_and_migrates(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=4U;config.global_ny=2U;config.global_nz=2U;config.domains_x=2U;config.domains_y=1U;config.domains_z=1U;config.length_x_m=4.0;config.length_y_m=2.0;config.length_z_m=2.0;config.periodic=true;
    const auto domains=make_pic_domain_grid_3d(config);
    ElectromagneticFieldBlocks3D fields;
    fields.electric_x_by_domain=make_component(domains,0.0);fields.electric_y_by_domain=make_component(domains,0.0);fields.electric_z_by_domain=make_component(domains,0.0);
    fields.magnetic_x_by_domain=make_component(domains,0.0);fields.magnetic_y_by_domain=make_component(domains,0.0);fields.magnetic_z_by_domain=make_component(domains,0.0);
    // Put a marker in Ex so the guard-field exchange can be checked independently
    // from the zero-field ballistic push.
    for(const auto& domain:domains){
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            fields.electric_x_by_domain[domain.domain_id][idx3(ix,iy,iz,domain.cells_x,domain.cells_y)]=100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz);
        }
    }
    std::vector<std::vector<PicParticle3D>> particles(domains.size());
    PicParticle3D left;left.position_m={1.75,0.6,0.6};left.velocity_m_per_s={0.50,0.0,0.0};left.mass_kg=1.0;left.charge_c=0.0;left.weight=1.0;particles[0].push_back(left);
    PicParticle3D right;right.position_m={3.85,0.5,0.5};right.velocity_m_per_s={0.40,0.0,0.0};right.mass_kg=1.0;right.charge_c=0.0;right.weight=1.0;particles[1].push_back(right);
    DistributedStaggeredPicStep3DConfig step;step.dt_s=1.0;step.guard_cells=1U;step.exterior_value=-7.0;
    const auto topology=make_pic_rank_topology_3d(domains,2U,1U,1U);
    const auto result=run_serialized_distributed_staggered_pic_step_3d(particles,fields,domains,config,step,topology);
    require(result.particle_count_before==2U&&result.particle_count_after_push==2U&&result.particle_count_after_migration==2U,"distributed step particle counts");
    require(result.max_particle_displacement_m>0.39&&result.max_particle_displacement_m<0.51,"distributed step records ballistic displacement");
    require(result.exchange_round.transport_diagnostics.particle_payload_count==2U,"distributed step moved both particles through migration messages");
    require(result.exchange_round.transport_diagnostics.remote_rank_messages>=1U,"distributed step remote rank message accounting");
    require(result.particles_by_domain[0].size()==1U&&result.particles_by_domain[1].size()==1U,"distributed step ownership after migration");
    const double wrapped_x=result.particles_by_domain[0][0].position_m.x;
    require(near(wrapped_x,0.25,1.0e-12),"periodic migrating particle wrapped into domain zero");
    require(result.guarded_fields.electric_x.size()==domains.size(),"distributed step returns EM guard fields");
    const auto& guarded0=result.guarded_fields.electric_x[0];
    require(guarded0.values.size()==guarded0.padded_nx()*guarded0.padded_ny()*guarded0.padded_nz(),"guarded field padded extents");
    require(result.exchange_round.serialized_message_count==result.exchange_round.transport_diagnostics.envelopes,"distributed step serialized diagnostic parity");
}

void distributed_step_rejects_bad_inputs(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=4U;config.global_ny=2U;config.global_nz=2U;config.domains_x=2U;config.length_x_m=4.0;config.length_y_m=2.0;config.length_z_m=2.0;
    const auto domains=make_pic_domain_grid_3d(config);
    ElectromagneticFieldBlocks3D fields;fields.electric_x_by_domain=make_component(domains,0.0);
    fields.electric_y_by_domain=make_component(domains,0.0);fields.electric_z_by_domain=make_component(domains,0.0);fields.magnetic_x_by_domain=make_component(domains,0.0);fields.magnetic_y_by_domain=make_component(domains,0.0);fields.magnetic_z_by_domain=make_component(domains,0.0);
    fields.electric_x_by_domain[0].pop_back();
    bool caught=false;
    try{
        DistributedStaggeredPicStep3DConfig step;step.dt_s=1.0;
        (void)run_serialized_distributed_staggered_pic_step_3d(std::vector<std::vector<PicParticle3D>>(domains.size()),fields,domains,config,step,make_pic_rank_topology_3d(domains,2U,1U,1U));
    }catch(const std::exception&){caught=true;}
    require(caught,"distributed step rejects malformed field block");
}
}

int main(){
    try{
        distributed_step_pushes_and_migrates();
        distributed_step_rejects_bad_inputs();
        std::cout<<"v0.10.7 distributed PIC step tests passed\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"v0.10.7 distributed PIC step test failure: "<<e.what()<<"\n";
        return 1;
    }
}
