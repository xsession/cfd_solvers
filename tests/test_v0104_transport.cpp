#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);} 
std::size_t flat(std::size_t ix,std::size_t iy,std::size_t iz,std::size_t nx,std::size_t ny){return (iz*ny+iy)*nx+ix;} 

double cell_value(std::size_t gx,std::size_t gy,std::size_t gz){return 100.0*static_cast<double>(gx)+10.0*static_cast<double>(gy)+static_cast<double>(gz);} 

cfd::particle::PicParticle3D make_particle(double x,double y,double z){
    cfd::particle::PicParticle3D particle;particle.position_m={x,y,z};particle.velocity_m_per_s={1.0,2.0,3.0};particle.charge_c=1.0;particle.mass_kg=1.0;particle.weight=1.0;return particle;
}

void rank_topology_and_envelopes(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=8U;config.global_ny=4U;config.global_nz=2U;config.domains_x=4U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=8.0;config.length_y_m=4.0;config.length_z_m=2.0;
    const auto domains=make_pic_domain_grid_3d(config);
    const auto topology=make_pic_rank_topology_3d(domains,2U,2U,1U);
    check(topology.rank_count==4U,"rank topology rank count");
    check(rank_for_pic_domain_3d(topology,0U)==0U,"domain zero on rank zero");
    check(rank_for_pic_domain_3d(topology,3U)==1U,"domain x-plus block on rank one");
    check(rank_for_pic_domain_3d(topology,4U)==2U,"domain y-plus block on rank two");

    std::vector<std::vector<PicParticle3D>> particles_by_domain(domains.size());
    particles_by_domain[0].push_back(make_particle(2.1,0.2,0.2)); // rank 0 -> rank 1
    particles_by_domain[1].push_back(make_particle(1.1,0.2,0.2)); // rank 0 -> rank 0, different domain
    particles_by_domain[7].push_back(make_particle(0.1,3.8,0.2)); // rank 3 -> rank 2/periodic position
    const auto migration=pack_particle_migration_messages_3d(particles_by_domain,domains,config);
    const auto envelopes=make_pic_transport_envelopes_3d(migration.messages,{},topology);
    const auto diagnostics=summarize_pic_transport_3d(envelopes);
    check(diagnostics.particle_messages==migration.messages.size(),"particle message diagnostics count");
    check(diagnostics.particle_payload_count==3U,"transport particle payload count");
    check(diagnostics.same_rank_messages>0U&&diagnostics.remote_rank_messages>0U,"transport distinguishes same and remote rank messages");

    InMemoryPicTransport3D transport(topology.rank_count);transport.post(envelopes);
    check(transport.pending_count()==envelopes.size(),"transport pending count after post");
    std::size_t received=0U;
    for(std::size_t rank=0;rank<topology.rank_count;++rank)received+=transport.receive_rank(rank).size();
    check(received==envelopes.size()&&transport.pending_count()==0U,"transport receive drains all inboxes");
}

void transported_exchange_matches_serial(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& domain:domains){
        values[domain.domain_id].resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            values[domain.domain_id][flat(ix,iy,iz,domain.cells_x,domain.cells_y)]=cell_value(domain.first_x+ix,domain.first_y+iy,domain.first_z+iz);
        }
    }
    std::vector<std::vector<PicParticle3D>> particles_by_domain(domains.size());
    particles_by_domain[0].push_back(make_particle(0.2,0.2,0.2));
    particles_by_domain[0].push_back(make_particle(2.2,0.2,0.2));
    particles_by_domain[1].push_back(make_particle(4.2,2.2,1.2));
    particles_by_domain[5].push_back(make_particle(6.1,0.2,0.2));
    const auto migration=pack_particle_migration_messages_3d(particles_by_domain,domains,config);
    const auto guard_messages=pack_scalar_guard_cell_messages_3d(values,domains,config,1U);
    const auto serial_particles=apply_particle_migration_messages_3d(migration.retained_by_domain,migration.messages);
    const auto serial_guards=apply_scalar_guard_cell_messages_3d(values,domains,1U,-7.0,guard_messages);
    const auto topology=make_pic_rank_topology_3d(domains,3U,1U,1U);
    const auto transported=exchange_pic_messages_in_memory_3d(migration.retained_by_domain,migration.messages,values,guard_messages,domains,1U,-7.0,topology);
    check(transported.particles_by_domain.size()==serial_particles.size(),"transport particle bucket size");
    std::size_t transported_count=0U,serial_count=0U;
    for(std::size_t domain=0;domain<serial_particles.size();++domain){
        transported_count+=transported.particles_by_domain[domain].size();serial_count+=serial_particles[domain].size();
    }
    check(transported_count==serial_count&&transported_count==4U,"transported migration preserves particle count");
    check(transported.guarded_blocks.size()==serial_guards.size(),"transport guarded block count");
    for(std::size_t domain=0;domain<serial_guards.size();++domain){
        check(transported.guarded_blocks[domain].values.size()==serial_guards[domain].values.size(),"transport guard value size");
        for(std::size_t i=0;i<serial_guards[domain].values.size();++i){
            check(std::abs(transported.guarded_blocks[domain].values[i]-serial_guards[domain].values[i])<1.0e-12,"transport guard parity");
        }
    }
    check(transported.diagnostics.scalar_guard_messages==guard_messages.size(),"transport scalar-message count");
    check(transported.diagnostics.scalar_guard_value_count>0U,"transport scalar payload count");
    check(transported.diagnostics.remote_rank_messages>0U,"rank transport has remote messages");
}
}

int main(){
    try{
        rank_topology_and_envelopes();
        transported_exchange_matches_serial();
    }catch(const std::exception& error){
        std::cerr<<"v0.10.4 transport test failed: "<<error.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.10.4 transport tests passed\n";
    return 0;
}
