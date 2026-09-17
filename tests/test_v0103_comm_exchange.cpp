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

void communicator_particle_packets(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=2U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=2.0;
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<PicParticle3D>> source(domains.size());
    auto make_particle=[](double x,double y,double z){PicParticle3D p;p.position_m={x,y,z};p.velocity_m_per_s={1.0,2.0,3.0};p.charge_c=1.0;p.mass_kg=1.0;p.weight=1.0;return p;};
    source[0].push_back(make_particle(0.25,0.25,0.25));     // retained in domain 0
    source[0].push_back(make_particle(2.25,0.25,0.25));     // migrates to +x neighbour
    source[0].push_back(make_particle(5.95,3.75,0.25));     // migrates to far wrapped domain by position
    source[5].push_back(make_particle(6.10,0.25,0.25));     // wraps to domain 0

    const auto packed=pack_particle_migration_messages_3d(source,domains,config);
    check(packed.outside_particles.empty(),"periodic communicator pack keeps all particles");
    check(packed.retained_by_domain[0].size()==1U,"one particle retained in source domain");
    std::size_t message_particles=0U;for(const auto& message:packed.messages){message_particles+=message.particles.size();check(message.particles.size()==message.source_indices.size(),"message source-index parity");check(message.source_domain_id!=message.destination_domain_id,"no self messages emitted");}
    check(message_particles==3U,"three particles emitted as communicator messages");
    const auto after=apply_particle_migration_messages_3d(packed.retained_by_domain,packed.messages);
    std::size_t total=0U;for(const auto& bucket:after)total+=bucket.size();check(total==4U,"apply messages preserves particle count");
    check(!after[0].empty(),"wrapped incoming particle reaches domain zero");
    bool found_wrapped=false;for(const auto& particle:after[0])if(particle.position_m.x<0.2)found_wrapped=true;check(found_wrapped,"periodic source particle is wrapped before send");

    config.periodic=false;
    const auto open=pack_particle_migration_messages_3d(source,domains,config);
    check(open.outside_particles.size()==1U,"nonperiodic communicator pack reports outside particle");
    check(open.outside_particles[0].source_domain_id==5U&&open.outside_particles[0].source_index==0U,"outside record identifies source domain and local index");
}

void communicator_guard_packets(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=5U;config.global_ny=4U;config.global_nz=3U;config.domains_x=2U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=5.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& d:domains){
        values[d.domain_id].resize(d.cells_x*d.cells_y*d.cells_z);
        for(std::size_t iz=0;iz<d.cells_z;++iz)for(std::size_t iy=0;iy<d.cells_y;++iy)for(std::size_t ix=0;ix<d.cells_x;++ix){
            values[d.domain_id][flat(ix,iy,iz,d.cells_x,d.cells_y)]=cell_value(d.first_x+ix,d.first_y+iy,d.first_z+iz);
        }
    }
    const auto reference=exchange_scalar_guard_cells_3d(values,domains,config,1U,-9.0);
    const auto messages=pack_scalar_guard_cell_messages_3d(values,domains,config,1U);
    check(!messages.empty(),"guard exchange produces communicator messages");
    bool has_corner=false,has_face=false;std::size_t entries=0U;
    for(const auto& message:messages){
        entries+=message.values.size();
        const int manhattan=std::abs(message.offset_x)+std::abs(message.offset_y)+std::abs(message.offset_z);
        if(manhattan==1)has_face=true; if(manhattan>=2)has_corner=true;
    }
    check(has_face&&has_corner&&entries>0U,"face and edge/corner guard messages are represented");
    const auto applied=apply_scalar_guard_cell_messages_3d(values,domains,1U,-9.0,messages);
    check(applied.size()==reference.size(),"applied guard block count");
    for(std::size_t d=0;d<reference.size();++d){
        check(applied[d].values.size()==reference[d].values.size(),"guard block size parity");
        for(std::size_t i=0;i<reference[d].values.size();++i)check(std::abs(applied[d].values[i]-reference[d].values[i])<1.0e-12,"communicator guard messages reproduce serial exchange");
    }

    config.periodic=false;
    const auto open_messages=pack_scalar_guard_cell_messages_3d(values,domains,config,1U);
    const auto open_applied=apply_scalar_guard_cell_messages_3d(values,domains,1U,-123.0,open_messages);
    const auto& first=domains.front();const auto& block=open_applied[first.domain_id];
    check(block.values[flat(0U,1U,1U,block.padded_nx(),block.padded_ny())]==-123.0,"exterior cells remain fill value when no message arrives");
}
}

int main(){
    try{
        communicator_particle_packets();
        communicator_guard_packets();
    }catch(const std::exception& error){
        std::cerr<<"v0.10.3 communicator exchange test failed: "<<error.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.10.3 communicator exchange tests passed\n";
    return 0;
}
