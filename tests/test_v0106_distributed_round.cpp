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

double global_value(std::size_t x,std::size_t y,std::size_t z,double offset=0.0){return offset+100.0*static_cast<double>(x)+10.0*static_cast<double>(y)+static_cast<double>(z);} 

std::vector<std::vector<double>> make_values(const std::vector<cfd::particle::PicDomain3D>& domains,double offset=0.0){
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& domain:domains){
        auto& block=values[domain.domain_id];
        block.resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz){
            for(std::size_t iy=0;iy<domain.cells_y;++iy){
                for(std::size_t ix=0;ix<domain.cells_x;++ix){
                    block[idx3(ix,iy,iz,domain.cells_x,domain.cells_y)]=global_value(domain.first_x+ix,domain.first_y+iy,domain.first_z+iz,offset);
                }
            }
        }
    }
    return values;
}

void field_vector_guard_exchange(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto domains=make_pic_domain_grid_3d(config);
    ElectromagneticFieldBlocks3D fields;
    fields.electric_x_by_domain=make_values(domains,1000.0);
    fields.electric_y_by_domain=make_values(domains,2000.0);
    fields.electric_z_by_domain=make_values(domains,3000.0);
    fields.magnetic_x_by_domain=make_values(domains,4000.0);
    fields.magnetic_y_by_domain=make_values(domains,5000.0);
    fields.magnetic_z_by_domain=make_values(domains,6000.0);
    const auto guarded=exchange_electromagnetic_field_guard_cells_3d(fields,domains,config,1U,-99.0);
    const auto reference_ex=exchange_scalar_guard_cells_3d(fields.electric_x_by_domain,domains,config,1U,-99.0);
    require(guarded.electric_x.size()==domains.size()&&guarded.magnetic_z.size()==domains.size(),"six-component EM guard exchange block count");
    for(std::size_t domain=0;domain<domains.size();++domain){
        require(guarded.electric_x[domain].values==reference_ex[domain].values,"EM Ex guard equals scalar reference");
        const auto& d=domains[domain];
        const auto px=d.cells_x+2U,py=d.cells_y+2U;
        const double interior_ez=guarded.electric_z[domain].values[idx3(1U,1U,1U,px,py)];
        const double interior_bz=guarded.magnetic_z[domain].values[idx3(1U,1U,1U,px,py)];
        require(near(interior_ez,global_value(d.first_x,d.first_y,d.first_z,3000.0)),"EM Ez interior copy");
        require(near(interior_bz,global_value(d.first_x,d.first_y,d.first_z,6000.0)),"EM Bz interior copy");
    }
}

void serialized_distributed_round(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto domains=make_pic_domain_grid_3d(config);
    auto values=make_values(domains,0.0);
    std::vector<std::vector<PicParticle3D>> particles_by_domain(domains.size());
    auto add=[&](std::size_t source,double x,double y,double z){PicParticle3D p;p.position_m={x,y,z};p.velocity_m_per_s={1.0e3,2.0e3,3.0e3};p.mass_kg=1.0e-6;p.charge_c=1.0e-15;p.weight=1.0;particles_by_domain[source].push_back(p);};
    add(0U,0.25,0.25,0.25);       // retained
    add(0U,2.25,0.25,0.25);       // migrates from domain 0 to 1
    add(1U,4.25,2.25,1.25);       // migrates to upper-y/right-x domain
    add(domains.back().domain_id,6.05,0.25,0.25); // periodic wrap to domain 0
    const auto topology=make_pic_rank_topology_3d(domains,3U,1U,1U);
    const auto round=run_serialized_distributed_pic_exchange_round_3d(particles_by_domain,values,domains,config,1U,-1.0,topology);
    std::size_t particles=0U;for(const auto& bucket:round.particles_by_domain)particles+=bucket.size();
    require(particles==4U,"distributed round preserves all periodic particles");
    require(round.transport_diagnostics.envelopes>0U&&round.transport_diagnostics.remote_rank_messages>0U,"distributed round records remote transport");
    require(round.serialized_message_count==round.transport_diagnostics.envelopes,"serialized round message count parity");
    require(round.serialized_byte_count>0U,"serialized round byte accounting");
    const auto direct_guard=exchange_scalar_guard_cells_3d(values,domains,config,1U,-1.0);
    require(round.scalar_guarded_blocks.size()==direct_guard.size(),"serialized round guard block count");
    for(std::size_t domain=0;domain<direct_guard.size();++domain){
        require(round.scalar_guarded_blocks[domain].values==direct_guard[domain].values,"serialized round guard values equal direct exchange");
    }
}
}

int main(){
    try{
        field_vector_guard_exchange();
        serialized_distributed_round();
        std::cout<<"v0.10.6 distributed round tests passed\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"v0.10.6 distributed round test failure: "<<e.what()<<"\n";
        return 1;
    }
}
