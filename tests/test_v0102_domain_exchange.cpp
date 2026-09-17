#include "cfd/particle/electromagnetic.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
std::size_t flat(std::size_t ix,std::size_t iy,std::size_t iz,std::size_t nx,std::size_t ny){return (iz*ny+iy)*nx+ix;}

double cell_value(std::size_t gx,std::size_t gy,std::size_t gz){return 100.0*static_cast<double>(gx)+10.0*static_cast<double>(gy)+static_cast<double>(gz);}

void domain_layout_and_migration(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=7U;config.global_ny=5U;config.global_nz=4U;config.domains_x=3U;config.domains_y=2U;config.domains_z=2U;config.length_x_m=7.0;config.length_y_m=5.0;config.length_z_m=4.0;
    const auto domains=make_pic_domain_grid_3d(config);
    check(domains.size()==12U,"domain count");
    std::size_t cell_sum=0U;for(const auto& d:domains){check(d.cells_x>0U&&d.cells_y>0U&&d.cells_z>0U,"non-empty domain");cell_sum+=d.cells_x*d.cells_y*d.cells_z;}
    check(cell_sum==config.global_nx*config.global_ny*config.global_nz,"domain cells cover global grid exactly");
    check(locate_pic_domain_3d({0.25,0.25,0.25},config)==0U,"first cell maps to first domain");
    const auto far_domain=locate_pic_domain_3d({6.75,4.75,3.75},config);
    check(far_domain==domains.back().domain_id,"last cell maps to last domain");
    check(locate_pic_domain_3d({7.1,0.25,0.25},config)==0U,"periodic position wraps into first x domain");

    std::vector<PicParticle3D> particles(4U);
    particles[0].position_m={0.2,0.2,0.2};
    particles[1].position_m={3.4,3.2,1.2};
    particles[2].position_m={6.9,4.7,3.8};
    particles[3].position_m={7.05,0.1,0.1};
    for(auto& p:particles){p.mass_kg=1.0;p.charge_c=1.0;p.weight=1.0;}
    const auto migration=plan_particle_domain_migration_3d(particles,config);
    check(migration.outside_indices.empty(),"periodic migration keeps wrapped particle");
    check(migration.destination_domain_by_particle.size()==particles.size(),"destination map size");
    check(migration.destination_domain_by_particle[0]==0U,"particle 0 destination");
    check(migration.destination_domain_by_particle[3]==0U,"wrapped particle destination");
    std::size_t migrated=0U;for(const auto& bucket:migration.particles_by_domain)migrated+=bucket.size();
    check(migrated==particles.size(),"all periodic particles bucketed");

    config.periodic=false;
    const auto nonperiodic=plan_particle_domain_migration_3d(particles,config);
    check(nonperiodic.outside_indices.size()==1U&&nonperiodic.outside_indices[0]==3U,"nonperiodic outside particle reported");
}

void scalar_guard_exchange(){
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
    const auto guarded=exchange_scalar_guard_cells_3d(values,domains,config,1U,-1.0);
    check(guarded.size()==domains.size(),"guarded block count");
    for(const auto& d:domains){
        const auto& block=guarded[d.domain_id];
        const std::size_t px=block.padded_nx(),py=block.padded_ny(),pz=block.padded_nz();
        for(std::size_t iz=0;iz<pz;++iz)for(std::size_t iy=0;iy<py;++iy)for(std::size_t ix=0;ix<px;++ix){
            long gx=static_cast<long>(d.first_x)+static_cast<long>(ix)-1L;
            long gy=static_cast<long>(d.first_y)+static_cast<long>(iy)-1L;
            long gz=static_cast<long>(d.first_z)+static_cast<long>(iz)-1L;
            const auto wrap=[](long v,long n){long r=v%n;if(r<0)r+=n;return r;};
            gx=wrap(gx,static_cast<long>(config.global_nx));gy=wrap(gy,static_cast<long>(config.global_ny));gz=wrap(gz,static_cast<long>(config.global_nz));
            const double expected=cell_value(static_cast<std::size_t>(gx),static_cast<std::size_t>(gy),static_cast<std::size_t>(gz));
            const double actual=block.values[flat(ix,iy,iz,px,py)];
            check(std::abs(actual-expected)<1.0e-12,"periodic guard value copied from correct neighbour");
        }
    }

    config.periodic=false;
    const auto open=exchange_scalar_guard_cells_3d(values,domains,config,1U,-123.0);
    const auto& first=domains.front();
    const auto& block=open[first.domain_id];
    check(block.values[flat(0U,1U,1U,block.padded_nx(),block.padded_ny())]==-123.0,"nonperiodic x-minus exterior fill");
    check(block.values[flat(1U,1U,1U,block.padded_nx(),block.padded_ny())]==cell_value(0U,0U,0U),"interior preserved in guarded block");
}
}

int main(){
    try{
        domain_layout_and_migration();
        scalar_guard_exchange();
    }catch(const std::exception& error){
        std::cerr<<"v0.10.2 domain exchange test failed: "<<error.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.10.2 domain exchange tests passed\n";
    return 0;
}
