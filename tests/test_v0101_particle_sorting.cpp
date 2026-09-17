#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

std::vector<cfd::particle::PicParticle3D> make_particles(){
    using cfd::particle::PicParticle3D;
    std::vector<PicParticle3D> particles;
    const double charges[] = {1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0};
    for(std::size_t i=0;i<9U;++i){
        PicParticle3D p;
        p.mass_kg=1.0e-6;
        p.charge_c=charges[i]*1.0e-15;
        p.weight=1.0+0.1*static_cast<double>(i);
        p.position_m={0.91-0.083*static_cast<double>(i),0.07+0.121*static_cast<double>((i*3U)%7U),0.03+0.097*static_cast<double>((i*5U)%8U)};
        p.velocity_m_per_s={2.0e3,-1.0e3,1.5e3};
        particles.push_back(p);
    }
    particles.front().position_m.x=-0.02;       // periodic wrapping path
    particles.back().position_m.z=1.04;        // periodic wrapping path
    return particles;
}

void sorted_ranges_are_stable_and_recoverable(){
    using namespace cfd::particle;
    constexpr std::size_t nx=5U,ny=4U,nz=3U;
    auto particles=make_particles();
    const auto sorted=sort_particles_by_cell_3d(particles,nx,ny,nz,1.0,0.8,0.6);
    require(sorted.sorted_particles.size()==particles.size(),"sorted particle count is preserved");
    require(sorted.original_indices.size()==particles.size(),"original index count is preserved");
    require(sorted.cell_offsets.size()==nx*ny*nz+1U,"cell offset table has sentinel");
    require(sorted.cell_offsets.front()==0U&&sorted.cell_offsets.back()==particles.size(),"offset table spans all particles");
    std::size_t counted=0U,occupied=0U;
    for(std::size_t cell=0;cell<nx*ny*nz;++cell){
        const auto first=sorted.cell_offsets[cell];
        const auto last=sorted.cell_offsets[cell+1U];
        require(last>=first,"cell offsets are monotone");
        require(last-first==sorted.cell_counts[cell],"cell count matches offsets");
        counted+=sorted.cell_counts[cell];
        if(sorted.cell_counts[cell]>0U)++occupied;
        for(std::size_t pos=first;pos<last;++pos){
            require(sorted.original_indices[pos]<particles.size(),"original index is valid");
            const auto recovered=particles[sorted.original_indices[pos]];
            require(std::abs(recovered.charge_c-sorted.sorted_particles[pos].charge_c)<1.0e-30,"stable sort keeps particle payload");
            const auto sorted_cell=particle_cell_index_3d(sorted.sorted_particles[pos].position_m,nx,ny,nz,1.0,0.8,0.6);
            require(sorted_cell==cell,"sorted particle lies in its cell range");
        }
    }
    require(counted==particles.size(),"all particles counted by cell ranges");
    require(occupied==sorted.occupied_cells,"occupied cell diagnostic matches counts");
}

void guard_halo_classification_handles_faces_edges_and_outside(){
    using namespace cfd::particle;
    std::vector<PicParticle3D> particles(5U);
    for(auto& p:particles){p.mass_kg=1.0;p.charge_c=1.0;p.weight=1.0;}
    particles[0].position_m={0.02,0.50,0.50}; // x-
    particles[1].position_m={0.98,0.98,0.50}; // x+ and y+
    particles[2].position_m={0.50,0.50,0.50}; // interior
    particles[3].position_m={-0.01,0.20,0.20}; // outside for nonperiodic
    particles[4].position_m={0.40,0.02,0.98}; // y- and z+
    ParticleGuardHalo3DConfig config;config.nx=10U;config.ny=10U;config.nz=10U;config.guard_cells=1U;config.periodic=false;
    const auto halo=classify_particle_guard_halos_3d(particles,config);
    require(std::find(halo.x_minus.begin(),halo.x_minus.end(),0U)!=halo.x_minus.end(),"x- halo classified");
    require(std::find(halo.x_plus.begin(),halo.x_plus.end(),1U)!=halo.x_plus.end(),"x+ halo classified");
    require(std::find(halo.y_plus.begin(),halo.y_plus.end(),1U)!=halo.y_plus.end(),"edge particle appears in y+ halo");
    require(std::find(halo.interior.begin(),halo.interior.end(),2U)!=halo.interior.end(),"interior particle classified");
    require(std::find(halo.outside_domain.begin(),halo.outside_domain.end(),3U)!=halo.outside_domain.end(),"nonperiodic outside particle classified");
    require(std::find(halo.y_minus.begin(),halo.y_minus.end(),4U)!=halo.y_minus.end()&&std::find(halo.z_plus.begin(),halo.z_plus.end(),4U)!=halo.z_plus.end(),"corner halo classified on both faces");
}

void staggered_pic_sorts_particles_after_step(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic3DConfig config;
    config.nx=6U;config.ny=5U;config.nz=4U;
    config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;
    config.current_deposition=CurrentDeposition3DMode::local_finite_volume;
    config.sort_particles_by_cell=true;
    config.particle_sort_interval=1U;
    const double dx=config.length_x_m/static_cast<double>(config.nx);
    const double dy=config.length_y_m/static_cast<double>(config.ny);
    const double dz=config.length_z_m/static_cast<double>(config.nz);
    config.dt_s=0.02/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
    StaggeredElectromagneticPic3D solver(config);
    solver.set_particles(make_particles());
    solver.step(1U);
    const auto diag=solver.diagnostics();
    require(diag.particle_count==make_particles().size(),"sorting preserves particle count");
    require(diag.particle_sort_passes==1U,"one sort pass recorded");
    require(diag.occupied_particle_cells>0U,"occupied cell diagnostic set");
    std::size_t previous=0U;
    bool first=true;
    for(const auto& particle:solver.particles()){
        const auto cell=particle_cell_index_3d(particle.position_m,config.nx,config.ny,config.nz,config.length_x_m,config.length_y_m,config.length_z_m);
        if(!first)require(cell>=previous,"particles are sorted by cell after the step");
        first=false;previous=cell;
    }
    require(diag.charge_continuity_linf_residual<1.0e-6,"sorting mode preserves local current continuity");
}
}

int main(){
    try{
        sorted_ranges_are_stable_and_recoverable();
        guard_halo_classification_handles_faces_edges_and_outside();
        staggered_pic_sorts_particles_after_step();
        std::cout<<"v0.10.1 particle sorting and guard-halo tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.10.1 particle sorting test failure: "<<error.what()<<'\n';
        return 1;
    }
}
