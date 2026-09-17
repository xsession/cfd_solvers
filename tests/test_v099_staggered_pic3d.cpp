#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

cfd::particle::StaggeredElectromagneticPic3DConfig compact_config(){
    cfd::particle::StaggeredElectromagneticPic3DConfig config;
    config.nx=8U;config.ny=7U;config.nz=6U;config.length_x_m=0.8;config.length_y_m=0.7;config.length_z_m=0.6;
    const double dx=config.length_x_m/static_cast<double>(config.nx);
    const double dy=config.length_y_m/static_cast<double>(config.ny);
    const double dz=config.length_z_m/static_cast<double>(config.nz);
    config.dt_s=0.05/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
    return config;
}

void staggered_vacuum_mode_is_bounded(){
    using namespace cfd::particle;
    auto config=compact_config();
    StaggeredElectromagneticPic3D solver(config);
    solver.initialize_z_polarized_mode(0.4,1U,0U,0U);
    const auto before=solver.diagnostics();
    require(before.field_energy_j>0.0,"staggered 3-D Yee EM-PIC initializes field energy");
    const auto sample=solver.gather_field({0.21,0.19,0.17});
    require(std::isfinite(sample.electric_v_per_m.z)&&std::isfinite(sample.magnetic_t.y),"staggered 3-D Yee EM-PIC gathers finite staggered fields");
    solver.step(5U);
    const auto after=solver.diagnostics();
    require(after.field_energy_j>0.05*before.field_energy_j&&after.field_energy_j<6.0*before.field_energy_j,"staggered 3-D Yee EM-PIC vacuum mode remains bounded");
    require(after.particle_count==0U,"staggered vacuum mode has no particles");
}

void moving_beam_reconstructs_current(){
    using namespace cfd::particle;
    auto config=compact_config();config.dt_s=5.0e-13;
    StaggeredElectromagneticPic3D solver(config);
    std::vector<PicParticle3D> particles;
    for(std::size_t i=0;i<5U;++i){
        PicParticle3D particle;particle.mass_kg=1.0e-6;particle.charge_c=(i%2U? -1.0:1.0)*1.0e-15;particle.position_m={0.09+0.11*static_cast<double>(i),0.08+0.06*static_cast<double>(i),0.05+0.04*static_cast<double>(i)};particle.velocity_m_per_s={2.0e4,-1.5e4,1.0e4};particles.push_back(particle);
    }
    solver.set_particles(particles);solver.step(1U);
    double max_current=0.0;for(double value:solver.current_x())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_y())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_z())max_current=std::max(max_current,std::abs(value));
    const auto diag=solver.diagnostics();
    require(diag.particle_count==particles.size(),"staggered 3-D Yee EM-PIC keeps periodic beam particles");
    require(max_current>0.0,"staggered 3-D Yee EM-PIC couples reconstructed current");
    require(diag.charge_continuity_linf_residual<1.0e-6,"staggered 3-D Yee EM-PIC spectral current satisfies continuity");
}

void wall_and_sponge_primitives_work(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic3DConfig config;config.nx=5U;config.ny=5U;config.nz=5U;config.length_x_m=1.0;config.length_y_m=1.0;config.length_z_m=1.0;config.dt_s=5.0e-13;config.periodic_particles=false;config.field_boundary.mode=GridBoundaryMode3D::electric_wall;config.particle_boundary=ParticleWallMode::absorb;config.enable_secondary_yield_report=true;
    StaggeredElectromagneticPic3D solver(config);
    const std::size_t n=config.nx*config.ny*config.nz;std::vector<double> ones(n,1.0),zeros(n,0.0);
    solver.set_fields(ones,zeros,ones,zeros,ones,zeros);
    require(solver.electric_x().front()==0.0&&solver.electric_z().front()==0.0,"3-D electric-wall boundary clamps edge electric fields");
    PicParticle3D p;p.mass_kg=9.1093837139e-31;p.charge_c=-1.602176634e-19;p.position_m={0.999999,0.5,0.5};p.velocity_m_per_s={1.0e7,0.0,0.0};solver.set_particles({p});
    solver.step(1U);const auto diag=solver.diagnostics();
    require(diag.absorbed_particles>0U&&diag.particle_count==0U,"staggered 3-D PIC absorbing wall removes outgoing particles");
    require(diag.secondary_macro_weight>0.0,"staggered 3-D PIC reports secondary-emission macro-weight");

    config.field_boundary.mode=GridBoundaryMode3D::absorbing_sponge;config.periodic_particles=true;StaggeredElectromagneticPic3D damped(config);damped.set_fields(ones,ones,ones,ones,ones,ones);
    require(damped.electric_x().front()<1.0&&damped.magnetic_y().front()<1.0,"3-D absorbing sponge damps boundary fields");
}

void cfl_guard_rejects_large_step(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic3DConfig config;config.nx=4U;config.ny=4U;config.nz=4U;config.dt_s=1.0e-6;
    bool rejected=false;try{StaggeredElectromagneticPic3D bad(config);(void)bad;}catch(const std::exception&){rejected=true;}
    require(rejected,"staggered 3-D Yee EM-PIC rejects CFL-unsafe timestep");
}
}

int main(){
    try{
        staggered_vacuum_mode_is_bounded();
        moving_beam_reconstructs_current();
        wall_and_sponge_primitives_work();
        cfl_guard_rejects_large_step();
        std::cout<<"v0.9.9 staggered 3-D Yee EM-PIC tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.9.9 staggered PIC3D test failure: "<<error.what()<<'\n';
        return 1;
    }
}
