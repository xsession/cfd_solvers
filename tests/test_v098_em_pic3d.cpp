#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

void vacuum_mode_remains_bounded(){
    using namespace cfd::particle;
    ElectromagneticPic3DConfig config;config.nx=10U;config.ny=8U;config.nz=6U;config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;config.solve_longitudinal_poisson=false;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);const double dz=config.length_z_m/static_cast<double>(config.nz);
    config.dt_s=0.035/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
    ElectromagneticPic3D solver(config);solver.initialize_z_polarized_mode(0.75,1U,0U,0U);
    const double energy0=solver.diagnostics().field_energy_j;require(energy0>0.0,"3-D EM-PIC initialized vacuum mode has field energy");
    const auto sample=solver.gather_field({0.23,0.17,0.11});
    require(std::isfinite(sample.electric_v_per_m.z)&&std::isfinite(sample.magnetic_t.y),"3-D EM-PIC trilinear E/B gather is finite");
    solver.step(6U);const auto diag=solver.diagnostics();
    require(diag.field_energy_j>0.1*energy0&&diag.field_energy_j<4.0*energy0,"3-D EM-PIC vacuum mode remains bounded for compact timestep");
    require(diag.particle_count==0U,"empty 3-D EM-PIC vacuum mode keeps zero particles");
}

void beam_motion_deposits_charge_conserving_current(){
    using namespace cfd::particle;
    ElectromagneticPic3DConfig config;config.nx=8U;config.ny=7U;config.nz=6U;config.length_x_m=1.0;config.length_y_m=0.7;config.length_z_m=0.6;config.solve_longitudinal_poisson=false;config.dt_s=1.0e-12;
    ElectromagneticPic3D solver(config);std::vector<PicParticle3D> particles;
    for(std::size_t i=0;i<6U;++i){PicParticle3D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={0.07+0.13*static_cast<double>(i),0.09+0.07*static_cast<double>(i),0.05+0.045*static_cast<double>(i)};p.velocity_m_per_s={2.0e4,-1.0e4,3.0e4};particles.push_back(p);} 
    solver.set_particles(std::move(particles));solver.step(1U);const auto diag=solver.diagnostics();
    double max_current=0.0;for(double value:solver.current_x())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_y())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_z())max_current=std::max(max_current,std::abs(value));
    require(diag.particle_count==6U,"3-D EM-PIC keeps periodic beam particles");
    require(max_current>0.0,"3-D EM-PIC deposits nonzero motion current");
    require(diag.charge_continuity_linf_residual<1.0e-6,"3-D EM-PIC spectral current satisfies continuity");
}

void poisson_correction_and_collision_coupling(){
    using namespace cfd::particle;
    ElectromagneticPic3DConfig config;config.nx=7U;config.ny=6U;config.nz=5U;config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;config.solve_longitudinal_poisson=true;config.dt_s=5.0e-13;
    ElectromagneticPic3D solver(config);PicParticle3D p;p.mass_kg=9.1093837139e-31;p.charge_c=-1.602176634e-19;p.position_m={0.31,0.27,0.19};p.velocity_m_per_s={2.0e5,1.0e5,-1.0e5};solver.set_particles({p});
    solver.deposit_sources();double peak_e=0.0;for(double value:solver.electric_x())peak_e=std::max(peak_e,std::abs(value));for(double value:solver.electric_y())peak_e=std::max(peak_e,std::abs(value));for(double value:solver.electric_z())peak_e=std::max(peak_e,std::abs(value));
    require(peak_e>0.0,"3-D EM-PIC Poisson correction produces longitudinal field");
    NeutralCollisionModel gas;gas.neutral_density_per_m3=1.0e31;gas.elastic_cross_section_m2=1.0e-18;gas.random_seed=123U;solver.step(1U,&gas);const auto diag=solver.diagnostics();
    require(diag.collision_statistics.elastic_events>0U,"3-D EM-PIC couples Monte-Carlo neutral collisions");
    require(diag.particle_count==1U,"elastic collision leaves particle count unchanged");
}

void rejected_bad_timestep(){
    using namespace cfd::particle;
    ElectromagneticPic3DConfig config;config.nx=4U;config.ny=4U;config.nz=4U;config.dt_s=1.0e-6;
    bool rejected=false;try{ElectromagneticPic3D bad(config);(void)bad;}catch(const std::exception&){rejected=true;}
    require(rejected,"3-D EM-PIC rejects CFL-unsafe timestep");
}
}

int main(){
    try{
        vacuum_mode_remains_bounded();
        beam_motion_deposits_charge_conserving_current();
        poisson_correction_and_collision_coupling();
        rejected_bad_timestep();
        std::cout<<"v0.9.8 3-D EM-PIC baseline tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.9.8 EM-PIC3D test failure: "<<error.what()<<'\n';
        return 1;
    }
}
