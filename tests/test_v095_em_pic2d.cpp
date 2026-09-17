#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

double max_abs(const std::vector<double>& values){double out=0.0;for(double v:values)out=std::max(out,std::abs(v));return out;}

void vacuum_tm_wave_is_stable(){
    using namespace cfd::particle;
    ElectromagneticPic2DConfig config;config.nx=28U;config.ny=18U;config.length_x_m=1.0;config.length_y_m=0.75;config.solve_longitudinal_poisson=false;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);
    config.dt_s=0.05*std::min(dx,dy)/299792458.0;
    ElectromagneticPic2D solver(config);solver.initialize_tm_z_mode(1.0,1U,0U);
    const double initial=solver.diagnostics().field_energy_j;
    require(initial>0.0,"2-D EM-PIC TM mode has finite initial field energy");
    solver.step(12U);
    const auto diag=solver.diagnostics();
    require(diag.particle_count==0U,"vacuum 2-D EM-PIC keeps empty particle set empty");
    require(std::isfinite(diag.field_energy_j)&&diag.field_energy_j>0.70*initial&&diag.field_energy_j<1.30*initial,"short 2-D EM-PIC vacuum TM wave keeps bounded energy");
}

void current_deposition_and_longitudinal_continuity(){
    using namespace cfd::particle;
    ElectromagneticPic2DConfig config;config.nx=20U;config.ny=16U;config.length_x_m=1.0;config.length_y_m=0.8;config.dt_s=1.0e-9;config.solve_longitudinal_poisson=false;
    ElectromagneticPic2D solver(config);std::vector<PicParticle2D> particles;
    for(std::size_t i=0;i<6U;++i){PicParticle2D p;p.mass_kg=1.0;p.charge_c=(i%2U==0U?1.0:-0.8)*1.0e-15;p.weight=1.0;p.position_m={0.11+0.12*static_cast<double>(i),0.13+0.07*static_cast<double>(i)};p.velocity_m_per_s={1.5e4-700.0*static_cast<double>(i),-8.0e3+500.0*static_cast<double>(i),0.0};particles.push_back(p);}    
    solver.set_particles(std::move(particles));solver.step(1U);const auto diag=solver.diagnostics();
    require(diag.charge_continuity_linf_residual<1.0e-8,"2-D/3-V EM-PIC in-plane current satisfies spectral continuity");
    require(max_abs(solver.current_x())>0.0&&max_abs(solver.current_y())>0.0,"2-D/3-V EM-PIC deposits in-plane currents");
    require(std::isfinite(diag.total_energy_j),"2-D/3-V EM-PIC diagnostics remain finite after current-coupled step");
}

void transverse_current_couples_to_ez(){
    using namespace cfd::particle;
    ElectromagneticPic2DConfig config;config.nx=18U;config.ny=14U;config.length_x_m=1.0;config.length_y_m=0.7;config.dt_s=2.0e-12;config.solve_longitudinal_poisson=false;
    ElectromagneticPic2D solver(config);std::vector<PicParticle2D> particles(2U);
    particles[0].mass_kg=1.0e-3;particles[0].charge_c=1.0e-11;particles[0].position_m={0.23,0.31};particles[0].velocity_m_per_s={0.0,0.0,2.0e5};
    particles[1].mass_kg=1.0e-3;particles[1].charge_c=-1.0e-11;particles[1].position_m={0.71,0.47};particles[1].velocity_m_per_s={0.0,0.0,-1.0e5};
    solver.set_particles(std::move(particles));solver.step(1U);
    require(max_abs(solver.current_z())>0.0,"2-D/3-V EM-PIC deposits transverse Jz");
    require(max_abs(solver.electric_z())>0.0,"2-D/3-V EM-PIC transverse current couples into Ez");
}

void field_gather_and_collision_coupling(){
    using namespace cfd::particle;
    ElectromagneticPic2DConfig config;config.nx=12U;config.ny=10U;config.dt_s=1.0e-12;config.solve_longitudinal_poisson=false;
    ElectromagneticPic2D solver(config);solver.initialize_tm_z_mode(0.5,1U,1U);
    PicParticle2D p;p.mass_kg=9.1093837139e-31;p.charge_c=-1.602176634e-19;p.position_m={0.37,0.29};p.velocity_m_per_s={1.0e6,2.0e5,1.0e5};solver.set_particles({p});
    const auto field=solver.gather_field(p.position_m);require(std::isfinite(field.electric_v_per_m.z)&&std::isfinite(field.magnetic_t.x),"2-D/3-V EM-PIC bilinear field gather returns finite values");
    NeutralCollisionModel gas;gas.neutral_density_per_m3=1.0e32;gas.elastic_cross_section_m2=1.0e-18;gas.ionization_cross_section_m2=0.0;gas.random_seed=123U;
    solver.step(2U,&gas);const auto diag=solver.diagnostics();
    require(diag.collision_statistics.elastic_events>0U,"2-D/3-V EM-PIC invokes Monte-Carlo collision coupling");
}
}

int main(){
    try{
        vacuum_tm_wave_is_stable();
        current_deposition_and_longitudinal_continuity();
        transverse_current_couples_to_ez();
        field_gather_and_collision_coupling();
        std::cout<<"v0.9.5 2-D/3-V EM-PIC tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.9.5 EM-PIC2D test failure: "<<error.what()<<'\n';
        return 1;
    }
}
