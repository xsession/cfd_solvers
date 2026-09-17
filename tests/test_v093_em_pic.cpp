#include "cfd/particle/electromagnetic.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

double mean(const std::vector<double>& values){return std::accumulate(values.begin(),values.end(),0.0)/static_cast<double>(values.size());}

double maximum_abs(const std::vector<double>& values){double result=0.0;for(double value:values)result=std::max(result,std::abs(value));return result;}

void charge_conserving_deposition(){
    using namespace cfd::particle;
    std::vector<ChargedParticle> before(3U),after(3U);
    for(std::size_t i=0;i<before.size();++i){
        before[i].mass_kg=1.0;before[i].charge_c=(i==2U?-0.5:1.0)*1.0e-12;before[i].weight=1.0;
        before[i].position_m={0.12+0.21*static_cast<double>(i),0.0,0.0};before[i].velocity_m_per_s={0.0,0.0,0.0};
        after[i]=before[i];after[i].position_m.x+=0.035+0.01*static_cast<double>(i);
    }
    const auto deposited=deposit_charge_conserving_current_1d(before,after,64U,1.0,1.0e-9);
    require(deposited.charge_density_old_c_per_m3.size()==64U&&deposited.current_x_a_per_m2.size()==64U,"charge-conserving deposition returns grid-sized arrays");
    require(deposited.continuity_linf_residual<1.0e-12,"spectral current reconstruction satisfies charge continuity");
}

void vacuum_yee_wave(){
    using namespace cfd::particle;
    ElectromagneticPic1DConfig config;config.grid_points=128U;config.length_m=1.0;config.solve_longitudinal_poisson=false;
    const double c0=299792458.0;config.dt_s=0.35*(config.length_m/static_cast<double>(config.grid_points))/c0;
    ElectromagneticPic1D solver(config);solver.initialize_right_traveling_mode(1.0,1U);
    const double initial_energy=solver.diagnostics().field_energy_j;
    solver.step(160U);
    const double final_energy=solver.diagnostics().field_energy_j;
    require(initial_energy>0.0&&std::abs(final_energy-initial_energy)/initial_energy<0.03,"source-free Yee electromagnetic PIC field update conserves a periodic wave energy envelope");
}

void transverse_beam_current_drives_field(){
    using namespace cfd::particle;
    ElectromagneticPic1DConfig config;config.grid_points=8U;config.length_m=1.0;config.solve_longitudinal_poisson=false;config.dt_s=1.0e-10;
    ElectromagneticPic1D solver(config);std::vector<ChargedParticle> particles(config.grid_points);
    for(std::size_t i=0;i<particles.size();++i){particles[i].mass_kg=1.0e-10;particles[i].charge_c=1.0e-15;particles[i].weight=1.0;particles[i].position_m={static_cast<double>(i)/static_cast<double>(config.grid_points),0.0,0.0};particles[i].velocity_m_per_s={0.0,1.0e5,0.0};}
    solver.set_particles(particles);solver.step();
    const double expected_jy=static_cast<double>(particles.size())*particles.front().charge_c*particles.front().velocity_m_per_s.y/config.length_m;
    constexpr double epsilon0=8.8541878128e-12;const double expected_ey=-config.dt_s*expected_jy/epsilon0;
    require(std::abs(mean(solver.current_y())-expected_jy)<std::abs(expected_jy)*1.0e-12,"uniform transverse macro-particle beam deposits the expected current density");
    require(std::abs(mean(solver.electric_y())-expected_ey)<std::abs(expected_ey)*1.0e-10+1.0e-18,"Ampere update responds to deposited transverse beam current");
    require(maximum_abs(solver.magnetic_z())<1.0e-18,"spatially uniform current does not create a spurious curl-generated magnetic field");
}

void collision_source_couples_to_em_pic(){
    using namespace cfd::particle;
    ElectromagneticPic1DConfig config;config.grid_points=16U;config.length_m=1.0;config.solve_longitudinal_poisson=false;config.dt_s=1.0e-12;
    ElectromagneticPic1D solver(config);std::vector<ChargedParticle> particles(1U);
    particles[0].mass_kg=9.1093837139e-31;particles[0].charge_c=-1.602176634e-19;particles[0].weight=1.0;particles[0].position_m={0.4,0.0,0.0};particles[0].velocity_m_per_s={3.0e6,0.0,0.0};
    solver.set_particles(particles);
    NeutralCollisionModel model;model.neutral_density_per_m3=1.0e30;model.ionization_cross_section_m2=1.0e-12;model.ionization_energy_ev=5.0;model.random_seed=7ULL;
    solver.step(1U,&model);const auto diagnostics=solver.diagnostics();
    require(diagnostics.collision_statistics.ionization_events==1U&&diagnostics.particle_count==2U,"electromagnetic PIC step couples MCC ionization into the particle population");
    solver.step();
    require(solver.particles().size()==2U,"secondary macro-particles continue to participate in later EM-PIC steps");
}
}

int main(){
    try{
        charge_conserving_deposition();
        vacuum_yee_wave();
        transverse_beam_current_drives_field();
        collision_source_couples_to_em_pic();
        std::cout<<"v0.9.3 electromagnetic PIC tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.9.3 EM-PIC test failure: "<<error.what()<<'\n';
        return 1;
    }
}
