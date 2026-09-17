#include "cfd/em/system.hpp"
#include "cfd/particle/electromagnetic.hpp"
#include "cfd/multiphysics/bioheat.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 
double speed(const cfd::particle::Vec3& v){return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z);} 

void signal_integrity_and_emc(){
    using namespace cfd::em;
    constexpr std::size_t sps=16U;const std::array<int,10> bits{{0,1,0,1,1,0,0,1,0,1}};std::vector<double> waveform(bits.size()*sps);
    for(std::size_t k=0;k<bits.size();++k){
        for(std::size_t phase=0;phase<sps;++phase){
            const double ui=(static_cast<double>(phase)+0.5)/static_cast<double>(sps);
            const double level=bits[k]?1.0:0.0;
            waveform[k*sps+phase]=level+0.015*std::sin(2.0*std::numbers::pi*ui)+0.005*std::cos(static_cast<double>(k));
        }
    }
    const auto eye=analyze_nrz_eye(waveform,bits,sps);
    require(eye.symbol_count==bits.size()&&eye.one_mean>0.98&&eye.zero_mean<0.02,"eye analysis separates one/zero levels");
    require(eye.eye_height>0.94&&eye.eye_width_ui>0.8&&eye.estimated_ber<1.0e-100,"eye height/width/BER metrics are open for clean NRZ data");

    const auto esd=sample_double_exponential_pulse(100.0e-9,0.1e-9,30.0,0.8e-9,20.0e-9);
    const auto esd_metrics=probe_waveform(esd,0.1e-9);
    require(std::abs(esd_metrics.peak-30.0)<0.2,"double-exponential ESD waveform is peak-normalized");
    require(esd_metrics.impulse>0.0&&esd_metrics.energy>0.0,"EMC probe reports positive impulse/energy");
    const auto bci=sample_damped_sine(2.0e-6,1.0e-9,5.0,10.0e6,0.5e-6);
    const auto bci_metrics=probe_waveform(bci,1.0e-9);
    require(bci_metrics.peak>4.0&&bci_metrics.rms>0.2,"damped-sine BCI/lightning-style waveform is sampled and probed");
}

void pdn_tools(){
    using namespace cfd::em;using C=std::complex<double>;
    std::array<double,5> frequencies{{1.0e5,3.0e5,1.0e6,3.0e6,1.0e7}};std::array<double,5> target{{0.05,0.05,0.05,0.05,0.05}};
    std::array<DecouplingCapacitor,3> library{{{1.0e-6,0.02,0.8e-9,1U},{100.0e-9,0.015,0.4e-9,1U},{10.0e-9,0.01,0.25e-9,1U}}};
    const auto before=std::abs(pdn_parallel_impedance_ohm({},1.0e6,C{0.2,0.0}));
    const auto optimized=optimize_decoupling_greedy(frequencies,target,library,3U,C{0.2,0.0});
    require(!optimized.selected.empty()&&optimized.worst_ratio_after<optimized.worst_ratio_before,"greedy decoupling selection improves target-impedance margin");
    require(std::abs(pdn_parallel_impedance_ohm(optimized.selected,1.0e6,C{0.2,0.0}))<before,"selected decoupling lowers PDN impedance near resonance");

    PdnGridConfig grid;grid.nx=4U;grid.ny=3U;grid.horizontal_resistance_ohm=0.01;grid.vertical_resistance_ohm=0.02;std::vector<double> loads(grid.nx*grid.ny,0.0);loads[grid.nx*grid.ny-1U]=1.0;
    const std::array<PdnVoltageSource,2> sources{{{0U,1.0},{grid.nx,1.0}}};
    const auto ir=solve_pdn_ir_drop_grid(grid,loads,sources);
    require(ir.converged&&ir.minimum_voltage_v<1.0&&ir.maximum_drop_v>0.0,"PDN IR-drop grid solves loaded rail below source voltage");
    require(ir.voltage_v.back()==ir.minimum_voltage_v,"far loaded corner is the worst node in this regression grid");
}

void particles_and_collisions(){
    using namespace cfd::particle;constexpr double c0=299792458.0;
    ChargedParticle particle;particle.velocity_m_per_s={0.95*c0,0.02*c0,0.0};particle.charge_c=1.0;particle.mass_kg=1.0;const double gamma0=lorentz_gamma(particle);
    for(std::size_t i=0;i<5000U;++i)vay_push(particle,{{0,0,0},{0,0,0.5}},2.0e-4);
    require(std::abs(lorentz_gamma(particle)-gamma0)<5.0e-10,"Vay pusher preserves relativistic magnetic-only energy");

    std::vector<ChargedParticle> swarm(3U);for(auto& p:swarm){p.mass_kg=1.0;p.charge_c=1.0;p.velocity_m_per_s={2.0,0.0,0.0};}
    NeutralCollisionModel elastic;elastic.neutral_density_per_m3=1.0e20;elastic.elastic_cross_section_m2=1.0e-12;elastic.random_seed=1234ULL;
    const auto stats=apply_monte_carlo_collisions(swarm,1.0,elastic);
    require(stats.elastic_events==3U&&stats.ionization_events==0U,"Monte-Carlo neutral collision model scatters each high-probability elastic particle");
    for(const auto& p:swarm)require(std::abs(speed(p.velocity_m_per_s)-2.0)<1.0e-12,"elastic Monte-Carlo collisions preserve kinetic energy");

    std::vector<ChargedParticle> ionizing(1U);ionizing[0].mass_kg=9.1093837139e-31;ionizing[0].charge_c=-1.602176634e-19;ionizing[0].velocity_m_per_s={3.0e6,0.0,0.0};
    NeutralCollisionModel ion;ion.neutral_density_per_m3=1.0e22;ion.ionization_cross_section_m2=1.0e-12;ion.ionization_energy_ev=5.0;ion.random_seed=99ULL;
    const double initial_energy=kinetic_energy_j(ionizing[0]);const auto ion_stats=apply_monte_carlo_collisions(ionizing,1.0e-6,ion);
    require(ion_stats.ionization_events==1U&&ionizing.size()==2U&&kinetic_energy_j(ionizing[0])<initial_energy,"ionizing Monte-Carlo collision removes threshold energy and emits a secondary macro-particle");
}

void heterogeneous_bioheat(){
    using namespace cfd::multiphysics;
    VoxelTissueGrid3D grid;grid.nx=2U;grid.ny=2U;grid.nz=2U;grid.dx_m=0.01;grid.dy_m=0.01;grid.dz_m=0.01;
    grid.materials={VoxelTissueMaterial{1000.0,3600.0,0.5,1.0,0.01,0.0},VoxelTissueMaterial{500.0,2500.0,0.2,0.2,0.0,0.0}};
    grid.material_index={0,0,0,0,1,1,1,1};grid.electric_rms_v_per_m.assign(8U,100.0);
    const auto sar=voxel_sar_w_per_kg(grid);
    require(std::abs(sar[0]-10.0)<1.0e-12&&std::abs(sar[4]-4.0)<1.0e-12,"heterogeneous voxel SAR uses each tissue conductivity/density");
    const double max1g=max_mass_averaged_sar_w_per_kg(grid,sar,1.0e-3);
    require(std::abs(max1g-10.0)<1.0e-12,"1 g local SAR picks the hottest high-density voxel fraction");
    const auto projected=project_voxel_sar_to_pennes2d(grid,sar,2U,2U);
    for(double value:projected)require(value>7.9&&value<8.1,"mass-weighted voxel SAR projection averages stacked tissues into Pennes mesh");
    PennesBioheat2DConfig config;config.nx=3U;config.ny=3U;config.boundary=BioheatBoundary::periodic;PennesBioheat2D bio(config);bio.initialize(300.0);bio.set_sar(8.0);bio.step();
    require(bio.linear_result().converged&&bio.temperature_k()[4]>300.0,"projected SAR can drive the implicit Pennes solver");
}
}

int main(){
    try{signal_integrity_and_emc();pdn_tools();particles_and_collisions();heterogeneous_bioheat();std::cout<<"v0.9.2 CST workflow tests passed\n";return 0;}
    catch(const std::exception& error){std::cerr<<"v0.9.2 test failure: "<<error.what()<<'\n';return 1;}
}
