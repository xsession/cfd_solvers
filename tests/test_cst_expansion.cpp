#include "cfd/em/frequency_domain.hpp"
#include "cfd/em/edge_fem2d.hpp"
#include "cfd/fem/mesh2d.hpp"
#include "cfd/multiphysics/bioheat.hpp"
#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}

double speed(const cfd::particle::Vec3& v){return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z);}

void frequency_domain_maxwell(){
    using namespace cfd::em;
    FrequencyDomain1DConfig config;config.points=161U;config.length_m=1.0;config.frequency_hz=8.0e7;
    std::vector<Complex> current(config.points);
    for(std::size_t i=0;i<config.points;++i){const double x=config.length_m*static_cast<double>(i)/static_cast<double>(config.points-1U);current[i]=std::sin(std::numbers::pi*x/config.length_m);}
    const auto result=solve_pec_driven_maxwell_1d(config,current);
    constexpr double mu0=1.25663706212e-6,epsilon0=8.8541878128e-12;
    const double omega=2.0*std::numbers::pi*config.frequency_hz;
    const double k2=omega*omega*mu0*epsilon0;
    const double lambda=std::numbers::pi*std::numbers::pi/(config.length_m*config.length_m);
    const Complex exact=Complex{0.0,omega*mu0}/(k2-lambda);
    const Complex numerical=result.electric_v_per_m[config.points/2U];
    require(std::abs(numerical-exact)/std::abs(exact)<2.0e-4,"driven frequency-domain Maxwell manufactured sine response");
    require(result.electric_v_per_m.front()==Complex{}&&result.electric_v_per_m.back()==Complex{},"frequency-domain PEC endpoints");

    const auto modes=pec_cavity_eigenmodes_1d(101U,1.0,4.0,1.0,3U);
    constexpr double c0=299792458.0;
    const double exact_fundamental=c0/(2.0*std::sqrt(4.0));
    require(std::abs(modes[0].frequency_hz-exact_fundamental)/exact_fundamental<6.0e-5,"PEC cavity fundamental eigenfrequency");
    require(modes[0].frequency_hz<modes[1].frequency_hz&&modes[1].frequency_hz<modes[2].frequency_hz,"cavity modes sorted by frequency");
}



double edge_maxwell_error(std::size_t cells){
    using namespace cfd::em;
    const auto mesh=cfd::fem::make_rectangle_tri_mesh(cells,cells,1.0,0.8);
    EdgeMaxwell2DConfig config;config.frequency_hz=8.0e7;
    constexpr double mu0=1.25663706212e-6,epsilon0=8.8541878128e-12;
    const double omega=2.0*std::numbers::pi*config.frequency_hz;
    const double lambda=std::numbers::pi*std::numbers::pi;
    const double coefficient=lambda/(mu0*config.relative_permeability)-omega*omega*epsilon0*config.relative_permittivity;
    const auto source=[=](cfd::fem::Node2 point)->ComplexVec2{
        const double ey=std::sin(std::numbers::pi*point.x);
        return {Complex{},Complex{0.0,coefficient*ey/omega}};
    };
    const auto solution=solve_driven_edge_maxwell_2d(mesh,config,source);
    double error2=0.0,reference2=0.0;
    for(std::size_t element=0;element<mesh.triangles.size();++element){
        const auto& tri=mesh.triangles[element];double x=0.0,y=0.0;
        for(std::size_t k=0;k<3U;++k){x+=mesh.nodes[tri.node[k]].x/3.0;y+=mesh.nodes[tri.node[k]].y/3.0;}
        (void)y;
        const ComplexVec2 exact{Complex{},Complex{std::sin(std::numbers::pi*x),0.0}};
        const auto numerical=solution.electric_centroid_v_per_m[element];
        error2+=std::norm(numerical[0]-exact[0])+std::norm(numerical[1]-exact[1]);
        reference2+=std::norm(exact[0])+std::norm(exact[1]);
    }
    return std::sqrt(error2/reference2);
}

void edge_element_maxwell(){
    const double coarse=edge_maxwell_error(4U),fine=edge_maxwell_error(8U);
    require(coarse<0.22,"Nedelec edge-Maxwell manufactured coarse accuracy");
    require(fine<coarse*0.65,"Nedelec edge-Maxwell mesh refinement");
}

void charged_particles_and_pic(){
    using namespace cfd::particle;
    ChargedParticle particle;particle.velocity_m_per_s={1.0,0.25,-0.1};particle.charge_c=1.0;particle.mass_kg=1.0;
    const double initial_speed=speed(particle.velocity_m_per_s);
    const ElectromagneticField magnetic{{0.0,0.0,0.0},{0.0,0.0,1.0}};
    for(std::size_t i=0;i<20000U;++i)boris_push(particle,magnetic,5.0e-4);
    require(std::abs(speed(particle.velocity_m_per_s)-initial_speed)<2.0e-12,"Boris pusher magnetic-field energy conservation");

    constexpr std::size_t n=64U;constexpr double length=1.0,rho0=2.0e-10,epsilon0=8.8541878128e-12;
    std::vector<double> rho(n);for(std::size_t i=0;i<n;++i)rho[i]=rho0*std::cos(2.0*std::numbers::pi*static_cast<double>(i)/static_cast<double>(n));
    const auto field=periodic_electric_field_from_charge_density(rho,length,1.0,false);
    const double amplitude=rho0/(epsilon0*2.0*std::numbers::pi/length);double max_error=0.0;
    for(std::size_t i=0;i<n;++i){const double exact=amplitude*std::sin(2.0*std::numbers::pi*static_cast<double>(i)/static_cast<double>(n));max_error=std::max(max_error,std::abs(field[i]-exact));}
    require(max_error/std::abs(amplitude)<2.0e-12,"spectral periodic PIC Poisson field");

    ElectrostaticPic1DConfig pic_config;pic_config.grid_points=n;pic_config.length_m=length;pic_config.dt_s=1.0e-3;
    ElectrostaticPic1D pic(pic_config);std::vector<PicParticle1D> particles;particles.reserve(n);
    for(std::size_t i=0;i<n;++i)particles.push_back({length*static_cast<double>(i)/static_cast<double>(n),0.5,1.0e-12,1.0,1.0});
    pic.set_particles(particles);pic.deposit_and_solve();
    const double peak_field=*std::max_element(pic.electric_field().begin(),pic.electric_field().end(),[](double a,double b){return std::abs(a)<std::abs(b);});
    require(std::abs(peak_field)<1.0e-8,"neutralized uniform PIC distribution has zero field");
    pic.step(3U);for(const auto& p:pic.particles())require(std::abs(p.velocity_m_per_s-0.5)<1.0e-12,"uniform PIC drift preserves particle velocity");
}

void bioheat(){
    using namespace cfd::multiphysics;
    const double sar=sar_from_rms_electric_field(1.0,1000.0,100.0);
    require(std::abs(sar-10.0)<1.0e-13,"SAR from RMS electric field");
    PennesBioheat2DConfig config;config.nx=9U;config.ny=7U;config.width_m=0.08;config.height_m=0.06;config.dt_s=1.0;
    config.tissue_density_kg_per_m3=1000.0;config.tissue_specific_heat_j_per_kg_k=3600.0;config.thermal_conductivity_w_per_m_k=0.5;
    config.blood_density_kg_per_m3=1060.0;config.blood_specific_heat_j_per_kg_k=3770.0;config.blood_perfusion_per_s=0.01;config.blood_temperature_k=310.0;config.boundary=BioheatBoundary::periodic;
    PennesBioheat2D solver(config);solver.initialize(300.0);solver.set_sar(sar);solver.step();
    const double capacity=config.tissue_density_kg_per_m3*config.tissue_specific_heat_j_per_kg_k;
    const double perfusion=config.blood_density_kg_per_m3*config.blood_specific_heat_j_per_kg_k*config.blood_perfusion_per_s;
    const double expected=(capacity/config.dt_s*300.0+perfusion*config.blood_temperature_k+config.tissue_density_kg_per_m3*sar)/(capacity/config.dt_s+perfusion);
    require(solver.linear_result().converged,"Pennes bioheat implicit solve converges");
    for(double t:solver.temperature_k())require(std::abs(t-expected)<2.0e-10,"uniform SAR/perfusion Pennes analytical step");

    config.boundary=BioheatBoundary::fixed_temperature;config.fixed_boundary_temperature_k=295.0;config.blood_perfusion_per_s=0.0;
    PennesBioheat2D fixed(config);fixed.initialize(300.0);fixed.set_sar(5.0);fixed.step();
    for(std::size_t x=0;x<config.nx;++x){require(std::abs(fixed.temperature_k()[x]-295.0)<1.0e-12,"fixed bioheat lower boundary");require(std::abs(fixed.temperature_k()[(config.ny-1U)*config.nx+x]-295.0)<1.0e-12,"fixed bioheat upper boundary");}
}
}

int main(){
    try{frequency_domain_maxwell();edge_element_maxwell();charged_particles_and_pic();bioheat();std::cout<<"CST expansion tests passed\n";return 0;}
    catch(const std::exception& error){std::cerr<<"CST expansion test failure: "<<error.what()<<'\n';return 1;}
}
