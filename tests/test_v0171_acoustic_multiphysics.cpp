#include "cfd/acoustics/coupling.hpp"
#include "cfd/acoustics/kspace2d.hpp"
#include "cfd/fem/mesh2d.hpp"
#include "cfd/multiphysics/bioheat.hpp"
#include "cfd/solvers/fem/dynamic_bar1d.hpp"
#include "cfd/solvers/fem/piezoelectric_bar1d.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool v,const char* msg){if(!v)throw std::runtime_error(msg);}


void test_time_reversal(){
    cfd::acoustics::KSpace2DConfig cfg;cfg.nx=32;cfg.ny=32;cfg.dx=2.5e-4;cfg.dy=2.5e-4;cfg.cfl=0.1;cfg.pml_cells=0;
    cfd::acoustics::KSpaceAcoustic2D forward(cfg);forward.set_uniform_medium(1500.0,1000.0);
    std::vector<std::size_t> sensor_cells;
    for(std::size_t i=0;i<cfg.nx;++i){sensor_cells.push_back(i);sensor_cells.push_back((cfg.ny-1U)*cfg.nx+i);}
    for(std::size_t j=1;j+1U<cfg.ny;++j){sensor_cells.push_back(j*cfg.nx);sensor_cells.push_back(j*cfg.nx+cfg.nx-1U);}
    forward.add_sensor_array(sensor_cells);
    forward.initialize_gaussian_pressure(0.5,0.5,5.0e-4,1.0e5);
    forward.run(260);
    std::vector<std::vector<double>> traces;traces.reserve(forward.sensors().size());
    for(const auto& s:forward.sensors())traces.push_back(s.pressure);
    cfd::acoustics::KSpaceAcoustic2D reverse(cfg);reverse.set_uniform_medium(1500.0,1000.0);
    const auto recon=reverse.time_reversal_reconstruct(sensor_cells,traces,true);
    const auto max_it=std::max_element(recon.begin(),recon.end());const std::size_t q=static_cast<std::size_t>(max_it-recon.begin());
    const std::size_t ix=q%cfg.nx,iy=q/cfg.nx;
    require(std::abs(static_cast<int>(ix)-16)<=3&&std::abs(static_cast<int>(iy)-16)<=3,"time reversal focuses near original source");
    require(*max_it>100.0,"time reversal reconstructs non-trivial initial pressure");
}

void test_structural_and_pennes_coupling(){
    auto mesh=cfd::fem::make_rectangle_tri_mesh(2,2,1.0,0.5);
    std::vector<double> edge_p(mesh.boundary_edges.size(),0.0);
    for(std::size_t e=0;e<mesh.boundary_edges.size();++e)if(mesh.boundary_edges[e].patch==1)edge_p[e]=2000.0;
    const auto nodal=cfd::acoustics::pressure_to_structural_boundary_load(mesh,edge_p,1);
    double fx=0.0,fy=0.0;for(const auto& f:nodal){fx+=f.x;fy+=f.y;}
    require(std::abs(fx+1000.0)<1e-10&&std::abs(fy)<1e-12,"pressure traction conserves integrated FEM load");

    cfd::fem::DynamicBar1DConfig dc;dc.elements=12;dc.length=0.1;dc.area=1.0e-4;dc.young_modulus=2.0e7;dc.density=1000.0;dc.dt=2.0e-5;
    cfd::fem::DynamicBar1D bar(dc);bar.initialize([](double){return 0.0;});
    std::vector<double> trace(40,5.0e3);cfd::acoustics::drive_dynamic_bar_from_pressure_trace(bar,trace,dc.area);
    require(std::abs(bar.displacement().back())>1e-10,"acoustic pressure drives structural vibration");

    cfd::multiphysics::PennesBioheat2DConfig bc;bc.nx=6;bc.ny=6;bc.thermal_conductivity_w_per_m_k=0.0;bc.blood_perfusion_per_s=0.0;bc.metabolic_heat_w_per_m3=0.0;bc.dt_s=1.0;bc.boundary=cfd::multiphysics::BioheatBoundary::periodic;
    cfd::multiphysics::PennesBioheat2D bio(bc);bio.initialize(310.0);
    std::vector<double> intensity(36,1000.0),alpha(36,2.0);cfd::acoustics::apply_acoustic_heating_to_pennes(bio,intensity,alpha);bio.step();
    const double expected=310.0+4000.0/(bc.tissue_density_kg_per_m3*bc.tissue_specific_heat_j_per_kg_k);
    require(std::abs(bio.temperature_k()[0]-expected)<1e-9,"ultrasound heating feeds Pennes source exactly");
}

void test_acoustic_cfd_coupling(){
    const std::vector<double> intensity(4,2000.0),alpha(4,1.5),rho(4,1000.0),c(4,1500.0);
    const auto acceleration=cfd::acoustics::acoustic_streaming_acceleration(intensity,alpha,rho,c,{1.0,0.0,0.0});
    require(acceleration[0].x>0.0&&acceleration[0].y==0.0,"acoustic radiation force direction");
    auto mesh=cfd::fvm::make_cartesian_hexa_mesh(4,1,1,1.0,0.1,0.1);
    cfd::fvm::CollocatedIncompressibleConfig cfg;cfg.kinematic_viscosity=0.0;cfg.dt=1.0e-2;cfg.include_convection=false;cfg.pressure_correctors=1;cfg.nonorthogonal_correctors=0;cfg.momentum_iterations=100;cfg.pressure_iterations=200;
    cfd::fvm::CollocatedIncompressible flow(std::move(mesh),cfg);flow.initialize_uniform();flow.set_body_acceleration(acceleration);flow.step_piso();
    double mean=0.0;for(const auto& u:flow.velocity())mean+=u.x;mean/=static_cast<double>(flow.velocity().size());
    require(mean>0.0,"acoustic streaming acceleration drives FVM velocity");
}

void test_piezoelectric_fem(){
    cfd::fem::PiezoelectricBar1DConfig cfg;cfg.elements=8;cfg.length_m=1.0e-3;cfg.area_m2=1.0e-4;cfg.elastic_modulus_pa=6.0e10;cfg.piezoelectric_stress_c_per_m2=15.0;cfg.permittivity_f_per_m=1.5e-8;
    cfd::fem::PiezoelectricBar1D piezo(cfg);
    const double pressure=1.0e5,force=cfd::acoustics::acoustic_pressure_to_bar_tip_force(pressure,cfg.area_m2);piezo.solve_receive(force);
    const double ceff=cfg.elastic_modulus_pa+cfg.piezoelectric_stress_c_per_m2*cfg.piezoelectric_stress_c_per_m2/cfg.permittivity_f_per_m;
    const double strain=force/(cfg.area_m2*ceff);const double expected_v=cfg.piezoelectric_stress_c_per_m2/cfg.permittivity_f_per_m*strain*cfg.length_m;
    require(std::abs(std::abs(piezo.open_circuit_voltage_v())-std::abs(expected_v))<1e-8*std::max(1.0,std::abs(expected_v)),"piezo FEM open-circuit acoustic receive voltage");
    piezo.solve_transmit(20.0);require(std::abs(piezo.tip_displacement_m())>1e-12,"piezo FEM voltage produces mechanical displacement");
}

void test_inverse_coupling(){
    const std::vector<double> base{0.0,1.0,0.5,-0.25,0.0};
    std::vector<std::vector<double>> measured(2,base);for(auto& t:measured)for(double& x:t)x*=2.0;
    cfd::acoustics::AcousticTraceSimulator simulator=[&](std::span<const double> p){std::vector<std::vector<double>> out(2,base);for(auto& t:out)for(double& x:t)x*=p[0];return out;};
    const auto fit=cfd::acoustics::fit_acoustic_parameters(simulator,measured,{0.5},{0.5},30,1e-7);
    require(std::abs(fit.parameters[0]-2.0)<1e-6&&fit.misfit<1e-8,"acoustic inverse parameter fit recovers synthetic amplitude");
}
}

int main(){
    try{test_time_reversal();test_structural_and_pennes_coupling();test_acoustic_cfd_coupling();test_piezoelectric_fem();test_inverse_coupling();}
    catch(const std::exception& e){std::cerr<<"v0.17.1 acoustic multiphysics regression failed: "<<e.what()<<'\n';return 1;}
    std::cout<<"v0.17.1 acoustic multiphysics regression passed\n";return 0;
}
