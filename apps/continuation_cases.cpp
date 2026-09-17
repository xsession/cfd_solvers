#include "continuation_cases.hpp"
#include <numbers>
#include "cfd/em/frequency_domain.hpp"
#include "cfd/em/edge_fem2d.hpp"
#include "cfd/em/edge_fem3d.hpp"
#include "cfd/em/wave_port.hpp"
#include "cfd/fem/mesh3d.hpp"
#include "cfd/particle/electromagnetic.hpp"
#include "cfd/particle/plasma.hpp"
#include "cfd/particle/transport.hpp"
#include "cfd/multiphysics/bioheat.hpp"
#include "cfd/fem/mesh2d.hpp"
#include "cfd/circuit/analysis.hpp"
#include "cfd/circuit/spice.hpp"
#include "cfd/rf/network.hpp"
#include "cfd/rf/antenna.hpp"
#include "cfd/rf/thin_wire_mom.hpp"
#include "cfd/chemistry/aqueous_equilibrium.hpp"
#include "cfd/chemistry/implicit_reactor.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/solvers/fvm/scalar_transport.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fem/elasticity2d.hpp"
#include "cfd/solvers/lbm/one_step_pull.hpp"
#include "cfd/solvers/optics/gaussian_beam.hpp"
#include "cfd/multiphysics/electro_thermal.hpp"
#include "cfd/workflow/campaign.hpp"
#include "cfd/workflow/campaign_control.hpp"
#include "cfd/workflow/deploy.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}

int trt(){
    cfd::lbm::D3Q19PullSolver solver({24,24,24,0.8F},1.1F);
    solver.initialize_taylor_green();const double mass=solver.mass();const auto start=Clock::now();solver.step(30);
    const double seconds=elapsed(start),drift=std::abs(solver.mass()-mass)/mass;
    std::cout<<"case=lbm-trt lattice=d3q19 steps=30 seconds="<<seconds<<" relative_mass_drift="<<drift
        <<" population_bytes="<<solver.population_bytes()<<'\n';
    return drift<1.0e-5?0:1;
}
int fdtd_pml(){
    cfd::fdtd::Maxwell1DConfig config{400,1.0e-4,0.7,1.0,1.0,cfd::fdtd::Boundary1D::pml};config.pml_cells=32;
    cfd::fdtd::Maxwell1D solver(config);solver.initialize_gaussian(0.5,0.025);const double energy=solver.energy();
    const auto start=Clock::now();solver.step(800);const double ratio=solver.energy()/energy;
    std::cout<<"case=fdtd-pml1d steps=800 seconds="<<elapsed(start)<<" residual_field_energy_ratio="<<ratio<<'\n';
    return ratio<1.0e-4?0:1;
}
int reactor(){
    using namespace cfd::chemistry;
    ReactionNetwork network({{"A",0.01,0},{"B",0.01,0}});
    network.add_reaction({{{0,1.0}},{{1,1.0}},{1000.0,0.0,0.0},EquilibriumConstant{4.0,298.15,0.0}});
    std::vector<double> concentrations{1.0,0.0};ReactorConfig config;config.initial_step=0.01;
    const auto start=Clock::now();const auto result=integrate_isothermal(network,concentrations,298.15,0.05,config);
    std::cout<<"case=chemistry-reactor seconds="<<elapsed(start)<<" A="<<concentrations[0]<<" B="<<concentrations[1]
        <<" accepted_steps="<<result.accepted_steps<<" rejected_steps="<<result.rejected_steps<<'\n';
    return std::abs(concentrations[0]-0.2)<1.0e-5?0:1;
}
int equilibrium(){
    const std::vector<cfd::chemistry::AcidFamily> families{{0.1,0,{1.0e-5}}};
    const auto start=Clock::now();const auto result=cfd::chemistry::equilibrate_acids(families,0.05);
    std::cout<<"case=chemistry-equilibrium seconds="<<elapsed(start)<<" pH="<<result.ph
        <<" charge_residual_mol_per_litre="<<result.charge_residual_mol_per_litre<<'\n';
    return std::abs(result.ph-5.0)<2.0e-4?0:1;
}
int gaussian(){
    const auto start=Clock::now();cfd::optics::GaussianBeam beam(500.0e-9,1.0e-3);
    beam.thin_lens(0.1);const double focus=-beam.q().real();beam.propagate(focus);
    std::cout<<"case=optics-gaussian seconds="<<elapsed(start)<<" focus_m="<<focus<<" waist_radius_m="<<beam.radius()<<'\n';
    return beam.radius()>0.0&&beam.radius()<1.0e-4?0:1;
}
int thermal_structure(){
    using Type=cfd::fem::ScalarBoundaryType;
    cfd::multiphysics::JouleHeatingCoupler2D coupled(cfd::fem::make_rectangle_tri_mesh(12,8));
    coupled.set_electrical_boundary(0,{Type::dirichlet,[](auto){return 0.0;},{}});
    coupled.set_electrical_boundary(1,{Type::dirichlet,[](auto){return 1.0;},{}});
    coupled.initialize_temperature([](auto){return 300.0;});
    coupled.set_thermal_dirichlet([](auto,double t){return 300.0+t;});
    const auto start=Clock::now();coupled.solve_electrical([](auto){return 1.0;});coupled.thermal_run(10);
    cfd::fem::Elasticity2D solid(coupled.thermal().mesh());
    solid.set_dirichlet([](auto x){return x.x==0.0;},[](auto){return cfd::fem::Displacement2{};},true,false);
    solid.set_dirichlet([](auto x){return x.y==0.0;},[](auto){return cfd::fem::Displacement2{};},false,true);
    solid.solve_thermal(coupled.thermal().temperature(),1.0e-3,300.0);
    double error=0.0;
    for(std::size_t i=0;i<solid.mesh().node_count();++i)
        error=std::max(error,std::abs(solid.displacement()[i].x-1.0e-3*coupled.thermal().time()*solid.mesh().nodes[i].x));
    std::cout<<"case=multiphysics-thermoelastic seconds="<<elapsed(start)<<" displacement_max_error="<<error<<'\n';
    return error<1.0e-8?0:1;
}
int fvm_workspace(){
    using namespace cfd::fvm;
    const auto mesh=make_cartesian_hexa_mesh(32,24,8);
    std::vector<double> values(mesh.cell_count()),boundary(mesh.face_count()),flux(mesh.face_count()),out(mesh.cell_count());
    const auto value=[](Vec3 x){return 1.0+std::sin(4.0*x.x)*std::cos(2.0*x.y);};
    for(std::size_t c=0;c<values.size();++c)values[c]=value(mesh.cells()[c].center);
    for(std::size_t f=0;f<flux.size();++f){boundary[f]=value(mesh.faces()[f].center);flux[f]=mesh.faces()[f].area.x;}
    constexpr auto scheme=FaceInterpolationScheme::bounded_linear;
    ScalarSchemeWorkspace workspace;workspace.resize(mesh);
    convective_divergence_scalar_into(mesh,values,flux,out,workspace,scheme,boundary);
    auto start=Clock::now();
    for(int repeat=0;repeat<100;++repeat)out=convective_divergence_scalar(mesh,values,flux,scheme,boundary);
    const double reference_seconds=elapsed(start);const auto reference=out;start=Clock::now();
    for(int repeat=0;repeat<100;++repeat)convective_divergence_scalar_into(mesh,values,flux,out,workspace,scheme,boundary);
    const double workspace_seconds=elapsed(start);double error=0.0;
    for(std::size_t c=0;c<out.size();++c)error=std::max(error,std::abs(out[c]-reference[c]));
    double integral=0.0,net_boundary=0.0;
    for(std::size_t c=0;c<out.size();++c)integral+=out[c]*mesh.cells()[c].volume;
    for(std::size_t f=0;f<flux.size();++f)if(mesh.faces()[f].boundary())net_boundary+=flux[f]*workspace.face_values[f];
    std::cout<<"case=fvm-workspace cells="<<mesh.cell_count()<<" threads="<<cfd::core::cpu_thread_capacity()
        <<" repeats=100 reference_seconds="<<reference_seconds<<" workspace_seconds="<<workspace_seconds
        <<" workspace_bytes="<<workspace.allocated_bytes()<<" max_parity_error="<<error
        <<" conservation_error="<<std::abs(integral-net_boundary)<<'\n';
    const auto small=make_cartesian_hexa_mesh(16,8,2);ScalarTransport cached(small),rebuilt(small);
    cached.initialize([](Vec3 x){return 1.0+x.x;});rebuilt.initialize([](Vec3 x){return 1.0+x.x;});
    const std::vector<double> no_flux(small.face_count(),0.0);
    start=Clock::now();for(int i=0;i<50;++i){rebuilt.set_face_flux(no_flux);rebuilt.step();}
    const double rebuild_seconds=elapsed(start);start=Clock::now();cached.run(50);const double cache_seconds=elapsed(start);
    double cache_error=0.0;for(std::size_t i=0;i<small.cell_count();++i)cache_error=std::max(cache_error,std::abs(cached.values()[i]-rebuilt.values()[i]));
    std::cout<<"case=fvm-operator-cache steps=50 rebuild_seconds="<<rebuild_seconds<<" cached_seconds="<<cache_seconds
        <<" operator_assemblies="<<cached.operator_assemblies()<<" max_parity_error="<<cache_error<<'\n';
    return error<1.0e-12&&cache_error<1.0e-12&&std::abs(integral-net_boundary)<1.0e-12?0:1;
}


int em_frequency1d(){
    using namespace cfd::em;
    FrequencyDomain1DConfig config;config.points=121U;config.length_m=1.0;config.frequency_hz=8.0e7;
    std::vector<Complex> current(config.points);
    for(std::size_t i=0;i<current.size();++i){const double x=static_cast<double>(i)/static_cast<double>(current.size()-1U);current[i]=std::sin(std::numbers::pi*x);}
    const auto start=Clock::now();const auto result=solve_pec_driven_maxwell_1d(config,current);
    constexpr double mu0=1.25663706212e-6,epsilon0=8.8541878128e-12;
    const double omega=2.0*std::numbers::pi*config.frequency_hz,k2=omega*omega*mu0*epsilon0,lambda=std::numbers::pi*std::numbers::pi;
    const Complex exact=Complex{0.0,omega*mu0}/(k2-lambda);const double relative=std::abs(result.electric_v_per_m[config.points/2U]-exact)/std::abs(exact);
    const auto modes=pec_cavity_eigenmodes_1d(81U,1.0,4.0,1.0,2U);
    std::cout<<"case=em-frequency1d seconds="<<elapsed(start)<<" manufactured_relative_error="<<relative
             <<" cavity_f1_hz="<<modes[0].frequency_hz<<" cavity_f2_hz="<<modes[1].frequency_hz<<'\n';
    return relative<4.0e-4?0:1;
}
int em_edge2d(){
    using namespace cfd::em;const auto mesh=cfd::fem::make_rectangle_tri_mesh(6U,6U,1.0,0.8);
    EdgeMaxwell2DConfig config;config.frequency_hz=8.0e7;constexpr double mu0=1.25663706212e-6,epsilon0=8.8541878128e-12;
    const double omega=2.0*std::numbers::pi*config.frequency_hz,lambda=std::numbers::pi*std::numbers::pi;
    const double coefficient=lambda/mu0-omega*omega*epsilon0;
    const auto start=Clock::now();const auto result=solve_driven_edge_maxwell_2d(mesh,config,[=](cfd::fem::Node2 point)->ComplexVec2{return {Complex{},Complex{0.0,coefficient*std::sin(std::numbers::pi*point.x)/omega}};});
    double error2=0.0,reference2=0.0;for(std::size_t e=0;e<mesh.triangles.size();++e){double x=0.0;for(auto n:mesh.triangles[e].node)x+=mesh.nodes[n].x/3.0;const Complex exact{std::sin(std::numbers::pi*x),0.0};error2+=std::norm(result.electric_centroid_v_per_m[e][0])+std::norm(result.electric_centroid_v_per_m[e][1]-exact);reference2+=std::norm(exact);}
    const double relative=std::sqrt(error2/reference2);std::cout<<"case=em-edge2d seconds="<<elapsed(start)<<" edges="<<result.edges.size()<<" relative_rms_error="<<relative<<'\n';return relative<0.15?0:1;
}

int em_edge3d(){
    using namespace cfd::em;
    const auto mesh=cfd::fem::make_box_tet_mesh(2U,2U,2U,1.0,0.8,1.0);
    EdgeMaxwell3DConfig config;config.frequency_hz=8.0e7;config.linear.relative_tolerance=1.0e-8;config.linear.gmres_restart=300U;
    constexpr double mu0=1.25663706212e-6,epsilon0=8.8541878128e-12;
    const double omega=2.0*std::numbers::pi*config.frequency_hz,lambda=2.0*std::numbers::pi*std::numbers::pi;
    const double coefficient=lambda/mu0-omega*omega*epsilon0;
    const auto start=Clock::now();const auto result=solve_driven_edge_maxwell_3d(mesh,config,[=](cfd::fem::Point3 p)->ComplexVec3{
        return {Complex{},Complex{0.0,coefficient*std::sin(std::numbers::pi*p.x)*std::sin(std::numbers::pi*p.z)/omega},Complex{}};
    });
    double q=0.0;const auto quality=resonator_quality_3d(mesh,config,result,5.8e7);q=quality.quality_factor;
    std::cout<<"case=em-edge3d seconds="<<elapsed(start)<<" edges="<<result.edges.size()
             <<" gmres_iterations="<<result.linear_result.iterations<<" q_wall="<<q<<'\n';
    return result.linear_result.converged&&std::isfinite(q)&&q>0.0?0:1;
}
int em_waveport(){
    using namespace cfd::em;const double a=22.86e-3,b=10.16e-3;const auto start=Clock::now();
    const auto modes=rectangular_waveguide_modes(a,b,10.0e9,1.0,1.0,4U);const auto& fundamental=modes.front();
    const double power=rectangular_mode_power_w(fundamental,a,b,96U,64U);
    std::cout<<"case=em-waveport seconds="<<elapsed(start)<<" m="<<fundamental.m<<" n="<<fundamental.n
             <<" cutoff_hz="<<fundamental.cutoff_frequency_hz<<" beta="<<fundamental.propagation_constant_rad_per_m.real()
             <<" normalized_power_w="<<power<<'\n';
    return fundamental.m==1U&&fundamental.n==0U&&fundamental.propagating&&std::abs(power-1.0)<1.0e-10?0:1;
}

int particle_pic1d(){
    using namespace cfd::particle;ChargedParticle charged;charged.velocity_m_per_s={1.0,0.2,0.0};charged.charge_c=1.0;charged.mass_kg=1.0;const double initial=std::hypot(charged.velocity_m_per_s.x,charged.velocity_m_per_s.y);
    const auto start=Clock::now();for(std::size_t i=0;i<5000U;++i)boris_push(charged,{{0,0,0},{0,0,1}},1.0e-3);
    ElectrostaticPic1DConfig config;config.grid_points=32U;config.dt_s=1.0e-3;ElectrostaticPic1D pic(config);std::vector<PicParticle1D> particles;for(std::size_t i=0;i<config.grid_points;++i)particles.push_back({static_cast<double>(i)/static_cast<double>(config.grid_points),0.1,1.0e-12,1.0,1.0});pic.set_particles(std::move(particles));pic.step(2U);double peak=0.0;for(double e:pic.electric_field())peak=std::max(peak,std::abs(e));const double speed_error=std::abs(std::hypot(charged.velocity_m_per_s.x,charged.velocity_m_per_s.y)-initial);
    std::cout<<"case=particle-pic1d seconds="<<elapsed(start)<<" boris_speed_error="<<speed_error<<" neutralized_peak_field="<<peak<<'\n';return speed_error<1.0e-12&&peak<1.0e-8?0:1;
}

int particle_em_pic1d(){
    using namespace cfd::particle;
    ElectromagneticPic1DConfig config;config.grid_points=32U;config.length_m=1.0;config.solve_longitudinal_poisson=false;config.dt_s=0.25*(config.length_m/static_cast<double>(config.grid_points))/299792458.0;
    ElectromagneticPic1D solver(config);solver.initialize_right_traveling_mode(1.0,1U);const double energy0=solver.diagnostics().field_energy_j;
    std::vector<ChargedParticle> beam(config.grid_points);
    for(std::size_t i=0;i<beam.size();++i){beam[i].mass_kg=1.0e-10;beam[i].charge_c=1.0e-16;beam[i].position_m={static_cast<double>(i)/static_cast<double>(config.grid_points),0.0,0.0};beam[i].velocity_m_per_s={0.0,1.0e5,0.0};}
    solver.set_particles(std::move(beam));const auto start=Clock::now();solver.step(4U);const auto diag=solver.diagnostics();
    double jy=0.0;for(double value:solver.current_y())jy+=value;jy/=static_cast<double>(solver.current_y().size());
    std::cout<<"case=particle-em-pic1d seconds="<<elapsed(start)<<" particles="<<diag.particle_count
             <<" field_energy_j="<<diag.field_energy_j<<" mean_jy="<<jy
             <<" continuity_linf="<<diag.charge_continuity_linf_residual<<'\n';
    return diag.particle_count==config.grid_points&&diag.field_energy_j>0.5*energy0&&diag.charge_continuity_linf_residual<1.0e-9?0:1;
}


int particle_pic2d(){
    using namespace cfd::particle;
    ElectrostaticPic2DConfig config;config.nx=18U;config.ny=14U;config.length_x_m=1.0;config.length_y_m=0.8;config.dt_s=1.0e-5;
    ElectrostaticPic2D solver(config);std::vector<PicParticle2D> particles;
    for(std::size_t i=0;i<config.nx;++i){
        PicParticle2D p;p.mass_kg=1.0;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={config.length_x_m*(static_cast<double>(i)+0.3)/static_cast<double>(config.nx),0.31+0.02*static_cast<double>(i%config.ny)};particles.push_back(p);
    }
    solver.set_particles(std::move(particles));const auto start=Clock::now();solver.deposit_and_solve();double peak=0.0;
    for(double value:solver.electric_x())peak=std::max(peak,std::abs(value));for(double value:solver.electric_y())peak=std::max(peak,std::abs(value));
    std::vector<PicParticle2D> old_particles=solver.particles(),new_particles=old_particles;
    for(std::size_t i=0;i<new_particles.size();++i){new_particles[i].position_m.x+=0.01;new_particles[i].position_m.y-=0.005;}
    const auto current=deposit_charge_conserving_current_2d(old_particles,new_particles,config.nx,config.ny,config.length_x_m,config.length_y_m,1.0e-9);
    solver.step(1U);std::cout<<"case=particle-pic2d seconds="<<elapsed(start)<<" particles="<<solver.particles().size()
        <<" peak_field="<<peak<<" continuity_linf="<<current.continuity_linf_residual<<" time_s="<<solver.time_s()<<'\n';
    return peak>0.0&&current.continuity_linf_residual<1.0e-8&&std::abs(solver.time_s()-config.dt_s)<1.0e-15?0:1;
}


int particle_pic3d(){
    using namespace cfd::particle;
    ElectrostaticPic3DConfig config;config.nx=8U;config.ny=7U;config.nz=6U;config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;config.dt_s=1.0e-6;
    ElectrostaticPic3D solver(config);std::vector<PicParticle3D> particles;
    for(std::size_t i=0;i<10U;++i){PicParticle3D p;p.mass_kg=1.0;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={0.07+0.083*static_cast<double>(i),0.11+0.057*static_cast<double>(i),0.05+0.041*static_cast<double>(i)};p.velocity_m_per_s={10.0,-4.0,2.0};particles.push_back(p);} 
    solver.set_particles(std::move(particles));const auto start=Clock::now();solver.deposit_and_solve();double peak=0.0;for(double value:solver.electric_x())peak=std::max(peak,std::abs(value));for(double value:solver.electric_y())peak=std::max(peak,std::abs(value));for(double value:solver.electric_z())peak=std::max(peak,std::abs(value));
    auto old_particles=solver.particles();auto new_particles=old_particles;for(auto& particle:new_particles){particle.position_m.x+=0.006;particle.position_m.y-=0.004;particle.position_m.z+=0.003;}
    const auto current=deposit_charge_conserving_current_3d(old_particles,new_particles,config.nx,config.ny,config.nz,config.length_x_m,config.length_y_m,config.length_z_m,1.0e-9);
    solver.step(1U);std::cout<<"case=particle-pic3d seconds="<<elapsed(start)<<" particles="<<solver.particles().size()
        <<" peak_field="<<peak<<" continuity_linf="<<current.continuity_linf_residual<<" time_s="<<solver.time_s()<<'\n';
    return peak>0.0&&current.continuity_linf_residual<1.0e-8&&std::abs(solver.time_s()-config.dt_s)<1.0e-15?0:1;
}



int particle_em_pic3d(){
    using namespace cfd::particle;
    ElectromagneticPic3DConfig config;config.nx=7U;config.ny=6U;config.nz=5U;config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;config.solve_longitudinal_poisson=false;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);const double dz=config.length_z_m/static_cast<double>(config.nz);
    config.dt_s=0.03/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
    ElectromagneticPic3D solver(config);solver.initialize_z_polarized_mode(0.5,1U,0U,0U);const double energy0=solver.diagnostics().field_energy_j;
    std::vector<PicParticle3D> particles;for(std::size_t i=0;i<6U;++i){PicParticle3D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={0.08+0.11*static_cast<double>(i),0.10+0.06*static_cast<double>(i),0.07+0.04*static_cast<double>(i)};p.velocity_m_per_s={1.0e4,-6.0e3,2.0e4};particles.push_back(p);} 
    solver.set_particles(std::move(particles));const auto start=Clock::now();solver.step(3U);const auto diag=solver.diagnostics();double max_current=0.0;for(double value:solver.current_x())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_y())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_z())max_current=std::max(max_current,std::abs(value));
    std::cout<<"case=particle-em-pic3d seconds="<<elapsed(start)<<" particles="<<diag.particle_count
             <<" field_energy_j="<<diag.field_energy_j<<" max_current="<<max_current
             <<" continuity_linf="<<diag.charge_continuity_linf_residual<<'\n';
    return diag.particle_count==6U&&diag.field_energy_j>0.0&&diag.field_energy_j<energy0*4.0&&max_current>0.0&&diag.charge_continuity_linf_residual<1.0e-6?0:1;
}

int particle_em_pic2d(){
    using namespace cfd::particle;
    ElectromagneticPic2DConfig config;config.nx=18U;config.ny=14U;config.length_x_m=1.0;config.length_y_m=0.8;config.solve_longitudinal_poisson=false;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);config.dt_s=0.04*std::min(dx,dy)/299792458.0;
    ElectromagneticPic2D solver(config);solver.initialize_tm_z_mode(0.75,1U,1U);const double energy0=solver.diagnostics().field_energy_j;
    std::vector<PicParticle2D> particles;
    for(std::size_t i=0;i<8U;++i){PicParticle2D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={0.08+0.1*static_cast<double>(i),0.11+0.045*static_cast<double>(i)};p.velocity_m_per_s={1.0e4,-5.0e3,2.0e4};particles.push_back(p);} 
    solver.set_particles(std::move(particles));const auto start=Clock::now();solver.step(3U);const auto diag=solver.diagnostics();
    double jz=0.0;for(double value:solver.current_z())jz=std::max(jz,std::abs(value));
    std::cout<<"case=particle-em-pic2d seconds="<<elapsed(start)<<" particles="<<diag.particle_count
             <<" field_energy_j="<<diag.field_energy_j<<" max_jz="<<jz
             <<" continuity_linf="<<diag.charge_continuity_linf_residual<<'\n';
    return diag.particle_count==8U&&diag.field_energy_j>0.5*energy0&&jz>0.0&&diag.charge_continuity_linf_residual<1.0e-7?0:1;
}


int particle_staggered_em_pic2d(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic2DConfig config;config.nx=20U;config.ny=16U;config.length_x_m=1.0;config.length_y_m=0.8;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);
    config.dt_s=0.035/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)));
    config.field_boundary.mode=GridBoundaryMode2D::absorbing_sponge;config.field_boundary.sponge_cells=3U;config.field_boundary.sponge_strength=1.2;config.periodic_particles=false;config.particle_boundary=ParticleWallMode::specular_reflect;
    StaggeredElectromagneticPic2D solver(config);solver.initialize_tm_z_mode(0.5,1U,1U);const double energy0=solver.diagnostics().field_energy_j;
    std::vector<PicParticle2D> particles;for(std::size_t i=0;i<6U;++i){PicParticle2D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={0.18+0.09*static_cast<double>(i),0.18+0.05*static_cast<double>(i)};p.velocity_m_per_s={8.0e3,-4.0e3,1.0e4};particles.push_back(p);}solver.set_particles(std::move(particles));const auto start=Clock::now();solver.step(3U);const auto diag=solver.diagnostics();double jz=0.0;for(double value:solver.current_z())jz=std::max(jz,std::abs(value));
    std::cout<<"case=particle-staggered-em-pic2d seconds="<<elapsed(start)<<" particles="<<diag.particle_count
             <<" field_energy_j="<<diag.field_energy_j<<" max_jz="<<jz
             <<" continuity_linf="<<diag.charge_continuity_linf_residual<<" absorbed="<<diag.absorbed_particles<<'\n';
    return diag.particle_count==6U&&diag.field_energy_j>0.0&&diag.field_energy_j<energy0*1.2&&jz>0.0&&diag.charge_continuity_linf_residual<1.0e-7?0:1;
}


int particle_staggered_em_pic3d(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic3DConfig config;config.nx=7U;config.ny=6U;config.nz=5U;config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);const double dz=config.length_z_m/static_cast<double>(config.nz);
    config.dt_s=0.035/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
    config.field_boundary.mode=GridBoundaryMode3D::absorbing_sponge;config.field_boundary.sponge_cells=2U;config.field_boundary.sponge_strength=1.0;config.periodic_particles=false;config.particle_boundary=ParticleWallMode::specular_reflect;
    StaggeredElectromagneticPic3D solver(config);solver.initialize_z_polarized_mode(0.45,1U,0U,0U);const double energy0=solver.diagnostics().field_energy_j;
    std::vector<PicParticle3D> particles;for(std::size_t i=0;i<5U;++i){PicParticle3D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={0.16+0.09*static_cast<double>(i),0.16+0.055*static_cast<double>(i),0.13+0.04*static_cast<double>(i)};p.velocity_m_per_s={8.0e3,-4.0e3,1.0e4};particles.push_back(p);}solver.set_particles(std::move(particles));
    const auto start=Clock::now();solver.step(3U);const auto diag=solver.diagnostics();double max_current=0.0;for(double value:solver.current_x())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_y())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_z())max_current=std::max(max_current,std::abs(value));
    std::cout<<"case=particle-staggered-em-pic3d seconds="<<elapsed(start)<<" particles="<<diag.particle_count
             <<" field_energy_j="<<diag.field_energy_j<<" max_current="<<max_current
             <<" continuity_linf="<<diag.charge_continuity_linf_residual<<" absorbed="<<diag.absorbed_particles<<'\n';
    return diag.particle_count==5U&&diag.field_energy_j>0.0&&diag.field_energy_j<energy0*2.0&&max_current>0.0&&diag.charge_continuity_linf_residual<1.0e-6?0:1;
}


int particle_local_current3d(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic3DConfig config;config.nx=7U;config.ny=6U;config.nz=5U;config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;config.current_deposition=CurrentDeposition3DMode::local_finite_volume;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);const double dz=config.length_z_m/static_cast<double>(config.nz);
    config.dt_s=0.025/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
    StaggeredElectromagneticPic3D solver(config);std::vector<PicParticle3D> particles;
    for(std::size_t i=0;i<6U;++i){PicParticle3D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.position_m={0.11+0.08*static_cast<double>(i),0.13+0.05*static_cast<double>(i),0.09+0.035*static_cast<double>(i)};p.velocity_m_per_s={7.0e3,-3.5e3,5.0e3};particles.push_back(p);}solver.set_particles(std::move(particles));
    const auto start=Clock::now();solver.step(2U);const auto diag=solver.diagnostics();double max_current=0.0;for(double value:solver.current_x())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_y())max_current=std::max(max_current,std::abs(value));for(double value:solver.current_z())max_current=std::max(max_current,std::abs(value));
    std::cout<<"case=particle-local-current3d seconds="<<elapsed(start)<<" particles="<<diag.particle_count
             <<" max_current="<<max_current<<" continuity_linf="<<diag.charge_continuity_linf_residual<<'\n';
    return diag.particle_count==6U&&max_current>0.0&&diag.charge_continuity_linf_residual<1.0e-6?0:1;
}


int particle_sort_halo3d(){
    using namespace cfd::particle;
    std::vector<PicParticle3D> particles;
    for(std::size_t i=0;i<10U;++i){
        PicParticle3D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.weight=1.0+0.02*static_cast<double>(i);
        p.position_m={0.03+0.087*static_cast<double>(i),0.04+0.069*static_cast<double>((i*3U)%8U),0.02+0.061*static_cast<double>((i*5U)%9U)};
        p.velocity_m_per_s={3.0e3,-2.0e3,1.0e3};particles.push_back(p);
    }
    particles[0].position_m.x=0.01;particles[1].position_m.x=0.99;particles[2].position_m.y=0.01;particles[3].position_m.z=0.99;
    const auto start=Clock::now();
    const auto sorted=sort_particles_by_cell_3d(particles,8U,6U,5U,1.0,0.75,0.5);
    ParticleGuardHalo3DConfig halo_config;halo_config.nx=8U;halo_config.ny=6U;halo_config.nz=5U;halo_config.length_x_m=1.0;halo_config.length_y_m=0.75;halo_config.length_z_m=0.5;halo_config.guard_cells=1U;
    const auto halo=classify_particle_guard_halos_3d(particles,halo_config);
    std::cout<<"case=particle-sort-halo3d seconds="<<elapsed(start)
             <<" particles="<<sorted.sorted_particles.size()
             <<" occupied_cells="<<sorted.occupied_cells
             <<" x_minus="<<halo.x_minus.size()<<" x_plus="<<halo.x_plus.size()
             <<" y_minus="<<halo.y_minus.size()<<" z_plus="<<halo.z_plus.size()<<'\n';
    return sorted.sorted_particles.size()==particles.size()&&sorted.occupied_cells>0U&&(!halo.x_minus.empty())&&(!halo.x_plus.empty())&&(!halo.y_minus.empty())&&(!halo.z_plus.empty())?0:1;
}


int particle_domain_exchange3d(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=5U;config.global_ny=4U;config.global_nz=3U;config.domains_x=2U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=5.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto start=Clock::now();
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& domain:domains){
        values[domain.domain_id].resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            values[domain.domain_id][(iz*domain.cells_y+iy)*domain.cells_x+ix]=100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz);
        }
    }
    const auto guarded=exchange_scalar_guard_cells_3d(values,domains,config,1U,-1.0);
    std::vector<PicParticle3D> particles(3U);particles[0].position_m={0.2,0.2,0.2};particles[1].position_m={3.7,2.5,1.0};particles[2].position_m={5.05,0.2,0.2};
    for(auto& p:particles){p.mass_kg=1.0e-6;p.charge_c=1.0e-15;p.weight=1.0;}
    const auto migration=plan_particle_domain_migration_3d(particles,config);
    double checksum=0.0;for(const auto& block:guarded){checksum+=block.values.front()+block.values.back();}
    std::size_t bucketed=0U;for(const auto& bucket:migration.particles_by_domain)bucketed+=bucket.size();
    std::cout<<"case=particle-domain-exchange3d seconds="<<elapsed(start)
             <<" domains="<<domains.size()<<" guarded_blocks="<<guarded.size()
             <<" bucketed="<<bucketed<<" checksum="<<checksum<<'\n';
    return domains.size()==4U&&guarded.size()==domains.size()&&bucketed==particles.size()&&checksum>0.0?0:1;
}


int particle_comm_exchange3d(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto start=Clock::now();
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& domain:domains){
        values[domain.domain_id].resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            values[domain.domain_id][(iz*domain.cells_y+iy)*domain.cells_x+ix]=100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz);
        }
    }
    const auto guard_messages=pack_scalar_guard_cell_messages_3d(values,domains,config,1U);
    const auto guarded=apply_scalar_guard_cell_messages_3d(values,domains,1U,-1.0,guard_messages);
    std::vector<std::vector<PicParticle3D>> particles_by_domain(domains.size());
    auto add_particle=[&](std::size_t source,double x,double y,double z){PicParticle3D p;p.position_m={x,y,z};p.velocity_m_per_s={1.0e3,0.0,0.0};p.mass_kg=1.0e-6;p.charge_c=1.0e-15;p.weight=1.0;particles_by_domain[source].push_back(p);};
    add_particle(0U,0.25,0.25,0.25);
    add_particle(0U,2.25,0.25,0.25);
    add_particle(domains.back().domain_id,6.05,0.25,0.25);
    const auto migration=pack_particle_migration_messages_3d(particles_by_domain,domains,config);
    const auto migrated=apply_particle_migration_messages_3d(migration.retained_by_domain,migration.messages);
    std::size_t particle_count=0U;for(const auto& bucket:migrated)particle_count+=bucket.size();
    std::size_t guard_entries=0U;for(const auto& message:guard_messages)guard_entries+=message.values.size();
    double checksum=0.0;for(const auto& block:guarded)checksum+=block.values.front()+block.values.back();
    std::cout<<"case=particle-comm-exchange3d seconds="<<elapsed(start)
             <<" domains="<<domains.size()<<" guard_messages="<<guard_messages.size()
             <<" guard_entries="<<guard_entries<<" particle_messages="<<migration.messages.size()
             <<" particles="<<particle_count<<" checksum="<<checksum<<'\n';
    return domains.size()==6U&&(!guard_messages.empty())&&guard_entries>0U&&(!migration.messages.empty())&&particle_count==3U&&checksum>0.0?0:1;
}


int particle_transport3d(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto start=Clock::now();
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& domain:domains){
        values[domain.domain_id].resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            values[domain.domain_id][(iz*domain.cells_y+iy)*domain.cells_x+ix]=100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz);
        }
    }
    std::vector<std::vector<PicParticle3D>> particles_by_domain(domains.size());
    auto add_particle=[&](std::size_t source,double x,double y,double z){PicParticle3D p;p.position_m={x,y,z};p.velocity_m_per_s={1.0e3,2.0e3,3.0e3};p.mass_kg=1.0e-6;p.charge_c=1.0e-15;p.weight=1.0;particles_by_domain[source].push_back(p);};
    add_particle(0U,0.25,0.25,0.25);add_particle(0U,2.25,0.25,0.25);add_particle(1U,4.25,2.25,1.25);add_particle(domains.back().domain_id,6.05,0.25,0.25);
    const auto migration=pack_particle_migration_messages_3d(particles_by_domain,domains,config);
    const auto guard_messages=pack_scalar_guard_cell_messages_3d(values,domains,config,1U);
    const auto topology=make_pic_rank_topology_3d(domains,3U,1U,1U);
    const auto exchanged=exchange_pic_messages_in_memory_3d(migration.retained_by_domain,migration.messages,values,guard_messages,domains,1U,-1.0,topology);
    std::size_t particles=0U;for(const auto& bucket:exchanged.particles_by_domain)particles+=bucket.size();
    double checksum=0.0;for(const auto& block:exchanged.guarded_blocks)checksum+=block.values.front()+block.values.back();
    std::cout<<"case=particle-transport3d seconds="<<elapsed(start)
             <<" ranks="<<topology.rank_count<<" envelopes="<<exchanged.diagnostics.envelopes
             <<" remote="<<exchanged.diagnostics.remote_rank_messages
             <<" particles="<<particles<<" guard_values="<<exchanged.diagnostics.scalar_guard_value_count
             <<" checksum="<<checksum<<'\n';
    return topology.rank_count==3U&&exchanged.diagnostics.envelopes>0U&&exchanged.diagnostics.remote_rank_messages>0U&&particles==4U&&checksum>0.0?0:1;
}

int particle_serialized_transport3d(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto start=Clock::now();
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& domain:domains){
        values[domain.domain_id].resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            values[domain.domain_id][(iz*domain.cells_y+iy)*domain.cells_x+ix]=100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz);
        }
    }
    std::vector<std::vector<PicParticle3D>> particles_by_domain(domains.size());
    auto add_particle=[&](std::size_t source,double x,double y,double z){PicParticle3D p;p.position_m={x,y,z};p.velocity_m_per_s={2.0e3,1.0e3,0.5e3};p.mass_kg=1.0e-6;p.charge_c=1.0e-15;p.weight=1.0;particles_by_domain[source].push_back(p);};
    add_particle(0U,0.25,0.25,0.25);add_particle(0U,2.25,0.25,0.25);add_particle(1U,4.25,2.25,1.25);add_particle(domains.back().domain_id,6.05,0.25,0.25);
    const auto migration=pack_particle_migration_messages_3d(particles_by_domain,domains,config);
    const auto guard_messages=pack_scalar_guard_cell_messages_3d(values,domains,config,1U);
    const auto topology=make_pic_rank_topology_3d(domains,3U,1U,1U);
    const auto envelopes=make_pic_transport_envelopes_3d(migration.messages,guard_messages,topology);
    const auto serialized=serialize_pic_transport_envelopes_3d(envelopes);
    const auto plan=plan_serialized_pic_exchange_3d(envelopes,topology.rank_count);
    const auto decoded=deserialize_pic_transport_envelopes_3d(serialized);
    std::vector<ParticleMigrationMessage3D> particle_messages;std::vector<ScalarGuardCellMessage3D> scalar_messages;
    for(const auto& envelope:decoded){if(envelope.kind==PicTransportPayloadKind3D::particle_migration)particle_messages.push_back(envelope.particle_migration);else scalar_messages.push_back(envelope.scalar_guard);}
    const auto migrated=apply_particle_migration_messages_3d(migration.retained_by_domain,particle_messages);
    const auto guarded=apply_scalar_guard_cell_messages_3d(values,domains,1U,-1.0,scalar_messages);
    std::size_t particles=0U,planned_messages=0U,planned_bytes=0U;for(const auto& bucket:migrated)particles+=bucket.size();for(std::size_t rank=0;rank<plan.rank_count;++rank){planned_messages+=plan.message_count_by_destination_rank[rank];planned_bytes+=plan.byte_count_by_destination_rank[rank];}
    double checksum=0.0;for(const auto& block:guarded)checksum+=block.values.front()+block.values.back();
    const auto diag=summarize_pic_serialized_transport_3d(serialized);
    std::cout<<"case=particle-serialized-transport3d seconds="<<elapsed(start)
             <<" ranks="<<topology.rank_count<<" serialized="<<serialized.size()
             <<" planned_messages="<<planned_messages<<" planned_bytes="<<planned_bytes
             <<" particles="<<particles<<" guard_values="<<diag.scalar_guard_value_count
             <<" checksum="<<checksum<<'\n';
    return topology.rank_count==3U&&serialized.size()==envelopes.size()&&planned_messages==envelopes.size()&&planned_bytes>0U&&particles==4U&&checksum>0.0?0:1;
}


int particle_distributed_round3d(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto start=Clock::now();
    const auto domains=make_pic_domain_grid_3d(config);
    std::vector<std::vector<double>> values(domains.size());
    for(const auto& domain:domains){
        values[domain.domain_id].resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            values[domain.domain_id][(iz*domain.cells_y+iy)*domain.cells_x+ix]=100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz);
        }
    }
    std::vector<std::vector<PicParticle3D>> particles_by_domain(domains.size());
    auto add_particle=[&](std::size_t source,double x,double y,double z){PicParticle3D p;p.position_m={x,y,z};p.velocity_m_per_s={2.0e3,1.0e3,0.5e3};p.mass_kg=1.0e-6;p.charge_c=1.0e-15;p.weight=1.0;particles_by_domain[source].push_back(p);};
    add_particle(0U,0.25,0.25,0.25);add_particle(0U,2.25,0.25,0.25);add_particle(1U,4.25,2.25,1.25);add_particle(domains.back().domain_id,6.05,0.25,0.25);
    const auto topology=make_pic_rank_topology_3d(domains,3U,1U,1U);
    const auto round=run_serialized_distributed_pic_exchange_round_3d(particles_by_domain,values,domains,config,1U,-1.0,topology);
    std::size_t particles=0U;for(const auto& bucket:round.particles_by_domain)particles+=bucket.size();
    double checksum=0.0;for(const auto& block:round.scalar_guarded_blocks)checksum+=block.values.front()+block.values.back();
    std::cout<<"case=particle-distributed-round3d seconds="<<elapsed(start)
             <<" ranks="<<topology.rank_count<<" serialized_messages="<<round.serialized_message_count
             <<" serialized_bytes="<<round.serialized_byte_count
             <<" remote="<<round.transport_diagnostics.remote_rank_messages
             <<" particles="<<particles<<" checksum="<<checksum<<'\n';
    return round.serialized_message_count==round.transport_diagnostics.envelopes&&round.serialized_byte_count>0U&&round.transport_diagnostics.remote_rank_messages>0U&&particles==4U&&checksum>0.0?0:1;
}

int particle_field_guards3d(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=6U;config.global_ny=4U;config.global_nz=3U;config.domains_x=3U;config.domains_y=2U;config.domains_z=1U;config.length_x_m=6.0;config.length_y_m=4.0;config.length_z_m=3.0;
    const auto start=Clock::now();
    const auto domains=make_pic_domain_grid_3d(config);
    auto make_component=[&](double offset){std::vector<std::vector<double>> values(domains.size());for(const auto& domain:domains){auto& block=values[domain.domain_id];block.resize(domain.cells_x*domain.cells_y*domain.cells_z);for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){block[(iz*domain.cells_y+iy)*domain.cells_x+ix]=offset+100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz);}}return values;};
    ElectromagneticFieldBlocks3D fields;fields.electric_x_by_domain=make_component(1000.0);fields.electric_y_by_domain=make_component(2000.0);fields.electric_z_by_domain=make_component(3000.0);fields.magnetic_x_by_domain=make_component(4000.0);fields.magnetic_y_by_domain=make_component(5000.0);fields.magnetic_z_by_domain=make_component(6000.0);
    const auto guarded=exchange_electromagnetic_field_guard_cells_3d(fields,domains,config,1U,-99.0);
    std::size_t guarded_values=0U;double checksum=0.0;
    for(const auto& block:guarded.electric_x){guarded_values+=block.values.size();checksum+=block.values.front()+block.values.back();}
    for(const auto& block:guarded.magnetic_z){guarded_values+=block.values.size();checksum+=block.values.front()+block.values.back();}
    std::cout<<"case=particle-field-guards3d seconds="<<elapsed(start)
             <<" domains="<<domains.size()<<" guarded_values="<<guarded_values
             <<" checksum="<<checksum<<'\n';
    return guarded.electric_x.size()==domains.size()&&guarded.magnetic_z.size()==domains.size()&&guarded_values>0U&&checksum>0.0?0:1;
}


int particle_distributed_step3d(){
    using namespace cfd::particle;
    PicDomainGrid3DConfig config;config.global_nx=4U;config.global_ny=2U;config.global_nz=2U;config.domains_x=2U;config.domains_y=1U;config.domains_z=1U;config.length_x_m=4.0;config.length_y_m=2.0;config.length_z_m=2.0;config.periodic=true;
    const auto start=Clock::now();
    const auto domains=make_pic_domain_grid_3d(config);
    auto make_component=[&](double base){std::vector<std::vector<double>> values(domains.size());for(const auto& domain:domains){auto& block=values[domain.domain_id];block.assign(domain.cells_x*domain.cells_y*domain.cells_z,base);for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){block[(iz*domain.cells_y+iy)*domain.cells_x+ix]+=0.01*(100.0*static_cast<double>(domain.first_x+ix)+10.0*static_cast<double>(domain.first_y+iy)+static_cast<double>(domain.first_z+iz));}}return values;};
    ElectromagneticFieldBlocks3D fields;fields.electric_x_by_domain=make_component(0.0);fields.electric_y_by_domain=make_component(0.0);fields.electric_z_by_domain=make_component(0.0);fields.magnetic_x_by_domain=make_component(0.0);fields.magnetic_y_by_domain=make_component(0.0);fields.magnetic_z_by_domain=make_component(0.0);
    std::vector<std::vector<PicParticle3D>> particles(domains.size());
    PicParticle3D p;p.position_m={1.75,0.6,0.6};p.velocity_m_per_s={0.5,0.0,0.0};p.mass_kg=1.0;p.charge_c=0.0;p.weight=1.0;particles[0].push_back(p);
    p.position_m={3.85,0.5,0.5};p.velocity_m_per_s={0.4,0.0,0.0};particles[1].push_back(p);
    DistributedStaggeredPicStep3DConfig step;step.dt_s=1.0;step.guard_cells=1U;step.exterior_value=-9.0;
    const auto topology=make_pic_rank_topology_3d(domains,2U,1U,1U);
    const auto result=run_serialized_distributed_staggered_pic_step_3d(particles,fields,domains,config,step,topology);
    std::cout<<"case=particle-distributed-step3d seconds="<<elapsed(start)
             <<" particles_before="<<result.particle_count_before
             <<" particles_after="<<result.particle_count_after_migration
             <<" max_displacement_m="<<result.max_particle_displacement_m
             <<" transport_messages="<<result.exchange_round.transport_diagnostics.envelopes
             <<" serialized_bytes="<<result.exchange_round.serialized_byte_count
             <<" guarded_domains="<<result.guarded_fields.electric_x.size()<<'\n';
    return result.particle_count_before==2U&&result.particle_count_after_migration==2U&&result.max_particle_displacement_m>0.0&&result.exchange_round.transport_diagnostics.envelopes>0U&&result.guarded_fields.electric_x.size()==domains.size()?0:1;
}


int particle_geant4_transport(){
    using namespace cfd::particle;
    TransportWorld world;
    world.particles.push_back({"electron",-1.602176634e-19,9.1093837e-31});
    world.materials.push_back({"argon-gas",1.6,5.0});
    world.materials.push_back({"detector",2500.0,10.0});
    world.regions.push_back({"drift-region",{{0.0,0.0,0.0},{1.0,1.0,1.0}},0U});
    world.regions.push_back({"scoring-region",{{1.0,0.0,0.0},{1.5,1.0,1.0}},1U});
    world.processes.push_back({"continuous-ionization-loss",TransportProcessKind::continuous_energy_loss,0U,0U,150.0,0.0,0.0,0.0,0.0});
    world.processes.push_back({"hard-scatter-secondary",TransportProcessKind::discrete_interaction,0U,0U,0.0,0.32,0.15,0.2,0.0});
    TransportTrack track;track.position_m={0.08,0.5,0.5};track.direction={1.0,0.0,0.0};track.kinetic_energy_ev=1200.0;track.particle_index=0U;
    TransportConfig config;config.max_step_m=0.8;config.energy_cut_ev=2.0;
    const auto start=Clock::now();const auto result=transport_track(world,track,config,8U);
    std::cout<<"case=particle-geant4-transport seconds="<<elapsed(start)
             <<" steps="<<result.scoring.steps
             <<" energy_deposit_ev="<<result.scoring.total_energy_deposit_ev
             <<" track_length_m="<<result.scoring.total_track_length_m
             <<" secondaries="<<result.scoring.secondaries
             <<" final_energy_ev="<<result.primary.kinetic_energy_ev<<'\n';
    return result.scoring.steps>0U&&result.scoring.total_energy_deposit_ev>0.0&&result.scoring.total_track_length_m>0.0&&result.scoring.secondaries>0U?0:1;
}


int particle_transport_dose_bvh(){
    using namespace cfd::particle;using namespace cfd::multiphysics;
    TransportWorld world;world.particles.push_back({"electron",-1.602176634e-19,9.1093837e-31});world.materials.push_back({"gas",1.2,2.0});world.sensitive_detectors.push_back({"dose",0.0,false});
    for(std::size_t i=0;i<4U;++i)world.regions.push_back({"slab",{{0.25*static_cast<double>(i),0.0,0.0},{0.25*static_cast<double>(i+1U),1.0,1.0}},0U,0U});
    const auto bvh=build_region_bvh(world,2U);
    TransportPhysicsList physics;physics.name="em-standard-lite";physics.production_cut_energy_ev=1.0;physics.processes.push_back({"dEdx",TransportProcessKind::continuous_energy_loss,0U,0U,300.0,0.0,0.0,0.0,0.0});physics.processes.push_back({"ionization",TransportProcessKind::discrete_interaction,0U,0U,0.0,0.03,0.12,0.1,0.0});world=apply_physics_list(std::move(world),physics);
    TransportTrack track;track.position_m={0.05,0.5,0.5};track.direction={1.0,0.0,0.0};track.kinetic_energy_ev=2000.0;TransportConfig config;config.max_step_m=0.18;TransportRandom rng; rng.state=77U;
    const auto start=Clock::now();const auto result=transport_track_stochastic(world,track,config,rng,16U);const auto hits=collect_transport_hits(world,result.steps);auto grid=make_transport_dose_grid(4U,4U,4U,{{0.0,0.0,0.0},{1.0,1.0,1.0}},1000.0);score_hits_to_dose_grid(hits,grid);const auto sar=project_dose_grid_to_pennes2d_sar(grid,1.0,6U,6U);PennesBioheat2DConfig bio;bio.nx=6U;bio.ny=6U;bio.dt_s=0.01;bio.blood_temperature_k=310.0;PennesBioheat2D solver(bio);solver.initialize(310.0);solver.set_sar(sar);solver.step();
    const std::size_t lookup=locate_region_bvh(world,bvh,{0.62,0.5,0.5});double max_sar=0.0;for(double value:sar)max_sar=std::max(max_sar,value);
    std::cout<<"case=particle-transport-dose-bvh seconds="<<elapsed(start)<<" steps="<<result.scoring.steps<<" hits="<<hits.hits.size()<<" bvh_lookup="<<lookup<<" max_sar_w_per_kg="<<max_sar<<" bioheat_iterations="<<solver.linear_result().iterations<<'\n';
    return result.scoring.steps>0U&&!hits.hits.empty()&&lookup==2U&&max_sar>0.0&&solver.linear_result().converged?0:1;
}

int particle_campaign_doe(){
    using namespace cfd::workflow;
    const auto start=Clock::now();std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4"}},{"aoa",CampaignParameterKind::discrete,0.0,0.0,{"0","2"}}};const auto factorial=generate_factorial_campaign(factors);
    std::vector<CampaignParameter> lhs{{"re",CampaignParameterKind::continuous,1.0e5,2.0e5,{}},{"turbulence",CampaignParameterKind::discrete,0.0,0.0,{"sa","sst"}}};const auto latin=generate_latin_hypercube_campaign(lhs,6U,{123U});
    const auto rendered=render_campaign_template("MACH={mach}\n<!-- IF aoa=2 -->ANGLE=high\n<!-- ENDIF -->",factorial.back());
    SolverAdapterDescriptor adapter{"cfd_solvers","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::restart},{"stop","checkpoint"}};validate_solver_adapter_descriptor(adapter);
    auto objective=[](std::span<const double> x){return (x[0]-1.0)*(x[0]-1.0)+0.25*(x[1]+2.0)*(x[1]+2.0);};std::vector<double> x{0.0,-1.0};const auto gradient=finite_difference_gradient(objective,x);const auto next=gradient_descent_update(x,gradient,0.2);
    std::vector<CampaignRegistryEntry> registry;update_case_status(registry,{factorial[0].case_id,CampaignCaseStatus::done,1.0,10U,"ok"});update_case_status(registry,{factorial[1].case_id,CampaignCaseStatus::running,0.0,3U,"running"});const auto summary=summarize_campaign_registry(registry);
    std::cout<<"case=particle-campaign-doe seconds="<<elapsed(start)<<" factorial_cases="<<factorial.size()<<" lhs_cases="<<latin.size()<<" rendered_bytes="<<rendered.size()<<" gradient0="<<gradient[0]<<" best="<<summary.best_case_id<<'\n';
    return factorial.size()==4U&&latin.size()==6U&&solver_has_capability(adapter,SolverCapability::restart)&&objective(next)<objective(x)&&summary.done==1U&&summary.running==1U?0:1;
}


int particle_campaign_execution(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0110_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.3","0.6"}},{"aoa",CampaignParameterKind::discrete,0.0,0.0,{"0","4"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n<!-- IF aoa=4 -->\nangle={aoa}\n<!-- ENDIF -->"},{"run.cfg","case={mach}\n"}};
    CampaignWriteOptions options;options.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,options);
    auto registry=load_campaign_registry(root/"registry.tsv");
    SolverAdapterDescriptor adapter{"external-stub","stub-solver","stub-solver",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","checkpoint"}};
    SolverRunRequest request{adapter,root/"case0001",SolverRuntime::native,2U,1U,{"--dry-run"},{},{},{}};
    const auto plan=build_solver_command_plan(request);
    const std::string log="iter=1 residual=1e-2 equation=flow\niter=2 residual=1e-4 equation=flow\nperf.wall_s=0.25 s\nDONE\n";
    const auto residuals=parse_residual_history(log);const auto metrics=parse_performance_metrics(log);
    registry[0].status=detect_solver_outcome_from_log(log);registry[0].objective=residuals.back().residual;registry[0].iterations=residuals.back().iteration;registry[0].message="parsed external solver log";
    save_campaign_registry(root/"registry.tsv",registry);
    auto objective=[](std::span<const double> x){return (x[0]-0.75)*(x[0]-0.75)+(x[1]-0.1)*(x[1]-0.1);};
    CampaignOptimizationConfig opt;opt.max_iterations=40U;opt.step_size=0.25;const auto opt_result=run_gradient_descent_campaign(objective,{0.0,1.0},opt);
    const auto summary=summarize_campaign_registry(load_campaign_registry(root/"registry.tsv"));
    std::cout<<"case=particle-campaign-execution seconds="<<elapsed(start)
             <<" cases="<<persisted.cases_written
             <<" files="<<persisted.files_written
             <<" command='"<<plan.display_command<<"'"
             <<" residuals="<<residuals.size()
             <<" metrics="<<metrics.size()
             <<" best_case="<<summary.best_case_id
             <<" optimized_objective="<<opt_result.objective<<'\n';
    std::filesystem::remove_all(root);
    return persisted.cases_written==4U&&persisted.files_written==12U&&registry[0].status==CampaignCaseStatus::done&&residuals.size()==2U&&metrics.size()==1U&&opt_result.objective<1.0e-4?0:1;
}


int particle_campaign_local_runner(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0111_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.25","0.5"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    const auto script=root/"stub_solver.sh";
    {std::ofstream out(script);out<<"#!/bin/sh\n"
        <<"echo 'iter=1 residual=2e-2 equation=flow'\n"
        <<"echo 'iter=2 residual=2e-4 equation=flow'\n"
        <<"echo 'perf.wall_time 0.01 s'\n"
        <<"echo 'DONE'\n"
        <<"echo 'ITER 3 RESIDUAL 2e-5 FIELD flow' > residual_history.log\n"
        <<"echo 'perf.cells_per_s 500 cells/s' > timing.log\n"
        <<"exit 0\n";}
    std::filesystem::permissions(script,std::filesystem::perms::owner_exec|std::filesystem::perms::group_exec|std::filesystem::perms::others_exec,std::filesystem::perm_options::add);
    SolverAdapterDescriptor adapter{"stub",script.string(),"stub",{SolverCapability::residuals,SolverCapability::performance},{"stop"}};
    SolverRunRequest request{adapter,root/"case0001",SolverRuntime::native,1U,1U,{},{},{},{}};
    const auto doctor=doctor_solver_runtime(request);
    LocalCampaignRunOptions run_options;run_options.adapter=adapter;run_options.runtime=SolverRuntime::native;run_options.mpi_ranks=1U;run_options.threads=1U;
    const auto run=run_local_campaign(root,run_options);
    const auto registry=load_campaign_registry(root/"registry.tsv");
    std::size_t residuals=0U,performance=0U;for(const auto& item:run.case_results){residuals+=item.residuals.size();performance+=item.performance.size();}
    std::cout<<"case=particle-campaign-local-runner seconds="<<elapsed(start)
             <<" cases="<<persisted.cases_written
             <<" doctor_ok="<<doctor.ok
             <<" launched="<<run.launched
             <<" done="<<run.summary.done
             <<" residuals="<<residuals
             <<" performance="<<performance
             <<" first_status="<<campaign_status_name(registry.front().status)<<'\n';
    std::filesystem::remove_all(root);
    return doctor.ok&&persisted.cases_written==2U&&run.launched==2U&&run.summary.done==2U&&residuals>=6U&&performance>=4U?0:1;
}

int particle_campaign_control(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0112_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"case_param",CampaignParameterKind::discrete,0.0,0.0,{"A","B"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","case={case_param}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    SolverAdapterDescriptor adapter{"stub","/bin/echo","stub",{SolverCapability::control,SolverCapability::residuals,SolverCapability::performance},{"stop","extend","checkpoint","flush"}};
    const auto directive=write_campaign_control_directive(root,adapter,{"case0001",CampaignControlAction::checkpoint,"after-iteration=10","smoke test"});
    const auto directives=discover_campaign_control_directives(root/"case0001");
    auto registry=load_campaign_registry(root/"registry.tsv");
    const auto jobs=parse_slurm_queue_table("JOBID CASE STATE ELAPSED\n17 case0001 RUNNING 4\n18 case0002 COMPLETED 9\n");
    const auto scheduler_summary=apply_scheduler_records_to_registry(registry,jobs);
    {std::ofstream out(root/"case0001"/"solver.stdout.log");out<<"ITER 4 RESIDUAL 3e-4 FIELD flow\nperf.wall_time 0.2 s\nDONE\n";}
    const auto refreshed=refresh_campaign_status_from_outputs(root,registry);
    save_campaign_registry(root/"registry.tsv",registry);
    const bool ok=std::filesystem::exists(directive.path)&&persisted.cases_written==2U&&directives.size()>=2U&&jobs.size()==2U&&refreshed.done==2U&&registry.front().iterations==4U;
    std::cout<<"case=particle-campaign-control seconds="<<elapsed(start)
             <<" directive_exists="<<std::filesystem::exists(directive.path)
             <<" directives="<<directives.size()
             <<" scheduler_jobs="<<jobs.size()
             <<" scheduler_running="<<scheduler_summary.running
             <<" refreshed_done="<<refreshed.done
             <<" first_status="<<campaign_status_name(registry.front().status)
             <<" first_iterations="<<registry.front().iterations<<'\n';
    std::filesystem::remove_all(root);
    return ok?0:1;
}


int particle_multiserver_deploy(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0113_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4","0.6"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    const auto registry=load_campaign_registry(root/"registry.tsv");
    SolverAdapterDescriptor adapter{"stub","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a","10.0.0.11","cfd","/srv/cfd/campaigns",CampaignServerRuntime::docker,2U,{"gpu","linux"},"cfd-solvers:local",true},
        {"cpu-a","localhost","","/tmp/cfd_remote",CampaignServerRuntime::native,2U,{"cpu","linux"},"",true}
    };
    const auto plan=plan_multi_server_campaign(root,registry,adapter,servers);
    const auto written=write_multiserver_plan_files(root,plan);
    DockerDeployConfig deploy;deploy.image_name="cfd-solvers:local";deploy.worker_replicas=2U;deploy.worker_cpus=2U;deploy.worker_memory="2g";
    const auto deploy_result=write_docker_deploy_system(root/"docker",deploy,true);
    const bool compose_exists=std::filesystem::exists(root/"docker"/"compose.yaml");
    std::cout<<"case=particle-multiserver-deploy seconds="<<elapsed(start)
             <<" cases="<<persisted.cases_written
             <<" assignments="<<plan.assignments.size()
             <<" commands="<<plan.commands.size()
             <<" active_servers="<<plan.active_servers
             <<" deploy_files="<<deploy_result.files_written
             <<" compose_exists="<<compose_exists
             <<" commands_file="<<written.commands_path.filename().string()<<'\n';
    std::filesystem::remove_all(root);
    return persisted.cases_written==3U&&plan.assignments.size()==3U&&plan.commands.size()==9U&&deploy_result.files_written>=6U&&compose_exists?0:1;
}


int particle_multiserver_execution(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0114_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4","0.6","0.8"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    auto registry=load_campaign_registry(root/"registry.tsv");
    SolverAdapterDescriptor adapter{"stub","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a","10.0.0.11","cfd","/srv/cfd/campaigns",CampaignServerRuntime::docker,2U,{"gpu","linux"},"cfd-solvers:local",true},
        {"cpu-a","localhost","","/tmp/cfd_remote",CampaignServerRuntime::native,2U,{"cpu","linux"},"",true}
    };
    const auto campaign_plan=plan_multi_server_campaign(root,registry,adapter,servers);
    const auto execution_plan=plan_multiserver_execution(campaign_plan,servers);
    const auto written=write_multiserver_execution_files(root,execution_plan);
    const auto statuses=parse_remote_job_status_table("case_id\tserver\tstate\texit_code\ncase0001\tcpu-a\trunning\t0\ncase0002\tgpu-a\tdone\t0\ncase0003\tgpu-a\tfailed\t2\ncase0004\tcpu-a\tcancelled\t0\n");
    const auto summary=apply_remote_job_status_to_registry(registry,statuses);
    save_campaign_registry(root/"registry.tsv",registry);
    const bool health_exists=std::filesystem::exists(written.health_script_path);
    const bool launch_exists=std::filesystem::exists(written.launch_script_path);
    std::cout<<"case=particle-multiserver-execution seconds="<<elapsed(start)
             <<" cases="<<persisted.cases_written
             <<" jobs="<<execution_plan.jobs.size()
             <<" health_checks="<<execution_plan.health_checks.size()
             <<" docker_jobs="<<execution_plan.docker_jobs
             <<" native_jobs="<<execution_plan.native_jobs
             <<" remote_jobs="<<execution_plan.remote_jobs
             <<" summary_done="<<summary.done
             <<" summary_failed="<<summary.failed
             <<" health_script="<<health_exists
             <<" launch_script="<<launch_exists<<'\n';
    std::filesystem::remove_all(root);
    return persisted.cases_written==4U&&execution_plan.jobs.size()==4U&&execution_plan.health_checks.size()==6U&&summary.done==1U&&summary.failed==1U&&health_exists&&launch_exists?0:1;
}


int particle_multiserver_supervision(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0115_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4","0.6","0.8"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    const auto registry=load_campaign_registry(root/"registry.tsv");
    SolverAdapterDescriptor adapter{"stub","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a","10.0.0.11","cfd","/srv/cfd/campaigns",CampaignServerRuntime::docker,2U,{"gpu","linux"},"cfd-solvers:local",true},
        {"cpu-a","localhost","","/tmp/cfd_remote",CampaignServerRuntime::native,2U,{"cpu","linux"},"",true}
    };
    const auto campaign_plan=plan_multi_server_campaign(root,registry,adapter,servers);
    const auto execution_plan=plan_multiserver_execution(campaign_plan,servers);
    RemoteSupervisionOptions options;options.max_attempts=3U;options.retry_delay_seconds=1U;options.log_tail_lines=32U;
    const auto supervision=plan_multiserver_supervision(execution_plan,servers,options);
    const auto written=write_multiserver_supervision_files(root,supervision);
    const auto secret=redact_sensitive_command_display("docker run -e TOKEN=abcd -e PASSWORD=hunter2 cfd",options.sensitive_markers);
    const bool access_exists=std::filesystem::exists(written.access_script_path);
    const bool retry_exists=std::filesystem::exists(written.retry_launch_script_path);
    const bool dashboard_exists=std::filesystem::exists(written.dashboard_json_path);
    std::cout<<"case=particle-multiserver-supervision seconds="<<elapsed(start)
             <<" cases="<<persisted.cases_written
             <<" probes="<<supervision.access_probes.size()
             <<" retry_launches="<<supervision.retry_launches.size()
             <<" tail_pairs="<<supervision.log_tails.size()
             <<" dashboard_cases="<<supervision.dashboard_cases.size()
             <<" access_script="<<access_exists
             <<" retry_script="<<retry_exists
             <<" dashboard="<<dashboard_exists
             <<" redacted="<<(secret.find("hunter2")==std::string::npos)<<'\n';
    std::filesystem::remove_all(root);
    return persisted.cases_written==4U&&supervision.retry_launches.size()==4U&&supervision.log_tails.size()==4U&&supervision.dashboard_cases.size()==4U&&access_exists&&retry_exists&&dashboard_exists&&secret.find("hunter2")==std::string::npos?0:1;
}


int particle_multiserver_controller(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0116_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4","0.6","0.8"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    const auto registry=load_campaign_registry(root/"registry.tsv");
    SolverAdapterDescriptor adapter{"stub","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","checkpoint","flush"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a","10.0.0.11","cfd","/srv/cfd/campaigns",CampaignServerRuntime::docker,2U,{"gpu","linux"},"cfd-solvers:local",true},
        {"cpu-a","localhost","","/tmp/cfd_remote",CampaignServerRuntime::native,2U,{"cpu","linux"},"",true}
    };
    const auto campaign_plan=plan_multi_server_campaign(root,registry,adapter,servers);
    const auto execution_plan=plan_multiserver_execution(campaign_plan,servers);
    RemoteSupervisionOptions sup_options;sup_options.log_tail_lines=32U;
    const auto supervision=plan_multiserver_supervision(execution_plan,servers,sup_options);
    RemoteControllerConfig controller;controller.port=9090U;controller.log_tail_lines=32U;controller.control_actions={"stop","checkpoint","flush"};
    const auto plan=plan_multiserver_controller(supervision,controller);
    const auto written=write_multiserver_controller_files(root,plan);
    const bool api_exists=std::filesystem::exists(written.openapi_path);
    const bool script_exists=std::filesystem::exists(written.controller_script_path);
    const bool status_exists=std::filesystem::exists(written.status_path);
    std::cout<<"case=particle-multiserver-controller seconds="<<elapsed(start)
             <<" cases="<<persisted.cases_written
             <<" routes="<<plan.routes.size()
             <<" controller_cases="<<plan.cases.size()
             <<" api="<<api_exists
             <<" script="<<script_exists
             <<" status="<<status_exists
             <<" token_route="<<(plan.openapi_json.find("requires-token")!=std::string::npos)<<'\n';
    std::filesystem::remove_all(root);
    return persisted.cases_written==4U&&plan.routes.size()==12U&&plan.cases.size()==4U&&api_exists&&script_exists&&status_exists?0:1;
}

int particle_multiserver_dashboard(){
    using namespace cfd::workflow;
    const auto start=Clock::now();
    const auto root=std::filesystem::temp_directory_path()/"cfd_solvers_v0117_cli_campaign";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"angle",CampaignParameterKind::discrete,0.0,0.0,{"0","5","10"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","angle={angle}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    const auto registry=load_campaign_registry(root/"registry.tsv");
    SolverAdapterDescriptor adapter{"stub","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","extend","checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-dashboard","10.0.0.42","cfd","/srv/cfd/campaigns",CampaignServerRuntime::docker,2U,{"gpu","linux"},"cfd-solvers:dash",true},
        {"local-dashboard","localhost","","/tmp/cfd_dashboard",CampaignServerRuntime::native,1U,{"cpu","linux"},"",true}
    };
    const auto campaign_plan=plan_multi_server_campaign(root,registry,adapter,servers);
    const auto execution_plan=plan_multiserver_execution(campaign_plan,servers);
    const auto supervision=plan_multiserver_supervision(execution_plan,servers);
    RemoteControllerConfig controller;controller.port=9191U;controller.dashboard_title="CFD campaign dashboard";controller.dashboard_refresh_seconds=7U;controller.control_actions={"stop","extend","checkpoint"};
    const auto plan=plan_multiserver_controller(supervision,controller);
    const auto written=write_multiserver_controller_files(root,plan);
    const bool html_exists=std::filesystem::exists(written.index_path);
    const bool js_exists=std::filesystem::exists(written.dashboard_js_path);
    const bool css_exists=std::filesystem::exists(written.dashboard_css_path);
    const bool events_exists=std::filesystem::exists(written.events_path);
    std::cout<<"case=particle-multiserver-dashboard seconds="<<elapsed(start)
             <<" cases="<<persisted.cases_written
             <<" routes="<<plan.routes.size()
             <<" assets="<<written.static_assets_written
             <<" html="<<html_exists
             <<" js="<<js_exists
             <<" css="<<css_exists
             <<" events="<<events_exists<<'\n';
    std::filesystem::remove_all(root);
    return persisted.cases_written==3U&&plan.routes.size()==12U&&written.static_assets_written==7U&&html_exists&&js_exists&&css_exists&&events_exists?0:1;
}

int particle_plasma_chemistry(){
    using namespace cfd::particle;
    constexpr double qe=1.602176634e-19;
    PlasmaState state;state.species={{"e",-qe,9.1093837015e-31},{"Ar",0.0,6.6335209e-26},{"Ar+",qe,6.6335209e-26}};state.number_density_m3={1.0e15,1.0e20,1.0e15};state.electron_temperature_ev=8.0;
    PlasmaReaction ionization;ionization.reactants={{0U,1.0},{1U,1.0}};ionization.products={{0U,2.0},{2U,1.0}};ionization.rate_coefficient=1.0e-20;
    const auto start=Clock::now();const auto diagnostics=advance_plasma_chemistry(state,std::span<const PlasmaReaction>(&ionization,1U),1.0e-7,10U);
    std::cout<<"case=particle-plasma-chemistry seconds="<<elapsed(start)
             <<" ne="<<state.number_density_m3[0]<<" ni="<<state.number_density_m3[2]
             <<" charge_before="<<diagnostics.charge_density_before_c_m3
             <<" charge_after="<<diagnostics.charge_density_after_c_m3
             <<" max_relative_change="<<diagnostics.max_relative_density_change<<'\n';
    return state.number_density_m3[0]>1.0e15&&std::abs(diagnostics.charge_density_after_c_m3-diagnostics.charge_density_before_c_m3)<1.0e-12?0:1;
}

int particle_breakdown_threshold(){
    using namespace cfd::particle;
    PaschenGas air;air.townsend_a_per_m_pa=112.5;air.townsend_b_v_per_m_pa=2737.0;air.secondary_emission_yield=0.01;
    const auto start=Clock::now();const auto gas=estimate_gas_breakdown_threshold(air,101325.0,1.0e-3,4.0e6);
    const auto multipactor=estimate_parallel_plate_multipactor(2.0e9,1.0e-3,1U,9.1093837015e-31,-1.602176634e-19,1.4);
    std::cout<<"case=particle-breakdown-threshold seconds="<<elapsed(start)
             <<" paschen_voltage_v="<<gas.breakdown_voltage_v
             <<" paschen_field_v_m="<<gas.breakdown_field_v_m
             <<" multipactor_field_v_m="<<multipactor.resonant_field_v_m
             <<" multipactor_impact_ev="<<multipactor.impact_energy_ev
             <<" multipactor_sustains="<<multipactor.secondary_yield_sustains<<'\n';
    return gas.breakdown_voltage_v>0.0&&multipactor.resonant_field_v_m>0.0&&multipactor.secondary_yield_sustains?0:1;
}

int bioheat_sar(){
    using namespace cfd::multiphysics;PennesBioheat2DConfig config;config.nx=8U;config.ny=8U;config.dt_s=1.0;config.blood_perfusion_per_s=0.01;config.blood_temperature_k=310.0;PennesBioheat2D solver(config);solver.initialize(300.0);const double sar=sar_from_rms_electric_field(1.0,1000.0,100.0);solver.set_sar(sar);const auto start=Clock::now();solver.step();double min_t=solver.temperature_k().front(),max_t=min_t;for(double t:solver.temperature_k()){min_t=std::min(min_t,t);max_t=std::max(max_t,t);}std::cout<<"case=bioheat-sar seconds="<<elapsed(start)<<" sar_w_per_kg="<<sar<<" temperature_min_k="<<min_t<<" temperature_max_k="<<max_t<<" cg_iterations="<<solver.linear_result().iterations<<'\n';return solver.linear_result().converged&&max_t>300.0?0:1;
}

int rf_dipole(){
    const double frequency=1.0e9;
    cfd::rf::SinusoidalDipoleSolver dipole(0.5*299792458.0/frequency,frequency);
    const double resistance=dipole.radiation_resistance();
    const double directivity=dipole.directivity();
    std::cout<<"case=rf-dipole frequency_hz="<<frequency<<" radiation_resistance_ohm="<<resistance
             <<" directivity="<<directivity<<"\n";
    return std::abs(resistance-73.13)<0.5&&std::abs(directivity-1.64)<0.03?0:1;
}
int rf_microstrip(){
    const auto line=cfd::rf::microstrip_quasi_static(2.0e-3,1.0e-3,4.4,2.4e9);
    std::cout<<"case=rf-microstrip z0_ohm="<<line.characteristic_impedance
             <<" eps_eff="<<line.effective_permittivity<<" guided_wavelength_m="<<line.guided_wavelength<<"\n";
    return line.characteristic_impedance>40.0&&line.characteristic_impedance<60.0?0:1;
}
int spice_rc(){
    cfd::circuit::Circuit circuit;const auto input=circuit.node("in"),output=circuit.node("out");
    circuit.add_voltage_source("V1",input,0,0.0,{1.0,0.0});circuit.add_resistor("R1",input,output,1000.0);
    circuit.add_capacitor("C1",output,0,1.0e-6);
    const double fc=1.0/(2.0*std::numbers::pi*1.0e3*1.0e-6);const auto ac=circuit.ac(fc);const auto h=circuit.voltage(ac,"out");
    std::cout<<"case=spice-rc frequency_hz="<<fc<<" magnitude="<<std::abs(h)<<" phase_deg="<<std::arg(h)*180.0/std::numbers::pi<<"\n";
    return std::abs(std::abs(h)-1.0/std::sqrt(2.0))<1.0e-5?0:1;
}
int spice_diode(){
    cfd::circuit::Circuit circuit;const auto input=circuit.node("in"),output=circuit.node("out");
    circuit.add_voltage_source("V1",input,0,5.0);circuit.add_resistor("R1",input,output,1000.0);
    circuit.add_diode_model("DDEFAULT",{1.0e-12,1.0,300.0});circuit.add_diode("D1",output,0,"DDEFAULT");
    const auto op=circuit.dc_operating_point();const double v=circuit.voltage(op,"out");
    std::cout<<"case=spice-diode converged="<<op.converged<<" iterations="<<op.iterations<<" diode_voltage="<<v<<"\n";
    return op.converged&&v>0.45&&v<0.9?0:1;
}
int spice_adaptive(){
    using namespace cfd::circuit;
    Circuit circuit;const auto input=circuit.node("in"),output=circuit.node("out");
    circuit.add_voltage_source("VSTEP",input,0,1.0);circuit.add_resistor("R",input,output,1000.0);
    circuit.add_capacitor("C",output,0,1.0e-6);
    AdaptiveTransientConfig control;control.initial_step=1.0e-3;control.maximum_step=1.0e-3;
    control.minimum_step=1.0e-8;control.relative_tolerance=2.0e-5;control.absolute_tolerance=1.0e-8;
    const auto result=circuit.transient_adaptive(5.0e-3,control);
    const double exact=1.0-std::exp(-5.0),error=std::abs(result.points.back().node_voltage[output]-exact);
    control.method=TransientMethod::bdf2;
    const auto gear=circuit.transient_adaptive(5.0e-3,control);
    const double gear_error=std::abs(gear.points.back().node_voltage[output]-exact);
    std::cout<<"case=spice-adaptive accepted="<<result.accepted_steps<<" rejected="<<result.rejected_steps
        <<" dt_min="<<result.minimum_accepted_step<<" dt_max="<<result.maximum_accepted_step
        <<" final_error="<<error<<" bdf2_accepted="<<gear.accepted_steps
        <<" bdf2_rejected="<<gear.rejected_steps<<" bdf2_error="<<gear_error<<'\n';
    return result.rejected_steps>0U&&result.minimum_accepted_step<result.maximum_accepted_step
        &&error<3.0e-4&&gear.accepted_steps>1U&&gear_error<3.0e-4?0:1;
}
int spice_pss_pz(){
    using namespace cfd::circuit;
    constexpr double resistance=1000.0,capacitance=1.0e-6,frequency=100.0;
    Circuit periodic;const auto input=periodic.node("in"),output=periodic.node("out");
    periodic.add_voltage_source("VS",input,0,0.0,{1.0,0.0},SineWaveform{0.0,1.0,frequency,0.0,0.0,0.0});
    periodic.add_resistor("R",input,output,resistance);periodic.add_capacitor("C",output,0,capacitance);
    PeriodicSteadyStateConfig pss;pss.period_s=1.0/frequency;pss.samples_per_period=128U;pss.max_periods=20U;
    pss.relative_tolerance=2.0e-4;pss.absolute_tolerance=1.0e-7;
    const auto steady=periodic_steady_state(periodic,pss);const auto fundamental=fourier_measurement(steady.period,output,frequency);
    PoleZeroConfig pz;pz.start_hz=1.0;pz.stop_hz=1.0e5;pz.samples=40U;pz.denominator_order=1U;pz.numerator_order=0U;
    const auto roots=pole_zero_analysis(periodic,"out",pz);
    const double pole=roots.poles_rad_per_s.empty()?0.0:roots.poles_rad_per_s.front().real();
    std::cout<<"case=spice-pss-pz converged="<<steady.converged<<" periods="<<steady.periods
        <<" pss_residual="<<steady.normalized_residual<<" fundamental="<<fundamental.magnitude
        <<" pole_rad_per_s="<<pole<<" fit_error="<<roots.relative_rms_fit_error<<'\n';
    return steady.converged&&roots.poles_rad_per_s.size()==1U&&std::abs((pole+1000.0)/1000.0)<2.0e-4?0:1;
}
int rf_multiwire(){
    using namespace cfd::rf;
    constexpr double frequency=3.0e8;const double wavelength=299792458.0/frequency;
    ParallelWireMomConfig config;config.frequency_hz=frequency;config.quadrature_order=8U;
    config.wires={{-0.05*wavelength,0.0,0.0,0.47*wavelength,0.001*wavelength,21U},
                  { 0.05*wavelength,0.0,0.0,0.47*wavelength,0.001*wavelength,21U}};
    config.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    const auto open=solve_parallel_thin_wires(config);double driven=0.0,induced=0.0;
    for(std::size_t i=0;i<21U;++i)driven=std::max(driven,std::abs(open.current_a[i]));
    for(std::size_t i=21U;i<42U;++i)induced=std::max(induced,std::abs(open.current_a[i]));
    config.loads={{0U,10U,{1.0e6,0.0}}};const auto loaded=solve_parallel_thin_wires(config);
    const double inv=1.0/std::sqrt(2.0);OrientedWireMomConfig rotated;rotated.frequency_hz=frequency;rotated.quadrature_order=8U;
    for(double sign:{-1.0,1.0}){const WirePoint3 center{0.0,sign*0.05*wavelength,0.0};
        rotated.wires.push_back({{center.x-0.47*wavelength*0.5*inv,center.y,center.z-0.47*wavelength*0.5*inv},
                                 {center.x+0.47*wavelength*0.5*inv,center.y,center.z+0.47*wavelength*0.5*inv},
                                 0.001*wavelength,21U});}
    rotated.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};const auto arbitrary=solve_oriented_thin_wires(rotated);
    double rotation_error=0.0;for(std::size_t i=0;i<open.current_a.size();++i)rotation_error=std::max(rotation_error,std::abs(open.current_a[i]-arbitrary.current_a[i]));
    const double relative_rotation=rotation_error/std::max(driven,1.0e-30);
    std::cout<<"case=rf-multiwire coupling_ratio="<<induced/std::max(driven,1.0e-30)
        <<" unloaded_feed_current="<<open.feed_current_a[0]<<" loaded_feed_current="<<loaded.feed_current_a[0]
        <<" rotation_relative_error="<<relative_rotation<<'\n';
    return induced>driven*1.0e-6&&loaded.feed_current_a[0]<open.feed_current_a[0]&&relative_rotation<2.0e-10?0:1;
}

}

int run_continuation_case(std::string_view name){
    if(name=="lbm-trt")return trt();
    if(name=="fdtd-pml1d")return fdtd_pml();
    if(name=="chemistry-reactor")return reactor();
    if(name=="chemistry-equilibrium")return equilibrium();
    if(name=="optics-gaussian")return gaussian();
    if(name=="multiphysics-thermoelastic")return thermal_structure();
    if(name=="fvm-workspace")return fvm_workspace();
    if(name=="em-frequency1d")return em_frequency1d();
    if(name=="em-edge2d")return em_edge2d();
    if(name=="em-edge3d")return em_edge3d();
    if(name=="em-waveport")return em_waveport();
    if(name=="particle-pic1d")return particle_pic1d();
    if(name=="particle-em-pic1d")return particle_em_pic1d();
    if(name=="particle-pic2d")return particle_pic2d();
    if(name=="particle-em-pic2d")return particle_em_pic2d();
    if(name=="particle-staggered-em-pic2d")return particle_staggered_em_pic2d();
    if(name=="particle-pic3d")return particle_pic3d();
    if(name=="particle-em-pic3d")return particle_em_pic3d();
    if(name=="particle-staggered-em-pic3d")return particle_staggered_em_pic3d();
    if(name=="particle-local-current3d")return particle_local_current3d();
    if(name=="particle-sort-halo3d")return particle_sort_halo3d();
    if(name=="particle-domain-exchange3d")return particle_domain_exchange3d();
    if(name=="particle-comm-exchange3d")return particle_comm_exchange3d();
    if(name=="particle-transport3d")return particle_transport3d();
    if(name=="particle-serialized-transport3d")return particle_serialized_transport3d();
    if(name=="particle-distributed-round3d")return particle_distributed_round3d();
    if(name=="particle-field-guards3d")return particle_field_guards3d();
    if(name=="particle-distributed-step3d")return particle_distributed_step3d();
    if(name=="particle-geant4-transport")return particle_geant4_transport();
    if(name=="particle-transport-dose-bvh")return particle_transport_dose_bvh();
    if(name=="particle-campaign-doe")return particle_campaign_doe();
    if(name=="particle-campaign-execution")return particle_campaign_execution();
    if(name=="particle-campaign-local-runner")return particle_campaign_local_runner();
    if(name=="particle-campaign-control")return particle_campaign_control();
    if(name=="particle-multiserver-deploy")return particle_multiserver_deploy();
    if(name=="particle-multiserver-execution")return particle_multiserver_execution();
    if(name=="particle-multiserver-supervision")return particle_multiserver_supervision();
    if(name=="particle-multiserver-controller")return particle_multiserver_controller();
    if(name=="particle-multiserver-dashboard")return particle_multiserver_dashboard();
    if(name=="particle-plasma-chemistry")return particle_plasma_chemistry();
    if(name=="particle-breakdown-threshold")return particle_breakdown_threshold();
    if(name=="bioheat-sar")return bioheat_sar();
    if(name=="rf-dipole")return rf_dipole();
    if(name=="rf-microstrip")return rf_microstrip();
    if(name=="spice-rc")return spice_rc();
    if(name=="spice-diode")return spice_diode();
    if(name=="spice-adaptive")return spice_adaptive();
    if(name=="spice-pss-pz")return spice_pss_pz();
    if(name=="rf-multiwire")return rf_multiwire();
    return -1;
}
