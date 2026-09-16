#include "continuation_cases.hpp"
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
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
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
}

int run_continuation_case(std::string_view name){
    if(name=="lbm-trt")return trt();
    if(name=="fdtd-pml1d")return fdtd_pml();
    if(name=="chemistry-reactor")return reactor();
    if(name=="chemistry-equilibrium")return equilibrium();
    if(name=="optics-gaussian")return gaussian();
    if(name=="multiphysics-thermoelastic")return thermal_structure();
    if(name=="fvm-workspace")return fvm_workspace();
    return -1;
}
