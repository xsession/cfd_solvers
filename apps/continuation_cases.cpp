#include "continuation_cases.hpp"
#include <numbers>
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
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
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
    if(name=="rf-dipole")return rf_dipole();
    if(name=="rf-microstrip")return rf_microstrip();
    if(name=="spice-rc")return spice_rc();
    if(name=="spice-diode")return spice_diode();
    if(name=="spice-adaptive")return spice_adaptive();
    if(name=="spice-pss-pz")return spice_pss_pz();
    if(name=="rf-multiwire")return rf_multiwire();
    return -1;
}
