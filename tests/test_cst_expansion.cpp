#include "cfd/em/frequency_domain.hpp"
#include "cfd/em/edge_fem2d.hpp"
#include "cfd/em/edge_fem3d.hpp"
#include "cfd/em/adaptivity3d.hpp"
#include "cfd/em/wave_port.hpp"
#include "cfd/em/cable.hpp"
#include "cfd/rf/nport.hpp"
#include "cfd/fem/mesh3d.hpp"
#include "cfd/solvers/fem/magnetostatics2d.hpp"
#include "cfd/fem/mesh2d.hpp"
#include "cfd/multiphysics/bioheat.hpp"
#include "cfd/particle/electromagnetic.hpp"
#include "cfd/particle/wakefield.hpp"

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


double edge_maxwell_3d_error(std::size_t cells){
    using namespace cfd::em;
    const auto mesh=cfd::fem::make_box_tet_mesh(cells,cells,cells,1.0,0.8,1.0);
    EdgeMaxwell3DConfig config;config.frequency_hz=8.0e7;config.linear.relative_tolerance=1.0e-8;config.linear.max_iterations=12000U;config.linear.gmres_restart=300U;
    constexpr double mu0=1.25663706212e-6,epsilon0=8.8541878128e-12;
    const double omega=2.0*std::numbers::pi*config.frequency_hz;
    const double lambda=2.0*std::numbers::pi*std::numbers::pi;
    const double coefficient=lambda/(mu0*config.relative_permeability)-omega*omega*epsilon0*config.relative_permittivity;
    const auto solution=solve_driven_edge_maxwell_3d(mesh,config,[=](cfd::fem::Point3 point)->ComplexVec3{
        const double ey=std::sin(std::numbers::pi*point.x)*std::sin(std::numbers::pi*point.z);
        return {Complex{},Complex{0.0,coefficient*ey/omega},Complex{}};
    });
    require(solution.linear_result.converged,"3-D edge-Maxwell sparse solve converges");
    double error2=0.0,reference2=0.0;
    for(std::size_t element=0;element<mesh.tetrahedra.size();++element){
        cfd::fem::Point3 c{};for(auto node:mesh.tetrahedra[element].node){c.x+=mesh.nodes[node].x/4.0;c.y+=mesh.nodes[node].y/4.0;c.z+=mesh.nodes[node].z/4.0;}
        (void)c.y;const Complex exact=std::sin(std::numbers::pi*c.x)*std::sin(std::numbers::pi*c.z);const auto numerical=solution.electric_centroid_v_per_m[element];
        error2+=std::norm(numerical[0])+std::norm(numerical[1]-exact)+std::norm(numerical[2]);reference2+=std::norm(exact);
    }
    return std::sqrt(error2/reference2);
}

void edge_element_maxwell_3d_and_ports(){
    using namespace cfd::em;
    const double coarse=edge_maxwell_3d_error(2U),fine=edge_maxwell_3d_error(3U);
    require(coarse<0.55,"3-D Nedelec Maxwell manufactured coarse accuracy");
    require(fine<coarse*0.82,"3-D Nedelec Maxwell mesh refinement");

    const auto mesh=cfd::fem::make_box_tet_mesh(2U,2U,2U,1.0,0.8,1.0);EdgeMaxwell3DConfig config;config.frequency_hz=8.0e7;config.linear.relative_tolerance=1.0e-8;config.linear.max_iterations=12000U;config.linear.gmres_restart=300U;
    const auto probe=solve_driven_edge_maxwell_3d(mesh,config,[](cfd::fem::Point3)->ComplexVec3{return {Complex{},Complex{0.0,1.0},Complex{}};});
    std::size_t free_edge=probe.edges.size();for(std::size_t e=0;e<probe.edges.size();++e)if(!probe.boundary_edge[e]){free_edge=e;break;}
    require(free_edge<probe.edges.size(),"3-D Maxwell mesh exposes an interior edge port");
    const EdgeCurrentPort3D port{probe.edges[free_edge][0],probe.edges[free_edge][1],Complex{1.0e-6,0.0}};
    const auto port_solution=solve_driven_edge_maxwell_3d(mesh,config,{},std::span<const EdgeCurrentPort3D>(&port,1U));
    double peak=0.0;for(const auto v:port_solution.edge_voltage_v)peak=std::max(peak,std::abs(v));require(peak>0.0,"3-D edge-current port excites the sparse Maxwell system");

    const auto q1=resonator_quality_3d(mesh,config,probe,5.8e7);const auto q4=resonator_quality_3d(mesh,config,probe,4.0*5.8e7);
    require(q1.conductor_loss_w>0.0&&std::isfinite(q1.quality_factor)&&q1.quality_factor>0.0,"3-D resonator conductor-loss Q is finite");
    require(std::abs(q4.quality_factor/q1.quality_factor-2.0)<1.0e-10,"surface-resistance Q scales with sqrt(conductivity)");
}


void adaptive_edge_mesh_3d(){
    using namespace cfd::em;const auto mesh=cfd::fem::make_box_tet_mesh(2U,2U,2U,1.0,0.8,1.0);EdgeMaxwell3DConfig config;config.frequency_hz=8.0e7;config.linear.relative_tolerance=1.0e-8;config.linear.gmres_restart=300U;
    constexpr double mu0=1.25663706212e-6,epsilon0=8.8541878128e-12;const double omega=2.0*std::numbers::pi*config.frequency_hz,coefficient=2.0*std::numbers::pi*std::numbers::pi/mu0-omega*omega*epsilon0;
    const auto source=[=](cfd::fem::Point3 p)->ComplexVec3{return {Complex{},Complex{0.0,coefficient*std::sin(std::numbers::pi*p.x)*std::sin(std::numbers::pi*p.z)/omega},Complex{}};};
    const auto field=solve_driven_edge_maxwell_3d(mesh,config,source);const auto indicators=maxwell_face_jump_indicators_3d(mesh,field);double peak=0.0;for(double value:indicators)peak=std::max(peak,value);require(peak>0.0,"3-D Maxwell jump estimator detects nonuniform field");const auto marked=mark_maxwell_dorfler(indicators,0.45);std::size_t count=0U;for(auto value:marked)count+=value!=0U;require(count>0U&&count<marked.size(),"Dorfler marking selects a local Tet4 subset");const auto refined=refine_tet4_marked_longest_edges(mesh,marked);require(refined.element_count()>mesh.element_count(),"marked longest-edge refinement adds Tet4 elements");require(refined.element_count()<2U*mesh.element_count(),"local Tet4 refinement remains below global bisection size");const auto refined_field=solve_driven_edge_maxwell_3d(refined,config,source);require(refined_field.linear_result.converged,"refined 3-D Maxwell system converges");
}

void wave_port_modes(){
    using namespace cfd::em;constexpr double c0=299792458.0;const double a=22.86e-3,b=10.16e-3;
    const auto modes=rectangular_waveguide_modes(a,b,10.0e9,1.0,1.0,5U);require(!modes.empty(),"rectangular wave-port modes returned");
    const auto& te10=modes.front();require(te10.family==RectangularWaveguideMode::Family::te&&te10.m==1U&&te10.n==0U,"WR90 fundamental is TE10");
    const double exact_fc=c0/(2.0*a);require(std::abs(te10.cutoff_frequency_hz-exact_fc)/exact_fc<1.0e-9,"TE10 analytical cutoff frequency");require(te10.propagating,"TE10 propagates at 10 GHz");
    require(std::abs(rectangular_mode_power_w(te10,a,b,96U,64U)-1.0)<2.0e-12,"propagating wave-port mode normalized to one watt");
    bool found_evanescent=false;for(const auto& mode:modes)if(!mode.propagating)found_evanescent=true;require(found_evanescent,"wave-port extraction retains evanescent higher modes");
}


void adaptive_rational_rf_model(){
    using namespace cfd::rf;
    const double fp=1.0e8;const std::array<double,1> z0{{50.0}};
    AdaptiveNetworkSweepConfig config;config.initial_points=5U;config.maximum_points=65U;config.maximum_refinements=8U;config.relative_tolerance=2.0e-3;config.absolute_tolerance=1.0e-7;
    const auto sweep=adaptive_network_sweep(1.0e6,1.0e9,z0,[=](double frequency){
        ComplexMatrix s(1U);const Complex jw{0.0,frequency/fp};s(0,0)=1.0/(1.0+jw);return s;
    },config);
    require(sweep.evaluations>=config.initial_points&&sweep.points.size()>=config.initial_points,"adaptive RF sweep evaluates/refines samples");
    const auto model=fit_rational_network(sweep.points,4U);for(double frequency:{3.0e6,3.0e7,3.0e8}){const Complex exact=1.0/(1.0+Complex{0.0,frequency/fp});const Complex fitted=model.evaluate(frequency)(0,0);require(std::abs(fitted-exact)<2.0e-7,"stable rational RF reduced-order model reproduces one-pole network");}
    for(double pole:model.poles_rad_per_s)require(pole<0.0,"rational RF model poles are stable");
}


void cable_harness_and_shield(){
    using namespace cfd::em;using cfd::rf::ComplexMatrix;
    MulticonductorRlcg line;line.conductors=1U;line.resistance={0.0};line.inductance={250.0e-9};line.conductance={0.0};line.capacitance={100.0e-12};ComplexMatrix load(1U);load(0,0)=50.0;const std::array<Complex,1> drive{{1.0}};const auto result=solve_multiconductor_cable(0.1,400U,1.0e9,line,drive,load);require(std::abs(result.input_current_a[0]-Complex{0.02,0.0})<2.0e-6,"matched RLCG cable input current equals V/Z0");require(std::abs(std::abs(result.load_voltage_v[0])-1.0)<2.0e-5,"matched lossless cable preserves voltage magnitude");require(std::abs(result.load_voltage_v[0]+Complex{1.0,0.0})<4.0e-5,"matched line accumulates analytical pi phase");
    MulticonductorRlcg pair;pair.conductors=2U;pair.resistance.assign(4U,0.0);pair.conductance.assign(4U,0.0);pair.inductance={250e-9,50e-9,50e-9,250e-9};pair.capacitance={100e-12,-20e-12,-20e-12,100e-12};ComplexMatrix loads(2U);loads(0,0)=50.0;loads(1,1)=50.0;const std::array<Complex,2> pair_drive{{1.0,0.0}};const auto coupled=solve_multiconductor_cable(0.2,500U,200.0e6,pair,pair_drive,loads);require(std::abs(coupled.load_voltage_v[1])>1.0e-3,"coupled multiconductor RLCG line produces crosstalk on quiet conductor");
    CableShieldTransferModel shield;shield.dc_transfer_resistance_ohm_per_m=5.0e-3;shield.transfer_inductance_h_per_m=2.0e-9;shield.skin_corner_hz=1.0e6;const auto low=cable_shield_transfer_impedance_ohm_per_m(1.0e3,shield),high=cable_shield_transfer_impedance_ohm_per_m(100.0e6,shield);require(std::abs(high)>std::abs(low),"shield transfer impedance includes skin/inductive frequency dependence");require(std::abs(cable_shield_induced_voltage_v(1.0e6,Complex{2.0,0.0},3.0,shield)-6.0*cable_shield_transfer_impedance_ohm_per_m(1.0e6,shield))<1.0e-14,"shield transfer impedance maps shield current to induced cable voltage");require(field_to_cable_open_circuit_voltage_v(Complex{10.0,0.0},0.25,Complex{0.5,0.0})==Complex{1.25,0.0},"external field effective-height cable coupling");
}


void nonlinear_magnetics_coils_and_force(){
    using namespace cfd::fem;constexpr double mu0=1.25663706212e-6;
    PiecewiseLinearBHCurve curve({0.0,0.5,1.0,1.5},{0.0,0.5/(mu0*1000.0),1.0/(mu0*500.0),1.5/(mu0*100.0)});require(curve.reluctivity(0.25)>0.0&&curve.reluctivity(1.4)>curve.reluctivity(0.25),"piecewise B-H law increases reluctivity toward saturation");
    const auto mesh=make_rectangle_tri_mesh(10U,8U,1.0,0.8);Magnetostatics2D solver(mesh);for(int patch=0;patch<4;++patch)solver.set_boundary(patch,{ScalarBoundaryType::dirichlet,[](Node2){return 0.0;},{}});RectangularStrandedCoil2D coil{0.35,0.65,0.25,0.55,100.0,2.0};const auto result=solver.solve_nonlinear([&](double b){return curve.reluctivity(b);},[&](Node2 p){return coil.current_density_a_per_m2(p);},40U,2.0e-5,0.6);require(result.converged,"nonlinear B-H magnetostatic Picard solve converges");const double linkage=coil_flux_linkage_wb_turn(mesh,solver.vector_potential(),coil),inductance=coil_inductance_h(mesh,solver.vector_potential(),coil);require(std::isfinite(linkage)&&std::isfinite(inductance)&&inductance>0.0,"stranded coil exports positive flux linkage/inductance for circuit coupling");
    const auto stress_mesh=make_rectangle_tri_mesh(2U,2U,2.0,1.0);std::vector<MagneticFluxDensity2> uniform(stress_mesh.element_count(),{1.0,0.0});const auto force=maxwell_stress_boundary_force_2d(stress_mesh,uniform,1,2.0,{0.0,0.0});require(std::abs(force.force_x_n_per_m-0.25)<1.0e-12&&std::abs(force.force_y_n_per_m)<1.0e-12,"Maxwell stress recovers uniform-field boundary traction");require(std::abs(force.torque_z_n+0.125)<1.0e-12,"Maxwell stress torque about origin matches boundary resultant");
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


void relativistic_particles_boundaries_and_wakes(){
    using namespace cfd::particle;constexpr double c0=299792458.0;
    ChargedParticle relativistic;relativistic.velocity_m_per_s={0.8*c0,0.0,0.0};relativistic.charge_c=1.0;relativistic.mass_kg=1.0;const double gamma0=lorentz_gamma(relativistic);
    for(std::size_t i=0;i<10000U;++i)relativistic_boris_push(relativistic,{{0,0,0},{0,0,1}},1.0e-3);
    require(std::abs(lorentz_gamma(relativistic)-gamma0)<2.0e-11,"relativistic Boris magnetic field preserves gamma");

    AxisAlignedParticleBox box{{0,0,0},{1,1,1}};ChargedParticle reflected;reflected.position_m={1.02,0.5,0.5};reflected.velocity_m_per_s={2.0,0.0,0.0};const auto reflected_result=apply_particle_box_boundary(reflected,box,ParticleWallMode::specular_reflect);require(reflected_result.alive&&reflected_result.impacted&&reflected.velocity_m_per_s.x==-2.0&&std::abs(reflected.position_m.x-0.98)<1.0e-14,"specular particle wall mirrors position/normal velocity");
    constexpr double me=9.1093837139e-31,qe=1.602176634e-19;ChargedParticle absorbed;absorbed.position_m={-1.0e-6,0.5,0.5};absorbed.mass_kg=me;absorbed.charge_c=-qe;absorbed.weight=2.0;absorbed.velocity_m_per_s={std::sqrt(2.0*300.0*qe/me),0.0,0.0};SecondaryEmissionModel secondary;const auto absorbed_result=apply_particle_box_boundary(absorbed,box,ParticleWallMode::absorb,&secondary);require(!absorbed_result.alive&&absorbed_result.secondary_macro_weight>2.9&&absorbed_result.secondary_macro_weight<3.1,"absorbing wall reports secondary-electron macro weight near yield maximum");

    constexpr double fr=1.0e9,rs=1000.0,q=10.0;const double ds=c0/(100.0*fr);std::vector<double> wake(8192U);for(std::size_t i=0;i<wake.size();++i)wake[i]=resonator_longitudinal_wake_v_per_c(static_cast<double>(i)*ds,rs,fr,q);const auto zres=wake_impedance_ohm(wake,ds,fr),zlow=wake_impedance_ohm(wake,ds,0.5*fr),zhigh=wake_impedance_ohm(wake,ds,2.0*fr);require(std::abs(zres)>3.0*std::abs(zlow)&&std::abs(zres)>3.0*std::abs(zhigh),"wake impedance peaks near resonator frequency");std::vector<double> line(64U,0.0);line[0]=1.0e-9/ds;const auto potential=bunch_wake_potential_v(wake,line,ds);require(std::abs(potential[0]+wake[0]*1.0e-9)<1.0e-8*std::max(1.0,std::abs(potential[0])),"wake convolution preserves point-bunch normalization");
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
    try{frequency_domain_maxwell();edge_element_maxwell();edge_element_maxwell_3d_and_ports();adaptive_edge_mesh_3d();wave_port_modes();adaptive_rational_rf_model();cable_harness_and_shield();nonlinear_magnetics_coils_and_force();charged_particles_and_pic();relativistic_particles_boundaries_and_wakes();bioheat();std::cout<<"CST expansion tests passed\n";return 0;}
    catch(const std::exception& error){std::cerr<<"CST expansion test failure: "<<error.what()<<'\n';return 1;}
}
