#include "cfd/chemistry/aqueous_equilibrium.hpp"
#include "cfd/chemistry/implicit_reactor.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/multiphysics/electro_thermal.hpp"
#include "cfd/multiphysics/partitioned.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fdtd/port.hpp"
#include "cfd/solvers/fem/heat2d.hpp"
#include "cfd/solvers/fem/elasticity2d.hpp"
#include "cfd/solvers/fvm/scalar_transport.hpp"
#include "cfd/solvers/lbm/one_step_pull.hpp"
#include "cfd/solvers/optics/gaussian_beam.hpp"
#include "cfd/solvers/optics/sequential.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
template<class F> void rejects(F&& f,const char* message){
    bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,message);
}

void heat_and_coupling(){
    // T=300+t is spatially uniform: rho*c*dT/dt=2, including ramped boundaries.
    cfd::fem::Heat2DConfig config;config.volumetric_heat_capacity=2.0;config.dt=0.01;
    cfd::fem::Heat2D heat(cfd::fem::make_rectangle_tri_mesh(8,7),config);
    heat.initialize([](auto){return 300.0;});
    heat.set_dirichlet([](auto,double t){return 300.0+t;});
    const std::vector<double> source(heat.mesh().element_count(),2.0);
    for(int i=0;i<5;++i)heat.step_element_source(source);
    for(double t:heat.temperature())check(std::abs(t-300.05)<1.0e-7,"heat must retain old-time boundary mass contribution");

    cfd::multiphysics::JouleHeatingCoupler2D coupled(cfd::fem::make_rectangle_tri_mesh(8,7),config);
    rejects([&]{coupled.thermal_step();},"unsolved electrical state must not silently heat with zeros");
    using Type=cfd::fem::ScalarBoundaryType;
    coupled.set_electrical_boundary(0,{Type::dirichlet,[](auto){return 0.0;},{}});
    coupled.set_electrical_boundary(1,{Type::dirichlet,[](auto){return 1.0;},{}});
    coupled.initialize_temperature([](auto){return 300.0;});
    coupled.set_thermal_dirichlet([](auto,double t){return 300.0+t;});
    coupled.solve_electrical([](auto){return 2.0;});
    coupled.thermal_run(5);
    for(double t:coupled.thermal().temperature())check(std::abs(t-300.05)<1.0e-7,"Joule-heating analytical temperature ramp");
    for(auto mode:{cfd::fem::ElasticityMode2D::planeStress,cfd::fem::ElasticityMode2D::planeStrain}){
        cfd::fem::Elasticity2DConfig mechanics;mechanics.mode=mode;
        cfd::fem::Elasticity2D solid(coupled.thermal().mesh(),mechanics);
        solid.set_dirichlet([](auto x){return x.x==0.0;},[](auto){return cfd::fem::Displacement2{};},true,false);
        solid.set_dirichlet([](auto x){return x.y==0.0;},[](auto){return cfd::fem::Displacement2{};},false,true);
        solid.solve_thermal(coupled.thermal().temperature(),1.0e-3,300.0);
        const double strain=5.0e-5*(mode==cfd::fem::ElasticityMode2D::planeStrain?1.3:1.0);
        for(std::size_t i=0;i<solid.mesh().node_count();++i){
            check(std::abs(solid.displacement()[i].x-strain*solid.mesh().nodes[i].x)<1.0e-8
                &&std::abs(solid.displacement()[i].y-strain*solid.mesh().nodes[i].y)<1.0e-8,
                "electrical to thermal to structural analytical expansion");
        }
    }
    coupled.set_electrical_boundary(1,{Type::dirichlet,[](auto){return 2.0;},{}});
    rejects([&]{coupled.thermal_step();},"boundary change must invalidate electrical heating source");

    std::vector<double> state{0.0};
    cfd::multiphysics::PartitionedCouplerConfig control;control.max_iterations=1;control.initial_relaxation=1.0;
    const auto result=cfd::multiphysics::solve_partitioned_fixed_point(state,[](auto,auto out){out[0]=2.0;},control);
    check(result.converged&&result.iterations==1&&result.residual_rms==0.0,"fixed-point final update must be checked");
}

void fdtd_absorption_and_port(){
    using namespace cfd::fdtd;
    Maxwell1DConfig config{400,1.0e-4,0.7,1.0,1.0,Boundary1D::pml};config.pml_cells=32;
    Maxwell1D pml(config);config.boundary=Boundary1D::pec;Maxwell1D pec(config);
    pml.initialize_gaussian(0.5,0.025);pec.initialize_gaussian(0.5,0.025);
    const double initial=pec.energy();pml.step(800);pec.step(800);
    check(pml.energy()<initial*1.0e-4,"PML must absorb outgoing pulse");
    check(pec.energy()>initial*0.8,"PEC enclosure must retain the pulse");
    rejects([&]{pml.set_debye_material(0,40,2.0,1.0,1.0e-12);},"dispersive material cannot overlap baseline PML");
    rejects([&]{pml.set_material(50,100,0.01);},"material replacement must respect fixed CFL limit");
    rejects([]{Maxwell1D bad({0});},"zero cells must be rejected before allocation underflow");
    rejects([&]{pml.set_drude_material(50,100,1.0,10.0/pml.dt(),0.0);},"unstable ADE pole must be rejected");
    config.boundary=Boundary1D::mur1;Maxwell1D wave(config);wave.initialize_gaussian(0.3,0.03);
    constexpr double impedance=376.73031366686166;
    WavePort1D port(1.0e9,impedance);
    double forward=0.0,backward=0.0;
    const std::size_t probe=220;
    for(int i=0;i<250;++i){
        const double old_e=wave.electric()[probe];wave.step();
        // E averaged in time, H averaged in space: both at (x_probe,t-dt/2).
        const auto a=port.instantaneous(0.5*(old_e+wave.electric()[probe]),
            0.5*(wave.magnetic()[probe-1]+wave.magnetic()[probe]));
        forward+=a.forward*a.forward;backward+=a.backward*a.backward;
    }
    check(forward>1.0&&backward<forward*1.0e-4,"physical right-travelling pulse must be classified as forward");
    wave.initialize_gaussian(0.3,0.03);check(wave.time()==0.0,"Gaussian reinitialization resets time");
    for(double h:wave.magnetic())check(h==0.0,"Gaussian reinitialization clears H");
}

void fdtd_ade_response(){
    using namespace cfd::fdtd;
    Maxwell1D reference({512,1.0e-4,0.4});const double dt=reference.dt();
    for(auto model:{DispersionModel1D::debye,DispersionModel1D::drude,DispersionModel1D::lorentz}){
        double previous_error=0.0;
        for(int refinement=1;refinement<=2;++refinement){
            Maxwell1D solver({512,1.0e-4,0.4/static_cast<double>(refinement)});
            const double tau=40.0*dt,omega=0.01/dt;
            if(model==DispersionModel1D::debye)solver.set_debye_material(0,512,2.0,3.0,tau);
            if(model==DispersionModel1D::drude)solver.set_drude_material(0,512,2.0,omega,0.0);
            if(model==DispersionModel1D::lorentz)solver.set_lorentz_material(0,512,2.0,3.0,omega,0.0);
            for(std::size_t i=0;i<512;++i)solver.set_hard_source(i,1.0);
            solver.step(200U*static_cast<std::size_t>(refinement));
            const double t=solver.time();
            double exact=0.4+0.6*std::exp(-t/(tau*0.4));
            if(model==DispersionModel1D::drude)exact=std::cos(omega*t/std::sqrt(2.0));
            if(model==DispersionModel1D::lorentz)exact=0.4+0.6*std::cos(omega*t*std::sqrt(2.5));
            const double error=std::abs(solver.electric()[256]-exact);
            check(error<0.01,"ADE uniform-field analytical relaxation/oscillation");
            if(refinement==2)check(error<previous_error*0.6+1.0e-13,"ADE time refinement must reduce analytical error");
            previous_error=error;
        }
    }
}

void chemistry(){
    using namespace cfd::chemistry;
    ReactionNetwork network({{"A",0.01,0},{"B",0.01,0}});
    ElementaryReaction reaction{{{0,1.0}},{{1,1.0}},{1000.0,0.0,0.0},EquilibriumConstant{4.0,298.15,0.0}};
    network.add_reaction(reaction);
    const std::vector<double> equilibrium{0.2,0.8};
    check(std::abs(network.reaction_rates(equilibrium,298.15)[0])<1.0e-12,"reversible detailed balance");
    std::vector<double> state{1.0,0.0};ReactorConfig control;control.initial_step=0.01;control.relative_tolerance=1.0e-5;
    const auto result=integrate_isothermal(network,state,298.15,0.05,control);
    check(result.accepted_steps>0&&result.rejected_steps>0,"stiff reactor adaptive control exercised");
    check(std::abs(state[0]-0.2)<1.0e-6&&std::abs(state[1]-0.8)<1.0e-6,"stiff reversible reactor reaches equilibrium");
    check(std::abs(state[0]+state[1]-1.0)<1.0e-10,"stiff reaction preserves stoichiometric invariant");
    ReactionNetwork dimer({{"A",0.01,0},{"A2",0.02,0}});
    dimer.add_reaction({{{0,2.0}},{{1,1.0}},{1000.0,0.0,0.0},{}});
    std::vector<double> nonlinear{1.0,0.0};
    (void)integrate_isothermal(dimer,nonlinear,298.15,0.01,control);
    check(std::abs(nonlinear[0]-1.0/21.0)<2.0e-4,"nonlinear stiff dimerization analytical transient");
    check(std::abs(nonlinear[0]+2.0*nonlinear[1]-1.0)<1.0e-9,"nonlinear reaction atom balance");
    ReactorConfig exhausted=control;exhausted.max_steps=1;
    std::vector<double> unchanged{1.0,0.0};
    rejects([&]{(void)integrate_isothermal(dimer,unchanged,298.15,0.01,exhausted);},"reactor step-budget failure is reported");
    check(unchanged[0]==1.0&&unchanged[1]==0.0,"failed reactor integration leaves caller state unchanged");
    std::vector<double> scratch(1),source(2);
    network.source_terms(state,298.15,source,scratch);
    check(std::abs(source[0]+source[1])<1.0e-14,"allocation-free chemistry source conserves mass");
    const EquilibriumConstant endothermic{2.0,300.0,10000.0};
    check(endothermic.value(400.0)>2.0,"van't Hoff endothermic equilibrium trend");

    const auto water=equilibrate_acids({});check(std::abs(water.ph-7.0)<1.0e-12,"pure-water pH");
    const std::vector<AcidFamily> acid{{0.1,0,{1.0e-5}}};
    const auto weak=equilibrate_acids(acid);
    const double exact_h=2.0*1.0e-5*0.1/(1.0e-5+std::sqrt(1.0e-10+4.0e-6));
    check(std::abs(weak.hydrogen_mol_per_litre-exact_h)<1.0e-10,"weak-acid mass-action quadratic");
    const auto buffer=equilibrate_acids(acid,0.05);
    check(std::abs(buffer.ph-5.0)<2.0e-4,"half-neutralized acid buffer pH");
    check(std::abs(buffer.charge_residual_mol_per_litre)<1.0e-13,"speciation electroneutrality");
    check(std::abs(buffer.family_concentrations[0][0]+buffer.family_concentrations[0][1]-0.1)<1.0e-14,"acid family mass balance");
    const auto alkaline=equilibrate_acids({},0.01);check(std::abs(alkaline.ph-12.0)<1.0e-8,"strong-base limit");
}

void fvm_workspaces(){
    using namespace cfd::fvm;
    const auto mesh=make_cartesian_hexa_mesh(12,5,3,1.0,1.0,1.0);
    std::vector<double> values(mesh.cell_count()),boundary(mesh.face_count()),flux(mesh.face_count());
    const auto phi=[](Vec3 x){return 1.0+2.0*x.x-0.5*x.y+0.25*x.z;};
    for(std::size_t i=0;i<values.size();++i)values[i]=phi(mesh.cells()[i].center);
    for(std::size_t f=0;f<flux.size();++f){boundary[f]=phi(mesh.faces()[f].center);flux[f]=mesh.faces()[f].area.x;}
    ScalarSchemeWorkspace workspace;workspace.resize(mesh);
    const auto bytes=workspace.allocated_bytes();
    const auto* gradient_ptr=workspace.gradient.data();
    std::vector<double> faces(mesh.face_count()),divergence(mesh.cell_count());
    for(auto scheme:{FaceInterpolationScheme::linear,FaceInterpolationScheme::upwind,
                    FaceInterpolationScheme::bounded_linear,FaceInterpolationScheme::minmod,FaceInterpolationScheme::van_leer}){
        const auto reference=convective_divergence_scalar(mesh,values,flux,scheme,boundary);
        convective_divergence_scalar_into(mesh,values,flux,divergence,workspace,scheme,boundary);
        for(std::size_t c=0;c<values.size();++c)check(divergence[c]==reference[c],"FVM workspace/reference parity");
        double integral=0.0,boundary_sum=0.0;
        interpolate_scalar_to_faces_into(mesh,values,faces,workspace,scheme,flux,boundary);
        for(std::size_t c=0;c<values.size();++c)integral+=divergence[c]*mesh.cells()[c].volume;
        for(std::size_t f=0;f<faces.size();++f)if(mesh.faces()[f].boundary())boundary_sum+=flux[f]*faces[f];
        check(std::abs(integral-boundary_sum)<1.0e-12,"convection global flux balance");
        if(scheme==FaceInterpolationScheme::minmod||scheme==FaceInterpolationScheme::van_leer)
            for(std::size_t f=0;f<faces.size();++f)if(!mesh.faces()[f].boundary())
                check(std::abs(faces[f]-boundary[f])<1.0e-11,"TVD linear field exactness on orthogonal mesh");
    }
    check(workspace.allocated_bytes()==bytes&&workspace.gradient.data()==gradient_ptr,"FVM workspace reuses allocations");
    for(std::size_t c=0;c<values.size();++c)values[c]=mesh.cells()[c].center.x<0.5?0.0:1.0;
    for(auto scheme:{FaceInterpolationScheme::minmod,FaceInterpolationScheme::van_leer}){
        interpolate_scalar_to_faces_into(mesh,values,faces,workspace,scheme,flux);
        for(double f:faces)check(f>=0.0&&f<=1.0,"TVD jump must not overshoot");
    }
}

void fvm_corrected_operators(){
    using namespace cfd::fvm;
    const auto mesh=make_sheared_cartesian_hexa_mesh(8,7,3,1.0,1.0,1.0,0.35);
    const Vec3 exact{2.0,-0.5,0.25};
    const auto phi=[&](Vec3 p){return 1.0+dot(exact,p);};
    std::vector<double> cells(mesh.cell_count()),boundary(mesh.face_count());
    for(std::size_t c=0;c<cells.size();++c)cells[c]=phi(mesh.cells()[c].center);
    for(std::size_t f=0;f<boundary.size();++f)boundary[f]=phi(mesh.faces()[f].center);
    const auto gradient=least_squares_gradient_scalar(mesh,cells,boundary);
    const auto normal=corrected_face_normal_gradient_scalar(mesh,cells,boundary);
    const auto laplacian=corrected_laplacian_scalar(mesh,cells,1.0,boundary);
    for(const auto& g:gradient)check(magnitude(g-exact)<1.0e-11,"least-squares exact linear gradient on sheared mesh");
    for(std::size_t f=0;f<normal.size();++f)check(std::abs(normal[f]-dot(exact,mesh.faces()[f].area)/magnitude(mesh.faces()[f].area))<1.0e-11,
        "nonorthogonal face-normal gradient of linear field");
    for(double l:laplacian)check(std::abs(l)<1.0e-9,"corrected Laplacian of linear field vanishes");
}

void transport_cache(){
    using namespace cfd::fvm;
    auto mesh=make_cartesian_hexa_mesh(8,4,2,1.0,1.0,1.0);
    ScalarTransport cached(mesh),fresh(mesh);
    cached.initialize([](Vec3 p){return 1.0+p.x;});fresh.initialize([](Vec3 p){return 1.0+p.x;});
    std::vector<double> flux(mesh.face_count(),0.0);
    for(int step=0;step<5;++step){
        fresh.set_face_flux(flux); // force independent reconstruction for comparison
        cached.step();fresh.step();
        for(std::size_t i=0;i<mesh.cell_count();++i)check(cached.values()[i]==fresh.values()[i],"cached scalar operator numerical parity");
    }
    check(cached.operator_assemblies()==1&&fresh.operator_assemblies()==5,"scalar operator/ILU are reused between time steps");
    const double initial=cached.volume_integral();cached.set_source(2.0);cached.step();
    check(cached.operator_assemblies()==1,"source change must preserve cached operator");
    check(std::abs(cached.volume_integral()-initial-0.002)<1.0e-9,"cached scalar source conservation");
    cached.set_boundary(mesh.patches().front().name,ScalarBoundaryType::fixedValue,0.0);cached.step();
    check(cached.operator_assemblies()==2,"boundary change invalidates scalar operator");
    cached.set_face_flux(flux);cached.step();check(cached.operator_assemblies()==3,"flux change invalidates scalar operator");
}

template<class Descriptor> void trt(){
    cfd::lbm::InPlaceLbmConfig config;config.nx=24;config.ny=24;config.nz=24;config.tau=0.8F;
    cfd::lbm::OneStepPullSolver<Descriptor> bgk(config),equal(config,config.tau),different(config,1.1F);
    bgk.initialize_taylor_green();equal.initialize_taylor_green();different.initialize_taylor_green();
    const auto initial=different.compute_macroscopic();
    const double mass=different.mass();bgk.step(10);equal.step(10);different.step(10);
    const auto reference=bgk.compute_macroscopic(),parity=equal.compute_macroscopic(),fields=different.compute_macroscopic();
    check(reference.ux==parity.ux&&reference.rho==parity.rho,"TRT equal relaxation reduces exactly to BGK");
    check(std::abs(different.mass()-mass)/mass<3.0e-6,"TRT conservation");
    double error=0.0;for(std::size_t i=0;i<fields.ux.size();++i)error=std::max(error,std::abs(static_cast<double>(reference.ux[i]-fields.ux[i])));
    std::cout<<"TRT D"<<Descriptor::dimensions<<"Q"<<Descriptor::q<<" BGK max error="<<error<<'\n';
    check(error<5.0e-4,"TRT smooth-flow agreement with BGK viscosity");
    double projection=0.0,norm=0.0;
    for(std::size_t i=0;i<fields.ux.size();++i){
        projection+=static_cast<double>(fields.ux[i])*initial.ux[i]+static_cast<double>(fields.uy[i])*initial.uy[i];
        norm+=static_cast<double>(initial.ux[i])*initial.ux[i]+static_cast<double>(initial.uy[i])*initial.uy[i];
    }
    const double k=2.0*std::numbers::pi/24.0;
    const double decay=std::exp(-Descriptor::dimensions*k*k*(static_cast<double>(config.tau)-0.5)/3.0*10.0);
    check(std::abs(projection/norm-decay)<0.015,"TRT analytical Taylor-Green viscous decay");
    config.acceleration_x=1.0e-5F;
    cfd::lbm::OneStepPullSolver<Descriptor> forced(config,1.1F);forced.step(20);
    for(float u:forced.compute_macroscopic().ux)check(std::abs(u-20.5e-5F)<2.0e-6F,"TRT Guo forcing momentum increment");
}

void optics(){
    using namespace cfd::optics;
    GaussianBeam beam(500.0e-9,1.0e-3);
    const double zr=std::numbers::pi*1.0e-6/500.0e-9;
    beam.propagate(zr);
    check(std::abs(beam.radius()-std::sqrt(2.0)*1.0e-3)<1.0e-15,"Gaussian Rayleigh-range beam radius");
    check(std::abs(beam.curvature_radius()-2.0*zr)<1.0e-12,"Gaussian wavefront curvature");
    beam.propagate(-zr);beam.thin_lens(0.1);
    const double focus=-beam.q().real();beam.propagate(focus);
    const double exact=1.0e-3/std::sqrt(1.0+zr*zr/0.01);
    check(std::abs(beam.radius()-exact)<1.0e-15,"Gaussian thin-lens waist");
    SequentialSurface parabola{SurfaceType::conic,0.0,10.0,2.0,1.0,-1.0};
    check(std::abs(surface_sag(parabola,1.0)-0.05)<1.0e-14,"parabolic conic sag");
    SequentialOpticalSystem system;system.add_surface(parabola);
    const auto hit=system.trace({{1.0,0.0,-1.0},{0.0,0.0,1.0},550.0});
    check(hit.valid&&std::abs(hit.ray.origin.z-0.05)<1.0e-12,"conic axial ray intersection");
    SequentialSurface asphere=parabola;asphere.type=SurfaceType::even_asphere;asphere.even_coefficients[0]=0.001;
    SequentialOpticalSystem aspheric;aspheric.add_surface(asphere);
    const auto ahit=aspheric.trace({{1.0,0.0,-1.0},{0.0,0.0,1.0},550.0});
    check(ahit.valid&&std::abs(ahit.ray.origin.z-0.051)<1.0e-12,"even-asphere intersection");
    SequentialSurface sphere{SurfaceType::sphere,0.0,10.0,2.0,1.5};
    SequentialSurface conic=sphere;conic.type=SurfaceType::conic;
    SequentialOpticalSystem a,b;a.add_surface(sphere);b.add_surface(conic);
    const Ray ray{{0.7,0.2,-2.0},{0.02,0.01,1.0},550.0};
    const auto ar=a.trace(ray),br=b.trace(ray);
    check(ar.valid&&br.valid&&std::abs(ar.ray.origin.z-br.ray.origin.z)<1.0e-11
        &&std::abs(ar.ray.direction.x-br.ray.direction.x)<1.0e-12,"zero-conic matches independent spherical intersection/refraction");
}
}

int main(){
    try{
        heat_and_coupling();fdtd_absorption_and_port();fdtd_ade_response();chemistry();fvm_workspaces();fvm_corrected_operators();transport_cache();
        trt<cfd::lbm::D2Q9InPlaceDescriptor>();trt<cfd::lbm::D3Q19Descriptor>();trt<cfd::lbm::D3Q27Descriptor>();optics();
        std::cout<<"All continuation regressions passed\n";return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
