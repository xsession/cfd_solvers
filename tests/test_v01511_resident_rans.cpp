#include "cfd/fvm/resident_rans_sycl.hpp"
#include "cfd/fvm/resident_sycl.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
void require(bool condition,std::string_view message){if(!condition)throw std::runtime_error(std::string(message));}

void test_cpu_reference_models(){
    const auto mesh=cfd::fvm::make_cartesian_hexa_mesh(8U,2U,1U,1.0,0.25,0.1);
    cfd::fvm::KEpsilonConfig c{};c.transport.dt=1.0e-3;
    cfd::fvm::KEpsilonTransport ke(mesh,c);ke.initialize(0.1,0.01);ke.set_strain_rate(0.0);
    const auto r=ke.step();
    require(r.first.converged&&r.second.converged,"CPU k-epsilon reference remains convergent");
    for(double v:ke.k())require(v>=c.minimum_k&&std::isfinite(v),"CPU k positivity");
    for(double v:ke.epsilon())require(v>=c.minimum_epsilon&&std::isfinite(v),"CPU epsilon positivity");
}

#if defined(CFD_HAS_SYCL)
void set_zero_gradient(cfd::fvm::ResidentSpalartAllmarasSycl& sa,const cfd::fvm::PolyMesh& mesh){for(const auto& p:mesh.patches())sa.set_boundary(p.name,cfd::fvm::TurbulenceScalarBoundaryType::zeroGradient);}

void test_shared_resident_rans_and_sources(){
    const auto mesh=cfd::fvm::make_cartesian_hexa_mesh(6U,3U,2U,1.0,0.5,0.4);
    cfd::fvm::ResidentIncompressibleConfig flow_cfg{};flow_cfg.dt=1.0e-3;flow_cfg.include_convection=false;
    cfd::fvm::ResidentIncompressibleSycl flow(mesh,flow_cfg,0U);
    flow.initialize_uniform({0.02,0.0,0.0},0.0);
    cfd::fvm::ResidentFvmFieldRegistrySycl fields(flow.mesh());

    cfd::fvm::SpalartAllmarasConfig sa_cfg{};sa_cfg.transport.dt=flow_cfg.dt;
    cfd::fvm::ResidentSpalartAllmarasSycl sa(mesh,fields,sa_cfg,"test.sa");
    set_zero_gradient(sa,mesh);sa.set_wall_distance(0.05);sa.initialize_uniform(2.0e-5);
    const auto sar=sa.step(flow.face_flux_device(),flow.velocity_device());require(sar.converged,"resident SA converges");

    cfd::fvm::KEpsilonConfig ke_cfg{};ke_cfg.transport.dt=flow_cfg.dt;
    cfd::fvm::ResidentKEpsilonSycl ke(mesh,fields,ke_cfg,"test.ke");
    for(const auto& p:mesh.patches()){ke.set_k_boundary(p.name,cfd::fvm::TurbulenceScalarBoundaryType::zeroGradient);ke.set_epsilon_boundary(p.name,cfd::fvm::TurbulenceScalarBoundaryType::zeroGradient);}
    ke.initialize_uniform(0.05,0.005);const auto ker=ke.step(flow.face_flux_device(),flow.velocity_device());require(ker.first.converged&&ker.second.converged,"resident k-epsilon converges");

    cfd::fvm::KOmegaSSTConfig sst_cfg{};sst_cfg.transport.dt=flow_cfg.dt;sst_cfg.des_enabled=true;
    cfd::fvm::ResidentKOmegaSSTSycl sst(mesh,fields,sst_cfg,"test.sst");
    for(const auto& p:mesh.patches()){sst.set_k_boundary(p.name,cfd::fvm::TurbulenceScalarBoundaryType::zeroGradient);sst.set_omega_boundary(p.name,cfd::fvm::TurbulenceScalarBoundaryType::zeroGradient);}
    sst.set_wall_distance(0.05);sst.initialize_uniform(0.05,5.0);const auto sr=sst.step(flow.face_flux_device(),flow.velocity_device());require(sr.first.converged&&sr.second.converged,"resident SST converges");

    cfd::fvm::ResidentScalarEquationSycl temperature(mesh,fields,"test.temperature",flow_cfg.dt,200U,1.0e-10);
    cfd::fvm::ResidentScalarEquationSycl fuel(mesh,fields,"test.fuel",flow_cfg.dt,200U,1.0e-10);
    for(const auto& p:mesh.patches()){temperature.set_boundary(p.name,cfd::fvm::TurbulenceScalarBoundaryType::zeroGradient);fuel.set_boundary(p.name,cfd::fvm::TurbulenceScalarBoundaryType::zeroGradient);}
    temperature.initialize_uniform(600.0);fuel.initialize_uniform(0.2);
    flow.mesh().fill(temperature.diffusivity_device(),mesh.cell_count(),1.0e-4);flow.mesh().fill(fuel.diffusivity_device(),mesh.cell_count(),1.0e-5);
    flow.mesh().fill(temperature.sink_device(),mesh.cell_count(),0.0);flow.mesh().fill(fuel.sink_device(),mesh.cell_count(),0.0);
    cfd::fvm::ResidentThermoSpeciesSourceSycl chemistry(fields,"test.chem");chemistry.arrhenius_one_step(temperature.field_device(),fuel.field_device(),2.0,1500.0,10.0);chemistry.apply_sources(temperature,fuel);
    temperature.reset_transfer_stats();fuel.reset_transfer_stats();
    const auto tr=temperature.step(flow.face_flux_device(),1.0,5000.0);const auto yr=fuel.step(flow.face_flux_device(),0.0,1.0);
    require(tr.converged&&yr.converged,"resident thermal/species equations converge");
    require(temperature.hot_loop_host_transfer_bytes()==0U&&fuel.hot_loop_host_transfer_bytes()==0U,"resident coupled scalar hot loop has zero host field traffic");
    require(fields.resident_bytes()>0U,"shared resident registry owns fields");

    std::vector<double> t(mesh.cell_count()),y(mesh.cell_count()),nut(mesh.cell_count());temperature.download(t);fuel.download(y);flow.mesh().download_cell_scalar(sst.eddy_viscosity_device(),nut);
    for(double v:t)require(v>=600.0&&std::isfinite(v),"reacting heat source raises/preserves temperature");
    for(double v:y)require(v>=0.0&&v<=0.2&&std::isfinite(v),"reacting species source consumes bounded fuel");
    for(double v:nut)require(v>=0.0&&std::isfinite(v),"resident SST eddy viscosity finite/nonnegative");
}
#endif
}

int main(){try{test_cpu_reference_models();
#if defined(CFD_HAS_SYCL)
    test_shared_resident_rans_and_sources();
#endif
    std::cout<<"v0.15.11 resident RANS/source tests passed\n";return EXIT_SUCCESS;}catch(const std::exception& e){std::cerr<<"v0.15.11 test failure: "<<e.what()<<'\n';return EXIT_FAILURE;}}
