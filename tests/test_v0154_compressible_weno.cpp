#include "cfd/solvers/fvm/compressible1d.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}

std::vector<cfd::fvm::EulerPrimitive1D> entropy_wave(std::size_t cells,double shift){
    const double dx=1.0/static_cast<double>(cells);
    const double average_factor=std::sin(std::numbers::pi*dx)/(std::numbers::pi*dx);
    std::vector<cfd::fvm::EulerPrimitive1D> profile(cells);
    for(std::size_t i=0;i<cells;++i){
        const double x=(static_cast<double>(i)+0.5)*dx-shift;
        profile[i]={1.0+0.2*average_factor*std::sin(2.0*std::numbers::pi*x),1.0,1.0};
    }
    return profile;
}

double entropy_error(std::size_t cells,cfd::fvm::CompressibleReconstruction reconstruction,cfd::fvm::CompressibleTimeIntegrator integrator){
    cfd::fvm::Compressible1DConfig cfg;
    cfg.cells=cells;cfg.length=1.0;cfg.cfl=0.25;cfg.periodic=true;cfg.reconstruction=reconstruction;cfg.time_integrator=integrator;
    cfd::fvm::CompressibleEuler1D solver(cfg);
    solver.initialize_profile(entropy_wave(cells,0.0));
    const double m0=solver.mass(),q0=solver.momentum(),e0=solver.total_energy();
    solver.run(0.2);
    const auto exact=entropy_wave(cells,0.2);
    double error=0.0;
    for(std::size_t i=0;i<cells;++i)error+=std::abs(solver.primitive(i).density-exact[i].density);
    error/=static_cast<double>(cells);
    require(std::abs(solver.mass()-m0)<2.0e-12,"periodic WENO must conserve mass");
    require(std::abs(solver.momentum()-q0)<2.0e-12,"periodic WENO must conserve momentum");
    require(std::abs(solver.total_energy()-e0)<4.0e-12,"periodic WENO must conserve total energy");
    return error;
}

void test_characteristic_weno_smooth_accuracy(){
    const double e40=entropy_error(40U,cfd::fvm::CompressibleReconstruction::characteristic_weno5,cfd::fvm::CompressibleTimeIntegrator::ssprk3);
    const double e80=entropy_error(80U,cfd::fvm::CompressibleReconstruction::characteristic_weno5,cfd::fvm::CompressibleTimeIntegrator::ssprk3);
    const double e160=entropy_error(160U,cfd::fvm::CompressibleReconstruction::characteristic_weno5,cfd::fvm::CompressibleTimeIntegrator::ssprk3);
    const double rate40_80=std::log(e40/e80)/std::log(2.0);
    const double rate80_160=std::log(e80/e160)/std::log(2.0);
    std::cout<<"entropy-wave WENO errors "<<e40<<" "<<e80<<" "<<e160<<" rates "<<rate40_80<<" "<<rate80_160<<'\n';
    require(rate40_80>2.4&&rate80_160>2.4,"WENO5 plus SSPRK3 should show at least third-order-like smooth convergence");
    const double muscl=entropy_error(80U,cfd::fvm::CompressibleReconstruction::muscl_minmod,cfd::fvm::CompressibleTimeIntegrator::ssprk3);
    require(e80<muscl*0.25,"characteristic WENO should materially outperform MUSCL on a smooth entropy wave");
}

void test_uniform_preservation(){
    cfd::fvm::Compressible1DConfig cfg;cfg.cells=48U;cfg.periodic=true;cfg.reconstruction=cfd::fvm::CompressibleReconstruction::characteristic_weno5;cfg.time_integrator=cfd::fvm::CompressibleTimeIntegrator::ssprk3;
    cfd::fvm::CompressibleEuler1D solver(cfg);solver.initialize_uniform({1.2,0.7,2.5});solver.run(0.1);
    for(std::size_t i=0;i<cfg.cells;++i){const auto q=solver.primitive(i);require(std::abs(q.density-1.2)<2e-13&&std::abs(q.velocity-0.7)<2e-13&&std::abs(q.pressure-2.5)<5e-13,"WENO must exactly preserve a uniform state");}
}

void test_sod_positivity_and_shock_sensor(){
    cfd::fvm::Compressible1DConfig cfg;cfg.cells=240U;cfg.length=1.0;cfg.cfl=0.2;cfg.periodic=false;cfg.reconstruction=cfd::fvm::CompressibleReconstruction::characteristic_weno5;cfg.time_integrator=cfd::fvm::CompressibleTimeIntegrator::ssprk3;
    cfd::fvm::CompressibleEuler1D solver(cfg);solver.initialize({1.0,0.0,1.0},{0.125,0.0,0.1},0.5);
    const auto initial_sensor=solver.pressure_jump_sensor();
    require(*std::max_element(initial_sensor.begin(),initial_sensor.end())>0.7,"pressure jump sensor should detect the initial Sod discontinuity");
    solver.run(0.12);
    double minimum_density=1e300,minimum_pressure=1e300,maximum_sensor=0.0;
    const auto sensor=solver.pressure_jump_sensor();
    for(std::size_t i=0;i<cfg.cells;++i){const auto q=solver.primitive(i);minimum_density=std::min(minimum_density,q.density);minimum_pressure=std::min(minimum_pressure,q.pressure);maximum_sensor=std::max(maximum_sensor,sensor[i]);}
    require(minimum_density>0.0&&minimum_pressure>0.0,"characteristic WENO shock evolution must remain positive");
    require(maximum_sensor>0.01,"shock sensor should remain active on evolved discontinuities");
}
}

int main(){try{test_characteristic_weno_smooth_accuracy();test_uniform_preservation();test_sod_positivity_and_shock_sensor();std::cout<<"v0.15.4 characteristic WENO compressible tests passed\n";return 0;}catch(const std::exception& e){std::cerr<<"v0.15.4 characteristic WENO test failure: "<<e.what()<<'\n';return 1;}}
