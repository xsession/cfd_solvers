#include "common/benchmark.hpp"
#include "cfd/solvers/fem/dynamic_bar1d.hpp"
#include <cmath>
#include <numbers>
int main(int argc,char**argv){using namespace cfd::examples;const auto o=parse_options(argc,argv);cfd::fem::DynamicBar1DConfig cfg;cfg.elements=o.quick?32U:256U*o.scale;cfg.length=1.0;cfg.area=1e-4;cfg.young_modulus=70e9;cfg.density=2700;cfg.dt=2e-6;const std::size_t steps=o.quick?20U:1200U;Timer setup;cfd::fem::DynamicBar1D bar(cfg);bar.initialize([&](double x){return 1e-5*std::sin(0.5*std::numbers::pi*x/cfg.length);});const double setup_ms=setup.milliseconds();Timer sim;bar.run(steps);const double sim_ms=sim.milliseconds();const double checksum=bar.kinetic_energy()+bar.strain_energy()+bar.displacement().back();emit({"phase04_fem","dynamic_cantilever_bar","cpu","dof_steps",cfg.elements+1U,steps,setup_ms,sim_ms,double(cfg.elements+1U)*double(steps),checksum});return 0;}
