#include "common/benchmark.hpp"
#include "cfd/core/geometric_multigrid.hpp"
#include <cmath>
#include <numbers>
#include <vector>
int main(int argc,char**argv){using namespace cfd::examples;const auto o=parse_options(argc,argv);const std::size_t n=o.quick?31U:127U*o.scale;Timer setup;cfd::core::PoissonMultigrid2D mg(n);std::vector<double>f(n*n),u(n*n);for(std::size_t j=0;j<n;++j)for(std::size_t i=0;i<n;++i){const double x=double(i+1)/double(n+1),y=double(j+1)/double(n+1);f[j*n+i]=2.0*std::numbers::pi*std::numbers::pi*std::sin(std::numbers::pi*x)*std::sin(std::numbers::pi*y);}const double setup_ms=setup.milliseconds();Timer sim;const auto result=mg.solve(f,u,o.quick?12U:60U,1e-8);const double sim_ms=sim.milliseconds();double checksum=0;for(double v:u)checksum+=v;emit({"phase01_hpc_runtime","poisson_multigrid","cpu","unknowns_cycles",n*n,result.cycles,setup_ms,sim_ms,double(n*n)*double(result.cycles),checksum});return result.converged?0:1;}
