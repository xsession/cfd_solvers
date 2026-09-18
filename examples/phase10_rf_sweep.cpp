#include "common/benchmark.hpp"
#include "cfd/rf/network.hpp"
#include <cmath>
#include <complex>
int main(int argc,char**argv){using namespace cfd::examples;const auto o=parse_options(argc,argv);const std::size_t points=o.quick?200U:50000U*o.scale;Timer setup;const double f0=1e8,f1=6e9;const double setup_ms=setup.milliseconds();Timer sim;double checksum=0.0;for(std::size_t i=0;i<points;++i){const double f=f0+(f1-f0)*double(i)/double(points-1U);const auto m=cfd::rf::microstrip_quasi_static(2e-3,1e-3,4.4,f);cfd::rf::TemTransmissionLine line{m.characteristic_impedance,m.phase_velocity,0.02,0.05};const auto s=line.s_parameters(f,50.0);checksum+=std::abs(s.a21)+m.guided_wavelength;}const double sim_ms=sim.milliseconds();emit({"phase10_rf","microstrip_frequency_sweep","cpu","frequency_points",points,1U,setup_ms,sim_ms,double(points),checksum});return 0;}
