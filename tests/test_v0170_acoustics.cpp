#include "cfd/acoustics/kspace2d.hpp"
#include "cfd/acoustics/postprocess.hpp"
#include "cfd/spectral/fft.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool v,const char* msg){if(!v)throw std::runtime_error(msg);}

double max_abs(std::span<const double> v){double m=0;for(double x:v)m=std::max(m,std::abs(x));return m;}

double mode_amplitude(std::span<const double> v,std::size_t nx,std::size_t ny,std::size_t mx,std::size_t my){
    double cs=0.0,sn=0.0;
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i){
        const double phase=2.0*std::numbers::pi*(static_cast<double>(mx*i)/static_cast<double>(nx)+static_cast<double>(my*j)/static_cast<double>(ny));
        const double x=v[j*nx+i];cs+=x*std::cos(phase);sn+=x*std::sin(phase);
    }
    return 2.0*std::hypot(cs,sn)/static_cast<double>(nx*ny);
}

void test_fft_derivative(){
    constexpr std::size_t n=64;constexpr double length=2.0;const double dx=length/static_cast<double>(n);constexpr std::size_t mode=5;
    std::vector<double> f(n),exact(n);
    const double k=2.0*std::numbers::pi*static_cast<double>(mode)/length;
    for(std::size_t i=0;i<n;++i){const double x=static_cast<double>(i)*dx;f[i]=std::sin(k*x);exact[i]=k*std::cos(k*x);}
    const auto d=cfd::spectral::derivative_periodic_1d(f,dx);
    double err=0;for(std::size_t i=0;i<n;++i)err=std::max(err,std::abs(d[i]-exact[i]));
    require(err<2e-12,"radix-2 spectral derivative accuracy");
}

void test_fractional_laplacian_and_power_law(){
    constexpr std::size_t nx=32,ny=32;constexpr double dx=0.02,dy=0.025;constexpr std::size_t mx=2,my=3;
    std::vector<double> f(nx*ny);
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i){const double phase=2.0*std::numbers::pi*(static_cast<double>(mx*i)/nx+static_cast<double>(my*j)/ny);f[j*nx+i]=std::sin(phase);}
    constexpr double order=0.75;
    const auto fl=cfd::acoustics::fractional_laplacian_periodic_2d(f,nx,ny,dx,dy,order);
    const double kx=2.0*std::numbers::pi*static_cast<double>(mx)/(nx*dx),ky=2.0*std::numbers::pi*static_cast<double>(my)/(ny*dy);
    const double eigen=std::pow(kx*kx+ky*ky,order);
    double err=0;for(std::size_t q=0;q<f.size();++q)err=std::max(err,std::abs(fl[q]-eigen*f[q]));
    require(err<2e-9*eigen,"2-D fractional Laplacian eigenmode");

    std::vector<double> mixed(nx*ny);
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i){const double p1=2.0*std::numbers::pi*static_cast<double>(i)/nx;const double p5=10.0*std::numbers::pi*static_cast<double>(i)/nx;mixed[j*nx+i]=std::sin(p1)+std::sin(p5);}
    cfd::acoustics::PowerLawAbsorption2D a{2.0,1.5,2.0e5};
    const auto filtered=cfd::acoustics::power_law_attenuate_periodic_2d(mixed,nx,ny,dx,dy,1500.0,2e-5,a);
    const double low_ratio=mode_amplitude(filtered,nx,ny,1,0)/mode_amplitude(mixed,nx,ny,1,0);
    const double high_ratio=mode_amplitude(filtered,nx,ny,5,0)/mode_amplitude(mixed,nx,ny,5,0);
    require(low_ratio<1.0&&high_ratio<low_ratio,"power-law attenuation increases with frequency");
}

void test_kspace_propagation_and_sensors(){
    cfd::acoustics::KSpace2DConfig cfg;cfg.nx=64;cfg.ny=32;cfg.dx=2.5e-4;cfg.dy=2.5e-4;cfg.cfl=0.15;cfg.pml_cells=0;
    cfd::acoustics::KSpaceAcoustic2D solver(cfg);solver.set_uniform_medium(1500.0,1000.0);
    std::vector<double> p(cfg.nx*cfg.ny);
    constexpr std::size_t mode=3;
    for(std::size_t j=0;j<cfg.ny;++j)for(std::size_t i=0;i<cfg.nx;++i)p[j*cfg.nx+i]=1000.0*std::sin(2.0*std::numbers::pi*static_cast<double>(mode*i)/static_cast<double>(cfg.nx));
    solver.initialize_pressure(p);
    const double e0=solver.total_energy();
    const auto sensor=solver.add_sensor(cfg.nx/4U+(cfg.ny/2U)*cfg.nx);
    solver.run(120);
    const double e1=solver.total_energy();
    require(std::abs(e1-e0)/e0<3e-2,"lossless k-space acoustic bounded energy error");
    require(solver.sensors()[sensor].pressure.size()==120,"sensor records every propagated step");
    require(std::isfinite(max_abs(solver.pressure())),"finite acoustic pressure after propagation");
}

void test_heterogeneous_nonlinear_pml_photoacoustic(){
    cfd::acoustics::KSpace2DConfig cfg;cfg.nx=64;cfg.ny=64;cfg.dx=2e-4;cfg.dy=2e-4;cfg.cfl=0.1;cfg.pml_cells=10;
    cfd::acoustics::KSpaceAcoustic2D solver(cfg);
    std::vector<double> c(cfg.nx*cfg.ny,1500.0),rho(cfg.nx*cfg.ny,1000.0);
    for(std::size_t j=20;j<44;++j)for(std::size_t i=32;i<cfg.nx;++i){c[j*cfg.nx+i]=1600.0;rho[j*cfg.nx+i]=1080.0;}
    solver.set_medium(c,rho);solver.set_uniform_nonlinearity(6.0);
    solver.initialize_gaussian_pressure(0.35,0.5,7e-4,2e5);
    const auto s=solver.add_sensor(48+(cfg.ny/2U)*cfg.nx);
    solver.run(240);
    require(!solver.sensors()[s].pressure.empty()&&max_abs(solver.sensors()[s].pressure)>1.0,"heterogeneous photoacoustic pulse reaches sensor");
    require(std::isfinite(solver.total_energy())&&solver.total_energy()>0.0,"heterogeneous nonlinear run remains finite");

    cfd::acoustics::KSpaceAcoustic2D forced(cfg);forced.set_uniform_medium(1500,1000);forced.add_pressure_source((cfg.nx/2U)+(cfg.ny/2U)*cfg.nx,[](double t){return 2e3*std::sin(2.0*std::numbers::pi*2e5*t);});
    const auto fs=forced.add_sensor((cfg.nx/2U+5U)+(cfg.ny/2U)*cfg.nx);forced.run(80);
    require(max_abs(forced.sensors()[fs].pressure)>1e-4,"time-varying pressure source propagates to sensor");
}

void test_postprocess_and_arrays(){
    const double intensity=cfd::acoustics::plane_wave_intensity_w_m2(1.0e5,1000.0,1500.0);
    require(std::abs(intensity-6666.666666666667)<1e-9,"plane-wave acoustic intensity");
    require(std::abs(cfd::acoustics::radiation_pressure_pa(intensity,1500.0,1.0)-8.88888888888889)<1e-10,"reflected acoustic radiation pressure");
    require(std::abs(cfd::acoustics::acoustic_heating_w_m3(intensity,2.0)-26666.666666666668)<1e-8,"acoustic absorption heating");

    std::vector<double> a(64,0.0),b(64,0.0);a[20]=1.0;b[22]=1.0;
    std::vector<cfd::acoustics::BeamformSample> samples{{-0.001,0.0,a},{0.001,0.0,b}};
    const auto beam=cfd::acoustics::delay_and_sum_beamform(samples,0.001,0.003,1.0e-6,1500.0,48);
    require(max_abs(beam)>0.0,"delay-and-sum beamformer returns focused signal");

    cfd::acoustics::KSpace2DConfig cfg;cfg.nx=32;cfg.ny=16;cfg.dx=5e-4;cfg.dy=5e-4;cfg.pml_cells=2;cfg.cfl=0.1;
    cfd::acoustics::KSpaceAcoustic2D solver(cfg);
    solver.add_plane_source_x(3,[](double t){return 100.0*std::sin(2.0*std::numbers::pi*5e4*t);});
    const std::vector<std::size_t> cells{8+8*cfg.nx,8+9*cfg.nx};const std::vector<double> delay{0.0,1e-6};
    solver.add_delayed_source_array(cells,delay,[](double t){return 50.0*std::sin(2.0*std::numbers::pi*5e4*t);});
    const auto ids=solver.add_sensor_array(cells);solver.run(20);
    require(ids.size()==2&&solver.sensors()[ids[0]].pressure.size()==20,"source and sensor array helpers");
}
}

int main(){
    try{test_fft_derivative();test_fractional_laplacian_and_power_law();test_kspace_propagation_and_sensors();test_heterogeneous_nonlinear_pml_photoacoustic();test_postprocess_and_arrays();}
    catch(const std::exception& e){std::cerr<<"v0.17.0 acoustics regression failed: "<<e.what()<<'\n';return 1;}
    std::cout<<"v0.17.0 acoustics regression passed\n";return 0;
}
