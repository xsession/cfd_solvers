#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
constexpr double eps0=8.8541878128e-12;
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

double rms_error(const std::vector<double>& actual,const std::vector<double>& expected){
    double numerator=0.0,denominator=0.0;
    for(std::size_t i=0;i<actual.size();++i){const double delta=actual[i]-expected[i];numerator+=delta*delta;denominator+=expected[i]*expected[i];}
    return std::sqrt(numerator/std::max(denominator,1.0e-300));
}

void three_dimensional_cic_conserves_charge(){
    using namespace cfd::particle;
    constexpr std::size_t nx=7U,ny=6U,nz=5U;constexpr double lx=1.0,ly=0.8,lz=0.6;
    std::vector<PicParticle3D> particles(3U);
    particles[0].mass_kg=1.0;particles[0].charge_c=2.0e-15;particles[0].weight=1.0;particles[0].position_m={0.17,0.23,0.31};
    particles[1].mass_kg=1.0;particles[1].charge_c=-1.0e-15;particles[1].weight=2.0;particles[1].position_m={1.11,-0.04,0.58};
    particles[2].mass_kg=1.0;particles[2].charge_c=0.5e-15;particles[2].weight=3.0;particles[2].position_m={0.49,0.79,0.12};
    const auto rho=deposit_cic_charge_density_3d(particles,nx,ny,nz,lx,ly,lz);
    const double cell=(lx/static_cast<double>(nx))*(ly/static_cast<double>(ny))*(lz/static_cast<double>(nz));
    const double grid_charge=std::accumulate(rho.begin(),rho.end(),0.0)*cell;
    double particle_charge=0.0;for(const auto& particle:particles)particle_charge+=particle.charge_c*particle.weight;
    require(std::abs(grid_charge-particle_charge)<1.0e-29,"3-D CIC deposition conserves total macro-charge");
}

void three_dimensional_periodic_poisson(){
    using namespace cfd::particle;
    constexpr std::size_t nx=8U,ny=7U,nz=6U;constexpr double lx=1.0,ly=0.7,lz=0.5,epsr=1.7;
    const double kx=2.0*std::numbers::pi/lx;const double ky=2.0*std::numbers::pi/ly;const double kz=2.0*std::numbers::pi/lz;
    std::vector<double> rho(nx*ny*nz),ex(nx*ny*nz),ey(nx*ny*nz),ez(nx*ny*nz);
    for(std::size_t iz=0;iz<nz;++iz){
        const double z=lz*static_cast<double>(iz)/static_cast<double>(nz);
        for(std::size_t iy=0;iy<ny;++iy){
            const double y=ly*static_cast<double>(iy)/static_cast<double>(ny);
            for(std::size_t ix=0;ix<nx;++ix){
                const double x=lx*static_cast<double>(ix)/static_cast<double>(nx);const double phi=std::sin(kx*x)*std::cos(ky*y)*std::cos(kz*z);const std::size_t idx=(iz*ny+iy)*nx+ix;
                rho[idx]=eps0*epsr*(kx*kx+ky*ky+kz*kz)*phi;
                ex[idx]=-kx*std::cos(kx*x)*std::cos(ky*y)*std::cos(kz*z);
                ey[idx]= ky*std::sin(kx*x)*std::sin(ky*y)*std::cos(kz*z);
                ez[idx]= kz*std::sin(kx*x)*std::cos(ky*y)*std::sin(kz*z);
            }
        }
    }
    const auto field=periodic_electric_field_from_charge_density_3d(rho,nx,ny,nz,lx,ly,lz,epsr,true);
    require(rms_error(field.electric_x_v_per_m,ex)<1.0e-12,"3-D periodic Poisson Ex matches manufactured sinusoid");
    require(rms_error(field.electric_y_v_per_m,ey)<1.0e-12,"3-D periodic Poisson Ey matches manufactured sinusoid");
    require(rms_error(field.electric_z_v_per_m,ez)<1.0e-12,"3-D periodic Poisson Ez matches manufactured sinusoid");
}

void charge_conserving_current_3d(){
    using namespace cfd::particle;
    std::vector<PicParticle3D> before(6U),after(6U);
    for(std::size_t i=0;i<before.size();++i){
        before[i].mass_kg=1.0;before[i].charge_c=(i%2U==0U?1.0:-0.8)*1.0e-12;before[i].weight=1.0+0.1*static_cast<double>(i);
        before[i].position_m={0.08+0.13*static_cast<double>(i),0.11+0.09*static_cast<double>(i),0.07+0.05*static_cast<double>(i)};
        after[i]=before[i];after[i].position_m.x+=0.013+0.001*static_cast<double>(i);after[i].position_m.y-=0.009+0.002*static_cast<double>(i);after[i].position_m.z+=0.006+0.001*static_cast<double>(i);
    }
    const auto current=deposit_charge_conserving_current_3d(before,after,8U,7U,6U,1.0,0.8,0.6,2.0e-9);
    require(current.current_x_a_per_m2.size()==8U*7U*6U&&current.current_y_a_per_m2.size()==8U*7U*6U&&current.current_z_a_per_m2.size()==8U*7U*6U,"3-D current deposition returns grid-sized vectors");
    require(current.continuity_linf_residual<5.0e-9,"3-D spectral current reconstruction satisfies continuity");
}

void electrostatic_pic3d_and_gather(){
    using namespace cfd::particle;
    ElectrostaticPic3DConfig config;config.nx=8U;config.ny=7U;config.nz=6U;config.length_x_m=1.0;config.length_y_m=0.75;config.length_z_m=0.6;config.dt_s=2.0e-6;
    ElectrostaticPic3D solver(config);std::vector<PicParticle3D> particles(2U);
    particles[0].mass_kg=1.0;particles[0].charge_c=1.0e-15;particles[0].position_m={0.23,0.27,0.13};
    particles[1].mass_kg=1.0;particles[1].charge_c=-1.0e-15;particles[1].position_m={0.61,0.52,0.41};
    solver.set_particles(particles);solver.deposit_and_solve();double peak_field=0.0;for(double value:solver.electric_x())peak_field=std::max(peak_field,std::abs(value));for(double value:solver.electric_y())peak_field=std::max(peak_field,std::abs(value));for(double value:solver.electric_z())peak_field=std::max(peak_field,std::abs(value));
    require(peak_field>0.0,"3-D electrostatic PIC produces nonzero local field for separated neutral pairs");
    const auto field=solver.gather_electric_field({0.33,0.29,0.21});require(std::isfinite(field.x)&&std::isfinite(field.y)&&std::isfinite(field.z),"3-D field gather returns finite trilinear values");
    solver.step(3U);require(std::abs(solver.time_s()-6.0e-6)<1.0e-15,"3-D electrostatic PIC advances time");
}
}

int main(){
    try{
        three_dimensional_cic_conserves_charge();
        three_dimensional_periodic_poisson();
        charge_conserving_current_3d();
        electrostatic_pic3d_and_gather();
        std::cout<<"v0.9.7 3-D PIC foundation tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.9.7 PIC3D test failure: "<<error.what()<<'\n';
        return 1;
    }
}
