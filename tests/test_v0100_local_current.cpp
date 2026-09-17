#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

std::size_t idx3(std::size_t ix,std::size_t iy,std::size_t iz,std::size_t nx,std::size_t ny){return (iz*ny+iy)*nx+ix;}

double finite_volume_residual(const cfd::particle::ChargeConservingCurrent3D& current,
                              std::size_t nx,std::size_t ny,std::size_t nz,
                              double lx,double ly,double lz,double dt){
    const double dx=lx/static_cast<double>(nx),dy=ly/static_cast<double>(ny),dz=lz/static_cast<double>(nz);
    double residual=0.0;
    for(std::size_t iz=0;iz<nz;++iz){
        const std::size_t iz_back=(iz+nz-1U)%nz;
        for(std::size_t iy=0;iy<ny;++iy){
            const std::size_t iy_back=(iy+ny-1U)%ny;
            for(std::size_t ix=0;ix<nx;++ix){
                const std::size_t ix_back=(ix+nx-1U)%nx;
                const std::size_t idx=idx3(ix,iy,iz,nx,ny);
                const double delta=(current.charge_density_new_c_per_m3[idx]-current.charge_density_old_c_per_m3[idx])/dt;
                const double div=(current.current_x_a_per_m2[idx]-current.current_x_a_per_m2[idx3(ix_back,iy,iz,nx,ny)])/dx
                    +(current.current_y_a_per_m2[idx]-current.current_y_a_per_m2[idx3(ix,iy_back,iz,nx,ny)])/dy
                    +(current.current_z_a_per_m2[idx]-current.current_z_a_per_m2[idx3(ix,iy,iz_back,nx,ny)])/dz;
                residual=std::max(residual,std::abs(delta+div));
            }
        }
    }
    return residual;
}

std::vector<cfd::particle::PicParticle3D> make_particles(){
    using cfd::particle::PicParticle3D;
    std::vector<PicParticle3D> particles;
    for(std::size_t i=0;i<8U;++i){
        PicParticle3D p;p.mass_kg=1.0e-6;p.charge_c=(i%2U==0U?1.0:-1.0)*1.0e-15;p.weight=1.0+0.05*static_cast<double>(i);
        p.position_m={0.07+0.091*static_cast<double>(i),0.06+0.071*static_cast<double>(i%5U),0.05+0.053*static_cast<double>(i%4U)};
        p.velocity_m_per_s={4.0e3,-2.5e3,1.7e3};particles.push_back(p);
    }
    return particles;
}

void local_current_satisfies_finite_volume_continuity(){
    using namespace cfd::particle;
    constexpr std::size_t nx=7U,ny=6U,nz=5U;constexpr double lx=1.0,ly=0.8,lz=0.6,dt=2.0e-10;
    auto old_particles=make_particles();auto new_particles=old_particles;
    for(std::size_t i=0;i<new_particles.size();++i){
        new_particles[i].position_m.x+=0.011+1.0e-3*static_cast<double>(i);
        new_particles[i].position_m.y-=0.007;
        new_particles[i].position_m.z+=0.004;
    }
    const auto local=deposit_charge_conserving_current_3d_local(old_particles,new_particles,nx,ny,nz,lx,ly,lz,dt);
    const auto spectral=deposit_charge_conserving_current_3d(old_particles,new_particles,nx,ny,nz,lx,ly,lz,dt);
    double max_current=0.0;
    for(double value:local.current_x_a_per_m2)max_current=std::max(max_current,std::abs(value));
    for(double value:local.current_y_a_per_m2)max_current=std::max(max_current,std::abs(value));
    for(double value:local.current_z_a_per_m2)max_current=std::max(max_current,std::abs(value));
    require(max_current>0.0,"local finite-volume current is non-zero for moving particles");
    require(local.continuity_linf_residual<1.0e-8,"local finite-volume current reports small continuity residual");
    require(finite_volume_residual(local,nx,ny,nz,lx,ly,lz,dt)<1.0e-8,"local current satisfies periodic finite-volume continuity");
    require(spectral.continuity_linf_residual<1.0e-8,"spectral reference current still satisfies continuity");
    require(local.charge_density_new_c_per_m3.size()==spectral.charge_density_new_c_per_m3.size(),"local/spectral charge arrays have matching size");
}

void staggered_pic_uses_local_current_mode(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic3DConfig config;config.nx=7U;config.ny=6U;config.nz=5U;config.length_x_m=1.0;config.length_y_m=0.8;config.length_z_m=0.6;config.current_deposition=CurrentDeposition3DMode::local_finite_volume;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);const double dz=config.length_z_m/static_cast<double>(config.nz);
    config.dt_s=0.03/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
    StaggeredElectromagneticPic3D solver(config);solver.set_particles(make_particles());solver.step(1U);
    const auto diag=solver.diagnostics();
    double max_current=0.0;
    for(double value:solver.current_x())max_current=std::max(max_current,std::abs(value));
    for(double value:solver.current_y())max_current=std::max(max_current,std::abs(value));
    for(double value:solver.current_z())max_current=std::max(max_current,std::abs(value));
    require(diag.particle_count==make_particles().size(),"local-current staggered PIC preserves periodic particles");
    require(max_current>0.0,"local-current staggered PIC couples finite-volume current");
    require(diag.charge_continuity_linf_residual<1.0e-6,"local-current staggered PIC preserves continuity");
}

void sponge_polynomial_profile_is_configurable(){
    using namespace cfd::particle;
    constexpr std::size_t nx=7U,ny=7U,nz=7U;const std::size_t n=nx*ny*nz;
    std::vector<double> ex(n,1.0),ey(n,1.0),ez(n,1.0),bx(n,1.0),by(n,1.0),bz(n,1.0);
    GridBoundary3DConfig linear;linear.mode=GridBoundaryMode3D::absorbing_sponge;linear.sponge_cells=3U;linear.sponge_strength=2.0;linear.sponge_polynomial_order=1.0;
    auto ex_linear=ex,ey_linear=ey,ez_linear=ez,bx_linear=bx,by_linear=by,bz_linear=bz;
    apply_electromagnetic_field_boundary_3d(ex_linear,ey_linear,ez_linear,bx_linear,by_linear,bz_linear,nx,ny,nz,linear);
    GridBoundary3DConfig quartic=linear;quartic.sponge_polynomial_order=4.0;
    auto ex_quartic=ex,ey_quartic=ey,ez_quartic=ez,bx_quartic=bx,by_quartic=by,bz_quartic=bz;
    apply_electromagnetic_field_boundary_3d(ex_quartic,ey_quartic,ez_quartic,bx_quartic,by_quartic,bz_quartic,nx,ny,nz,quartic);
    const std::size_t near_edge=idx3(1U,3U,3U,nx,ny);
    require(ex_linear[near_edge]<ex_quartic[near_edge],"higher-order sponge damps less in the inner part of the layer");
    require(ex_linear.front()==ex_quartic.front(),"all sponge profiles share the same outer-wall damping");
}
}

int main(){
    try{
        local_current_satisfies_finite_volume_continuity();
        staggered_pic_uses_local_current_mode();
        sponge_polynomial_profile_is_configurable();
        std::cout<<"v0.10.0 local 3-D current/PIC boundary tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.10.0 local current test failure: "<<error.what()<<'\n';
        return 1;
    }
}
