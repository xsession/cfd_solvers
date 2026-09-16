#include "cfd/solvers/fdtd/maxwell3d.hpp"

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fdtd {
namespace {
constexpr double epsilon0=8.8541878128e-12;
constexpr double mu0=1.25663706212e-6;
}

Maxwell3D::Maxwell3D(Maxwell3DConfig c)
 :config_(c),epsilon_(epsilon0*c.epsilon_r),mu_(mu0*c.mu_r),
  ex_(c.nx*c.ny*c.nz,0.0),ey_(ex_.size(),0.0),ez_(ex_.size(),0.0),
  hx_(ex_.size(),0.0),hy_(ex_.size(),0.0),hz_(ex_.size(),0.0){
    if(c.nx<4U||c.ny<4U||c.nz<4U||!(c.dx>0.0)||!(c.dy>0.0)||!(c.dz>0.0)||
       !(c.courant>0.0&&c.courant<1.0)||!(c.epsilon_r>0.0)||!(c.mu_r>0.0))
        throw std::invalid_argument("invalid Maxwell3D configuration");
    const double wave=1.0/std::sqrt(epsilon_*mu_);
    dt_=c.courant/(wave*std::sqrt(1.0/(c.dx*c.dx)+1.0/(c.dy*c.dy)+1.0/(c.dz*c.dz)));
}
std::size_t Maxwell3D::idx(std::size_t i,std::size_t j,std::size_t k) const noexcept{
    return (k*config_.ny+j)*config_.nx+i;
}
void Maxwell3D::initialize_gaussian_ez(double amplitude,double wf){
    if(!(wf>0.0)) throw std::invalid_argument("Gaussian width must be positive");
    const double cx=0.5*static_cast<double>(config_.nx-1U),cy=0.5*static_cast<double>(config_.ny-1U),cz=0.5*static_cast<double>(config_.nz-1U);
    const double w=wf*static_cast<double>(std::min({config_.nx,config_.ny,config_.nz}));
    cfd::core::parallel_for(cell_count(),[&](std::size_t n){
        const std::size_t i=n%config_.nx,j=(n/config_.nx)%config_.ny,k=n/(config_.nx*config_.ny);
        const double r2=((static_cast<double>(i)-cx)*(static_cast<double>(i)-cx)+
                         (static_cast<double>(j)-cy)*(static_cast<double>(j)-cy)+
                         (static_cast<double>(k)-cz)*(static_cast<double>(k)-cz))/(w*w);
        ez_[n]=amplitude*std::exp(-0.5*r2);
    });
    enforce_boundary();
}
void Maxwell3D::add_soft_ez_source(std::size_t i,std::size_t j,std::size_t k,double value){
    if(i>=config_.nx||j>=config_.ny||k>=config_.nz||!std::isfinite(value)) throw std::out_of_range("invalid Maxwell3D soft source");
    ez_[idx(i,j,k)]+=value;
}
void Maxwell3D::step(std::size_t count){
    const double hdt=dt_/mu_,edt=dt_/epsilon_;
    const std::size_t nx=config_.nx,ny=config_.ny,nz=config_.nz;
    for(std::size_t s=0;s<count;++s){
        cfd::core::parallel_for(nx*(ny-1U)*(nz-1U),[&](std::size_t q){
            const std::size_t i=q%nx,j=(q/nx)%(ny-1U),k=q/(nx*(ny-1U)); const auto n=idx(i,j,k);
            hx_[n]-=hdt*((ez_[idx(i,j+1U,k)]-ez_[n])/config_.dy-(ey_[idx(i,j,k+1U)]-ey_[n])/config_.dz);
        });
        cfd::core::parallel_for((nx-1U)*ny*(nz-1U),[&](std::size_t q){
            const std::size_t i=q%(nx-1U),j=(q/(nx-1U))%ny,k=q/((nx-1U)*ny); const auto n=idx(i,j,k);
            hy_[n]-=hdt*((ex_[idx(i,j,k+1U)]-ex_[n])/config_.dz-(ez_[idx(i+1U,j,k)]-ez_[n])/config_.dx);
        });
        cfd::core::parallel_for((nx-1U)*(ny-1U)*nz,[&](std::size_t q){
            const std::size_t i=q%(nx-1U),j=(q/(nx-1U))%(ny-1U),k=q/((nx-1U)*(ny-1U)); const auto n=idx(i,j,k);
            hz_[n]-=hdt*((ey_[idx(i+1U,j,k)]-ey_[n])/config_.dx-(ex_[idx(i,j+1U,k)]-ex_[n])/config_.dy);
        });
        cfd::core::parallel_for((nx-2U)*(ny-2U)*(nz-2U),[&](std::size_t q){
            const std::size_t i=1U+q%(nx-2U),j=1U+(q/(nx-2U))%(ny-2U),k=1U+q/((nx-2U)*(ny-2U)); const auto n=idx(i,j,k);
            ex_[n]+=edt*((hz_[n]-hz_[idx(i,j-1U,k)])/config_.dy-(hy_[n]-hy_[idx(i,j,k-1U)])/config_.dz);
            ey_[n]+=edt*((hx_[n]-hx_[idx(i,j,k-1U)])/config_.dz-(hz_[n]-hz_[idx(i-1U,j,k)])/config_.dx);
            ez_[n]+=edt*((hy_[n]-hy_[idx(i-1U,j,k)])/config_.dx-(hx_[n]-hx_[idx(i,j-1U,k)])/config_.dy);
        });
        enforce_boundary();
    }
}
void Maxwell3D::enforce_boundary(){
    const std::size_t nx=config_.nx,ny=config_.ny,nz=config_.nz;
    cfd::core::parallel_for(cell_count(),[&](std::size_t n){
        const std::size_t i=n%nx,j=(n/nx)%ny,k=n/(nx*ny);
        if(i==0U||i+1U==nx||j==0U||j+1U==ny||k==0U||k+1U==nz){
            if(config_.boundary==Boundary3D::pec){ex_[n]=0.0;ey_[n]=0.0;ez_[n]=0.0;}
            else {hx_[n]=0.0;hy_[n]=0.0;hz_[n]=0.0;}
        }
    });
}
double Maxwell3D::energy() const{
    return cfd::core::parallel_sum(cell_count(),[&](std::size_t n){
        return 0.5*(epsilon_*(ex_[n]*ex_[n]+ey_[n]*ey_[n]+ez_[n]*ez_[n])+mu_*(hx_[n]*hx_[n]+hy_[n]*hy_[n]+hz_[n]*hz_[n]))*config_.dx*config_.dy*config_.dz;
    });
}
double Maxwell3D::max_field() const{
    double m=0.0; for(std::size_t n=0;n<cell_count();++n) m=std::max(m,std::sqrt(ex_[n]*ex_[n]+ey_[n]*ey_[n]+ez_[n]*ez_[n])); return m;
}
} // namespace cfd::fdtd
