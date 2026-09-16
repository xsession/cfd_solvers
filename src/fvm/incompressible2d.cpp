#include "cfd/solvers/fvm/incompressible2d.hpp"
#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <span>
#include <stdexcept>

namespace cfd::fvm {
namespace {

std::size_t wrap(std::ptrdiff_t value, std::size_t n) noexcept {
    const auto nn = static_cast<std::ptrdiff_t>(n);
    value %= nn;
    if (value < 0) value += nn;
    return static_cast<std::size_t>(value);
}

} // namespace

Incompressible2D::Incompressible2D(Incompressible2DConfig config)
    : config_(config),
      dx_(config.lx/static_cast<double>(config.nx)),
      dy_(config.ly/static_cast<double>(config.ny)),
      u_((config.nx+1U)*config.ny, 0.0),
      v_(config.nx*(config.ny+1U), 0.0),
      u_star_(u_.size(), 0.0),
      v_star_(v_.size(), 0.0),
      p_(config.nx*config.ny, 0.0),
      pressure_rhs_(p_.size(), 0.0) {
    if (config_.nx < 8U || config_.ny < 8U || !(config_.lx > 0.0) || !(config_.ly > 0.0) ||
        !(config_.density > 0.0) || !(config_.kinematic_viscosity >= 0.0) || !(config_.dt > 0.0) ||
        config_.pressure_iterations == 0U || !(config_.pressure_tolerance > 0.0)) {
        throw std::invalid_argument("Incompressible2D requires nx,ny>=8 and valid positive numerical parameters");
    }
}

std::size_t Incompressible2D::u_index(std::size_t i, std::size_t j) const noexcept {
    return j*(config_.nx+1U)+i;
}
std::size_t Incompressible2D::v_index(std::size_t i, std::size_t j) const noexcept {
    return j*config_.nx+i;
}
std::size_t Incompressible2D::p_index(std::size_t i, std::size_t j) const noexcept {
    return j*config_.nx+i;
}

double Incompressible2D::u_periodic(std::ptrdiff_t i, std::ptrdiff_t j) const noexcept {
    return u_[u_index(wrap(i,config_.nx),wrap(j,config_.ny))];
}
double Incompressible2D::v_periodic(std::ptrdiff_t i, std::ptrdiff_t j) const noexcept {
    return v_[v_index(wrap(i,config_.nx),wrap(j,config_.ny))];
}

void Incompressible2D::enforce_periodic(cfd::core::AlignedVector<double>& u,
                                        cfd::core::AlignedVector<double>& v) {
    for (std::size_t j=0;j<config_.ny;++j) u[u_index(config_.nx,j)] = u[u_index(0,j)];
    for (std::size_t i=0;i<config_.nx;++i) v[v_index(i,config_.ny)] = v[v_index(i,0)];
}

void Incompressible2D::initialize_taylor_green(double amplitude) {
    const double kx=2.0*std::numbers::pi/config_.lx;
    const double ky=2.0*std::numbers::pi/config_.ly;
    cfd::core::parallel_for(config_.nx*config_.ny,[&](std::size_t n) {
        const std::size_t i=n%config_.nx;
        const std::size_t j=n/config_.nx;
        const double xu=static_cast<double>(i)*dx_;
        const double yu=(static_cast<double>(j)+0.5)*dy_;
        const double xv=(static_cast<double>(i)+0.5)*dx_;
        const double yv=static_cast<double>(j)*dy_;
        u_[u_index(i,j)]= amplitude*std::sin(kx*xu)*std::cos(ky*yu);
        v_[v_index(i,j)]=-amplitude*std::cos(kx*xv)*std::sin(ky*yv);
    });
    enforce_periodic(u_,v_);
    std::fill(p_.begin(),p_.end(),0.0);
    time_=0.0;
    steps_=0U;
}

void Incompressible2D::momentum_predictor() {
    const double inv_dx=1.0/dx_;
    const double inv_dy=1.0/dy_;
    const double inv_dx2=inv_dx*inv_dx;
    const double inv_dy2=inv_dy*inv_dy;
    const double nu=config_.kinematic_viscosity;

    cfd::core::parallel_for(config_.nx*config_.ny,[&](std::size_t n) {
        const std::ptrdiff_t i=static_cast<std::ptrdiff_t>(n%config_.nx);
        const std::ptrdiff_t j=static_cast<std::ptrdiff_t>(n/config_.nx);
        const double uc=u_periodic(i,j);
        const double ue=0.5*(uc+u_periodic(i+1,j));
        const double uw=0.5*(u_periodic(i-1,j)+uc);
        const double un=0.5*(uc+u_periodic(i,j+1));
        const double us=0.5*(u_periodic(i,j-1)+uc);
        const double vn=0.5*(v_periodic(i-1,j+1)+v_periodic(i,j+1));
        const double vs=0.5*(v_periodic(i-1,j)+v_periodic(i,j));
        const double convection=(ue*ue-uw*uw)*inv_dx+(un*vn-us*vs)*inv_dy;
        const double lap=(u_periodic(i+1,j)-2.0*uc+u_periodic(i-1,j))*inv_dx2+
                         (u_periodic(i,j+1)-2.0*uc+u_periodic(i,j-1))*inv_dy2;
        u_star_[u_index(static_cast<std::size_t>(i),static_cast<std::size_t>(j))]=
            uc+config_.dt*(-convection+nu*lap);
    });

    cfd::core::parallel_for(config_.nx*config_.ny,[&](std::size_t n) {
        const std::ptrdiff_t i=static_cast<std::ptrdiff_t>(n%config_.nx);
        const std::ptrdiff_t j=static_cast<std::ptrdiff_t>(n/config_.nx);
        const double vc=v_periodic(i,j);
        const double ve=0.5*(vc+v_periodic(i+1,j));
        const double vw=0.5*(v_periodic(i-1,j)+vc);
        const double vn=0.5*(vc+v_periodic(i,j+1));
        const double vs=0.5*(v_periodic(i,j-1)+vc);
        const double ue=0.5*(u_periodic(i+1,j-1)+u_periodic(i+1,j));
        const double uw=0.5*(u_periodic(i,j-1)+u_periodic(i,j));
        const double convection=(ue*ve-uw*vw)*inv_dx+(vn*vn-vs*vs)*inv_dy;
        const double lap=(v_periodic(i+1,j)-2.0*vc+v_periodic(i-1,j))*inv_dx2+
                         (v_periodic(i,j+1)-2.0*vc+v_periodic(i,j-1))*inv_dy2;
        v_star_[v_index(static_cast<std::size_t>(i),static_cast<std::size_t>(j))]=
            vc+config_.dt*(-convection+nu*lap);
    });
    enforce_periodic(u_star_,v_star_);
}

void Incompressible2D::project_star_velocity() {
    const double inv_dx=1.0/dx_;
    const double inv_dy=1.0/dy_;
    const double scale=config_.density/config_.dt;
    const std::size_t cells=config_.nx*config_.ny;
    cfd::core::parallel_for(cells,[&](std::size_t n) {
        const std::size_t i=n%config_.nx;
        const std::size_t j=n/config_.nx;
        const double div=(u_star_[u_index(i+1U,j)]-u_star_[u_index(i,j)])*inv_dx+
                         (v_star_[v_index(i,j+1U)]-v_star_[v_index(i,j)])*inv_dy;
        pressure_rhs_[n]=-scale*div;
    });

    // Enforce the compatible mean-zero subspace explicitly to prevent roundoff
    // from exciting the periodic Poisson nullspace.
    double rhs_mean=cfd::core::parallel_sum(cells,[&](std::size_t n){return pressure_rhs_[n];});
    rhs_mean/=static_cast<double>(cells);
    cfd::core::parallel_for(cells,[&](std::size_t n){pressure_rhs_[n]-=rhs_mean;});

    const double idx2=inv_dx*inv_dx;
    const double idy2=inv_dy*inv_dy;
    std::fill(p_.begin(),p_.end(),0.0);
    auto apply_negative_laplacian=[&](std::span<const double> in,std::span<double> out) {
        cfd::core::parallel_for(cells,[&](std::size_t n) {
            const std::size_t i=n%config_.nx;
            const std::size_t j=n/config_.nx;
            const std::size_t im=(i+config_.nx-1U)%config_.nx;
            const std::size_t ip=(i+1U)%config_.nx;
            const std::size_t jm=(j+config_.ny-1U)%config_.ny;
            const std::size_t jp=(j+1U)%config_.ny;
            out[n]=2.0*(idx2+idy2)*in[n]
                  -idx2*(in[p_index(im,j)]+in[p_index(ip,j)])
                  -idy2*(in[p_index(i,jm)]+in[p_index(i,jp)]);
        });
    };
    pressure_result_=cfd::core::conjugate_gradient(
        std::span<const double>(pressure_rhs_.data(),pressure_rhs_.size()),
        std::span<double>(p_.data(),p_.size()),apply_negative_laplacian,
        pressure_workspace_,config_.pressure_iterations,config_.pressure_tolerance);

    double mean=cfd::core::parallel_sum(cells,[&](std::size_t n){return p_[n];});
    mean/=static_cast<double>(cells);
    cfd::core::parallel_for(cells,[&](std::size_t n){p_[n]-=mean;});

    const double correction=config_.dt/config_.density;
    cfd::core::parallel_for(cells,[&](std::size_t n) {
        const std::size_t i=n%config_.nx;
        const std::size_t j=n/config_.nx;
        const std::size_t im=(i+config_.nx-1U)%config_.nx;
        u_[u_index(i,j)]=u_star_[u_index(i,j)]-
            correction*(p_[p_index(i,j)]-p_[p_index(im,j)])*inv_dx;
        const std::size_t jm=(j+config_.ny-1U)%config_.ny;
        v_[v_index(i,j)]=v_star_[v_index(i,j)]-
            correction*(p_[p_index(i,j)]-p_[p_index(i,jm)])*inv_dy;
    });
    enforce_periodic(u_,v_);
}

void Incompressible2D::step() {
    momentum_predictor();
    project_star_velocity();
    time_+=config_.dt;
    ++steps_;
}

void Incompressible2D::run(std::size_t count) {
    for (std::size_t s=0;s<count;++s) step();
}

double Incompressible2D::divergence_l2() const {
    const double sum=cfd::core::parallel_sum(config_.nx*config_.ny,[&](std::size_t n) {
        const std::size_t i=n%config_.nx;
        const std::size_t j=n/config_.nx;
        const double div=(u_[u_index(i+1U,j)]-u_[u_index(i,j)])/dx_+
                         (v_[v_index(i,j+1U)]-v_[v_index(i,j)])/dy_;
        return div*div;
    });
    return std::sqrt(sum/static_cast<double>(config_.nx*config_.ny));
}

double Incompressible2D::kinetic_energy() const {
    const double sum=cfd::core::parallel_sum(config_.nx*config_.ny,[&](std::size_t n) {
        const std::size_t i=n%config_.nx;
        const std::size_t j=n/config_.nx;
        const double uc=0.5*(u_[u_index(i,j)]+u_[u_index(i+1U,j)]);
        const double vc=0.5*(v_[v_index(i,j)]+v_[v_index(i,j+1U)]);
        return 0.5*(uc*uc+vc*vc);
    });
    return sum*dx_*dy_;
}

double Incompressible2D::taylor_green_velocity_error(double initial_amplitude) const {
    const double kx=2.0*std::numbers::pi/config_.lx;
    const double ky=2.0*std::numbers::pi/config_.ly;
    const double decay=std::exp(-config_.kinematic_viscosity*(kx*kx+ky*ky)*time_);
    const double amplitude=initial_amplitude*decay;
    const double sum=cfd::core::parallel_sum(config_.nx*config_.ny,[&](std::size_t n) {
        const std::size_t i=n%config_.nx;
        const std::size_t j=n/config_.nx;
        const double xu=static_cast<double>(i)*dx_;
        const double yu=(static_cast<double>(j)+0.5)*dy_;
        const double xv=(static_cast<double>(i)+0.5)*dx_;
        const double yv=static_cast<double>(j)*dy_;
        const double exact_u= amplitude*std::sin(kx*xu)*std::cos(ky*yu);
        const double exact_v=-amplitude*std::cos(kx*xv)*std::sin(ky*yv);
        const double du=u_[u_index(i,j)]-exact_u;
        const double dv=v_[v_index(i,j)]-exact_v;
        return du*du+dv*dv;
    });
    return std::sqrt(sum/(2.0*static_cast<double>(config_.nx*config_.ny)));
}

} // namespace cfd::fvm
