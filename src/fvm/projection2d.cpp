#include "cfd/solvers/fvm/projection2d.hpp"
#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::fvm {

Projection2D::Projection2D(Projection2DConfig config)
    : config_(config),
      dx_(config.lx/static_cast<double>(config.nx)),
      dy_(config.ly/static_cast<double>(config.ny)),
      u_((config.nx+1U)*config.ny,0.0),
      v_(config.nx*(config.ny+1U),0.0),
      p_(config.nx*config.ny,0.0),
      rhs_(config.nx*config.ny,0.0),
      pressure_rhs_(config.nx*config.ny,0.0) {
    if (config_.nx<4 || config_.ny<4 || !(config_.lx>0.0) || !(config_.ly>0.0) ||
        !(config_.density>0.0) || !(config_.dt>0.0) || config_.pressure_iterations==0 ||
        !(config_.pressure_tolerance>0.0)) {
        throw std::invalid_argument("Projection2D requires nx,ny>=4 and positive physical/numerical parameters");
    }
}

std::size_t Projection2D::u_index(std::size_t i,std::size_t j) const noexcept { return j*(config_.nx+1U)+i; }
std::size_t Projection2D::v_index(std::size_t i,std::size_t j) const noexcept { return j*config_.nx+i; }
std::size_t Projection2D::p_index(std::size_t i,std::size_t j) const noexcept { return j*config_.nx+i; }

void Projection2D::enforce_periodic_faces() {
    for (std::size_t j=0;j<config_.ny;++j) u_[u_index(config_.nx,j)]=u_[u_index(0,j)];
    for (std::size_t i=0;i<config_.nx;++i) v_[v_index(i,config_.ny)]=v_[v_index(i,0)];
}

void Projection2D::initialize_divergent(double amplitude) {
    const double kx=2.0*std::numbers::pi/config_.lx;
    const double ky=2.0*std::numbers::pi/config_.ly;
    for (std::size_t j=0;j<config_.ny;++j) {
        const double y=(static_cast<double>(j)+0.5)*dy_;
        (void)y;
        for (std::size_t i=0;i<=config_.nx;++i) {
            const double x=static_cast<double>(i)*dx_;
            u_[u_index(i,j)]=amplitude*std::sin(kx*x);
        }
    }
    for (std::size_t j=0;j<=config_.ny;++j) {
        const double y=static_cast<double>(j)*dy_;
        for (std::size_t i=0;i<config_.nx;++i) {
            const double x=(static_cast<double>(i)+0.5)*dx_;
            (void)x;
            v_[v_index(i,j)]=amplitude*std::sin(ky*y);
        }
    }
    std::fill(p_.begin(),p_.end(),0.0);
    enforce_periodic_faces();
}

void Projection2D::initialize_taylor_green(double amplitude) {
    const double kx=2.0*std::numbers::pi/config_.lx;
    const double ky=2.0*std::numbers::pi/config_.ly;
    for (std::size_t j=0;j<config_.ny;++j) {
        const double y=(static_cast<double>(j)+0.5)*dy_;
        for (std::size_t i=0;i<=config_.nx;++i) {
            const double x=static_cast<double>(i)*dx_;
            u_[u_index(i,j)]=amplitude*std::sin(kx*x)*std::cos(ky*y);
        }
    }
    for (std::size_t j=0;j<=config_.ny;++j) {
        const double y=static_cast<double>(j)*dy_;
        for (std::size_t i=0;i<config_.nx;++i) {
            const double x=(static_cast<double>(i)+0.5)*dx_;
            v_[v_index(i,j)]=-amplitude*std::cos(kx*x)*std::sin(ky*y);
        }
    }
    std::fill(p_.begin(),p_.end(),0.0);
    enforce_periodic_faces();
}

void Projection2D::project() {
    const double inv_dx=1.0/dx_;
    const double inv_dy=1.0/dy_;
    const double scale=config_.density/config_.dt;
    cfd::core::parallel_for(config_.nx * config_.ny, [&](std::size_t n) {
        const std::size_t i = n % config_.nx;
        const std::size_t j = n / config_.nx;
        const double div=(u_[u_index(i+1U,j)]-u_[u_index(i,j)])*inv_dx +
                         (v_[v_index(i,j+1U)]-v_[v_index(i,j)])*inv_dy;
        rhs_[p_index(i,j)]=scale*div;
    });

    const double idx2=inv_dx*inv_dx;
    const double idy2=inv_dy*inv_dy;
    const std::size_t cells=config_.nx*config_.ny;
    std::fill(p_.begin(),p_.end(),0.0);
    cfd::core::parallel_for(cells,[&](std::size_t n){pressure_rhs_[n]=-rhs_[n];});

    auto apply_negative_laplacian = [&](std::span<const double> in, std::span<double> out) {
        cfd::core::parallel_for(cells, [&](std::size_t n) {
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
        std::span<double>(p_.data(),p_.size()),
        apply_negative_laplacian,pressure_workspace_,config_.pressure_iterations,config_.pressure_tolerance);

    double mean=cfd::core::parallel_sum(p_.size(), [&](std::size_t n) { return p_[n]; });
    mean/=static_cast<double>(p_.size());
    cfd::core::parallel_for(p_.size(), [&](std::size_t n) { p_[n]-=mean; });

    const double correction=config_.dt/config_.density;
    cfd::core::parallel_for(config_.nx * config_.ny, [&](std::size_t n) {
        const std::size_t i = n % config_.nx;
        const std::size_t j = n / config_.nx;
        const std::size_t im=(i+config_.nx-1U)%config_.nx;
        u_[u_index(i,j)]-=correction*(p_[p_index(i,j)]-p_[p_index(im,j)])*inv_dx;
    });
    cfd::core::parallel_for(config_.nx * config_.ny, [&](std::size_t n) {
        const std::size_t i = n % config_.nx;
        const std::size_t j = n / config_.nx;
        const std::size_t jm=(j+config_.ny-1U)%config_.ny;
        v_[v_index(i,j)]-=correction*(p_[p_index(i,j)]-p_[p_index(i,jm)])*inv_dy;
    });
    enforce_periodic_faces();
}

double Projection2D::divergence_l2() const {
    const double sum=cfd::core::parallel_sum(config_.nx * config_.ny, [&](std::size_t n) {
        const std::size_t i = n % config_.nx;
        const std::size_t j = n / config_.nx;
        const double div=(u_[u_index(i+1U,j)]-u_[u_index(i,j)])/dx_ +
                         (v_[v_index(i,j+1U)]-v_[v_index(i,j)])/dy_;
        return div*div;
    });
    return std::sqrt(sum/static_cast<double>(config_.nx*config_.ny));
}

double Projection2D::kinetic_energy() const {
    const double sum=cfd::core::parallel_sum(config_.nx * config_.ny, [&](std::size_t n) {
        const std::size_t i = n % config_.nx;
        const std::size_t j = n / config_.nx;
        const double uc=0.5*(u_[u_index(i,j)]+u_[u_index(i+1U,j)]);
        const double vc=0.5*(v_[v_index(i,j)]+v_[v_index(i,j+1U)]);
        return 0.5*(uc*uc+vc*vc);
    });
    return sum*dx_*dy_;
}

} // namespace cfd::fvm
