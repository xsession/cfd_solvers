#include "cfd/solvers/fem/poisson3d.hpp"

#include "cfd/core/csr_matrix.hpp"
#include "cfd/fem/reference_element.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace cfd::fem {

Poisson3D::Poisson3D(Mesh3D mesh,Poisson3DConfig config)
    :mesh_(std::move(mesh)),config_(config),solution_(mesh_.node_count(),0.0) {
    mesh_.validate();
    if(!(config_.diffusivity>0.0)||config_.max_iterations==0U||!(config_.relative_tolerance>0.0)) throw std::invalid_argument("invalid Poisson3D controls");
}

void Poisson3D::solve(const std::function<double(Point3)>& source,
                      const std::function<double(Point3)>& dirichlet_value) {
    if(!source||!dirichlet_value) throw std::invalid_argument("Poisson3D callback missing");
    const std::size_t n=mesh_.node_count();
    std::vector<std::size_t> map(n,static_cast<std::size_t>(-1));
    std::size_t nf=0U;
    for(std::size_t i=0;i<n;++i) {
        if(mesh_.boundary_node[i]) solution_[i]=dirichlet_value(mesh_.nodes[i]); else map[i]=nf++;
    }
    if(nf==0U) throw std::runtime_error("Poisson3D mesh has no free nodes");
    cfd::core::CsrBuilder builder(nf,nf);
    std::vector<double> rhs(nf,0.0);
    const auto qrule=gaussian_quadrature(ElementType::tet4,2U);
    for(const auto&t:mesh_.tetrahedra) {
        std::array<Point3,4> nodes{{mesh_.nodes[t.node[0]],mesh_.nodes[t.node[1]],mesh_.nodes[t.node[2]],mesh_.nodes[t.node[3]]}};
        double ke[4][4]{}; double fe[4]{};
        for(const auto&qp:qrule) {
            const auto shape=evaluate_shape(ElementType::tet4,qp.point);
            const auto iso=evaluate_isoparametric(ElementType::tet4,nodes,qp.point);
            const double f=source(iso.point); const double w=qp.weight*iso.measure;
            for(std::size_t i=0;i<4U;++i) {
                fe[i]+=f*shape.value[i]*w;
                for(std::size_t j=0;j<4U;++j) {
                    const auto&gi=iso.gradient_physical[i];const auto&gj=iso.gradient_physical[j];
                    ke[i][j]+=config_.diffusivity*(gi.x*gj.x+gi.y*gj.y+gi.z*gj.z)*w;
                }
            }
        }
        for(std::size_t i=0;i<4U;++i) {
            const auto ni=t.node[i]; if(mesh_.boundary_node[ni]) continue; const auto row=map[ni]; rhs[row]+=fe[i];
            for(std::size_t j=0;j<4U;++j) {
                const auto nj=t.node[j]; if(mesh_.boundary_node[nj]) rhs[row]-=ke[i][j]*solution_[nj]; else builder.add(row,map[nj],ke[i][j]);
            }
        }
    }
    const auto A=builder.build(); cfd::core::JacobiPreconditioner jacobi(A.diagonal()); cfd::core::KrylovWorkspace ws; std::vector<double>x(nf,0.0);
    linear_result_=cfd::core::preconditioned_conjugate_gradient(rhs,x,
        [&](std::span<const double>v,std::span<double>o){A.multiply(v,o);},
        [&](std::span<const double>r,std::span<double>z){jacobi(r,z);},ws,config_.max_iterations,config_.relative_tolerance);
    if(!linear_result_.converged) throw std::runtime_error("Poisson3D PCG did not converge");
    for(std::size_t i=0;i<n;++i) if(!mesh_.boundary_node[i]) solution_[i]=x[map[i]];
}

} // namespace cfd::fem
