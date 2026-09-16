#include "cfd/solvers/fem/poisson2d.hpp"

#include "cfd/core/csr_matrix.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cfd::fem {
namespace {

struct TriGeometry {
    double area{};
    std::array<std::array<double,2>,3> grad{};
    Node2 centroid{};
};

TriGeometry geometry(const Mesh2D& mesh, const Tri3& tri) {
    const auto& a=mesh.nodes[tri.node[0]];
    const auto& b=mesh.nodes[tri.node[1]];
    const auto& c=mesh.nodes[tri.node[2]];
    const double det=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    const double area=0.5*std::abs(det);
    if (!(area>0.0)) throw std::runtime_error("degenerate FEM triangle during assembly");
    const double inv_det=1.0/det;
    TriGeometry g;
    g.area=area;
    g.grad = {{{(b.y-c.y)*inv_det,(c.x-b.x)*inv_det},
               {(c.y-a.y)*inv_det,(a.x-c.x)*inv_det},
               {(a.y-b.y)*inv_det,(b.x-a.x)*inv_det}}};
    g.centroid={(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0};
    return g;
}

} // namespace

Poisson2D::Poisson2D(Mesh2D mesh, Poisson2DConfig config)
    : mesh_(std::move(mesh)), config_(config), solution_(mesh_.node_count(),0.0) {
    mesh_.validate();
    if (!(config_.diffusivity>0.0) || config_.max_iterations==0U || !(config_.relative_tolerance>0.0)) {
        throw std::invalid_argument("invalid Poisson2D controls");
    }
}

void Poisson2D::solve(const std::function<double(Node2)>& source,
                      const std::function<double(Node2)>& dirichlet_value) {
    const std::size_t n=mesh_.node_count();
    std::vector<std::size_t> free_index(n, static_cast<std::size_t>(-1));
    std::size_t free_count=0U;
    for (std::size_t i=0;i<n;++i) {
        if (mesh_.boundary_node[i]) solution_[i]=dirichlet_value(mesh_.nodes[i]);
        else free_index[i]=free_count++;
    }
    if (free_count==0U) throw std::runtime_error("Poisson2D mesh has no free nodes");

    cfd::core::CsrBuilder builder(free_count,free_count);
    std::vector<double> rhs(free_count,0.0);
    for (const auto& tri:mesh_.triangles) {
        const auto g=geometry(mesh_,tri);
        const double f=source(g.centroid);
        for (std::size_t i=0;i<3U;++i) {
            const auto ni=tri.node[i];
            if (mesh_.boundary_node[ni]) continue;
            const auto row=free_index[ni];
            rhs[row]+=f*g.area/3.0; // one-point centroid quadrature
            for (std::size_t j=0;j<3U;++j) {
                const auto nj=tri.node[j];
                const double kij=config_.diffusivity*g.area*(g.grad[i][0]*g.grad[j][0]+g.grad[i][1]*g.grad[j][1]);
                if (mesh_.boundary_node[nj]) rhs[row]-=kij*solution_[nj];
                else builder.add(row,free_index[nj],kij);
            }
        }
    }
    const auto matrix=builder.build();
    auto diagonal=matrix.diagonal();
    cfd::core::JacobiPreconditioner jacobi(diagonal);
    cfd::core::KrylovWorkspace workspace;
    std::vector<double> x(free_count,0.0);
    linear_result_=cfd::core::preconditioned_conjugate_gradient(
        rhs,x,
        [&](std::span<const double> v,std::span<double> out){matrix.multiply(v,out);},
        [&](std::span<const double> r,std::span<double> z){jacobi(r,z);},
        workspace,config_.max_iterations,config_.relative_tolerance);
    if (!linear_result_.converged) throw std::runtime_error("Poisson2D PCG did not converge");
    for (std::size_t i=0;i<n;++i) if (!mesh_.boundary_node[i]) solution_[i]=x[free_index[i]];
}

} // namespace cfd::fem
