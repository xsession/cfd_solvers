#include "cfd/solvers/fem/heat2d.hpp"

#include "cfd/core/csr_matrix.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace cfd::fem {
namespace {
struct Geometry {
    double area{};
    std::array<std::array<double,2>,3> grad{};
    Node2 centroid{};
};
Geometry geometry(const Mesh2D& mesh,const Tri3& tri) {
    const auto& a=mesh.nodes[tri.node[0]]; const auto& b=mesh.nodes[tri.node[1]]; const auto& c=mesh.nodes[tri.node[2]];
    const double det=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    const double area=0.5*std::abs(det); if (!(area>0.0)) throw std::runtime_error("degenerate heat triangle");
    const double inv=1.0/det;
    return {area, {{{(b.y-c.y)*inv,(c.x-b.x)*inv},{(c.y-a.y)*inv,(a.x-c.x)*inv},{(a.y-b.y)*inv,(b.x-a.x)*inv}}},
            {(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0}};
}
} // namespace

Heat2D::Heat2D(Mesh2D mesh, Heat2DConfig config)
    : mesh_(std::move(mesh)), config_(config), temperature_(mesh_.node_count(),0.0),
      boundary_temperature_([](Node2,double){return 0.0;}) {
    mesh_.validate();
    if (!(config_.conductivity>0.0) || !(config_.volumetric_heat_capacity>0.0) || !(config_.dt>0.0)
        || config_.max_iterations==0U || !(config_.relative_tolerance>0.0)) {
        throw std::invalid_argument("invalid Heat2D controls");
    }
}
void Heat2D::initialize(const std::function<double(Node2)>& temperature) {
    if (!temperature) throw std::invalid_argument("Heat2D requires an initial temperature");
    for (std::size_t i=0;i<mesh_.node_count();++i) temperature_[i]=temperature(mesh_.nodes[i]);
    time_=0.0;
}
void Heat2D::set_dirichlet(const std::function<double(Node2,double)>& boundary_temperature) {
    if (!boundary_temperature) throw std::invalid_argument("Heat2D requires a boundary temperature");
    boundary_temperature_=boundary_temperature;
}
void Heat2D::run(std::size_t steps,const std::function<double(Node2,double)>& source) {
    for (std::size_t i=0;i<steps;++i) step(source);
}
void Heat2D::step(const std::function<double(Node2,double)>& source) {
    step_impl([&](std::size_t,Node2 x,double t){ return source ? source(x,t) : 0.0; });
}
void Heat2D::step_element_source(std::span<const double> element_source) {
    if (element_source.size()!=mesh_.element_count()) throw std::invalid_argument("Heat2D element source size mismatch");
    for (double q : element_source) if (!std::isfinite(q)) throw std::invalid_argument("non-finite Heat2D source");
    step_impl([&](std::size_t element,Node2,double){ return element_source[element]; });
}
void Heat2D::step_impl(const ElementSource& source) {
    const std::size_t n=mesh_.node_count();
    std::vector<std::size_t> free_index(n,static_cast<std::size_t>(-1));
    std::size_t free_count=0U;
    const double next_time=time_+config_.dt;
    auto next_temperature=temperature_;
    for (std::size_t i=0;i<n;++i) {
        if (mesh_.boundary_node[i]) next_temperature[i]=boundary_temperature_(mesh_.nodes[i],next_time);
        else free_index[i]=free_count++;
    }
    cfd::core::CsrBuilder builder(free_count,free_count);
    std::vector<double> rhs(free_count,0.0);
    for (std::size_t element=0;element<mesh_.triangles.size();++element) {
        const auto& tri=mesh_.triangles[element];
        const auto g=geometry(mesh_,tri);
        const double q=source(element,g.centroid,next_time);
        for (std::size_t i=0;i<3U;++i) {
            const auto ni=tri.node[i]; if (mesh_.boundary_node[ni]) continue;
            const auto row=free_index[ni];
            rhs[row]+=q*g.area/3.0;
            for (std::size_t j=0;j<3U;++j) {
                const auto nj=tri.node[j];
                const double mass=config_.volumetric_heat_capacity*g.area*(i==j?2.0:1.0)/12.0;
                const double stiff=config_.conductivity*g.area*(g.grad[i][0]*g.grad[j][0]+g.grad[i][1]*g.grad[j][1]);
                const double aij=mass/config_.dt+stiff;
                rhs[row]+=mass/config_.dt*temperature_[nj];
                if (mesh_.boundary_node[nj]) rhs[row]-=aij*next_temperature[nj];
                else builder.add(row,free_index[nj],aij);
            }
        }
    }
    const auto matrix=builder.build();
    cfd::core::JacobiPreconditioner jacobi(matrix.diagonal());
    cfd::core::KrylovWorkspace workspace;
    std::vector<double> x(free_count,0.0);
    for (std::size_t i=0;i<n;++i) if (!mesh_.boundary_node[i]) x[free_index[i]]=temperature_[i];
    linear_result_=cfd::core::preconditioned_conjugate_gradient(
        rhs,x,[&](std::span<const double> v,std::span<double> out){matrix.multiply(v,out);},
        [&](std::span<const double> r,std::span<double> z){jacobi(r,z);},workspace,
        config_.max_iterations,config_.relative_tolerance);
    if (!linear_result_.converged) throw std::runtime_error("Heat2D PCG did not converge");
    for (std::size_t i=0;i<n;++i) if (!mesh_.boundary_node[i]) next_temperature[i]=x[free_index[i]];
    temperature_.swap(next_temperature);
    time_=next_time;
}
} // namespace cfd::fem
