#include "cfd/solvers/fem/electromagnetics2d.hpp"

#include "cfd/fem/reference_element.hpp"

#include <array>
#include <stdexcept>

namespace cfd::fem {
namespace {

Vector2 triangle_gradient(const Mesh2D& mesh,
                          const Tri3& tri,
                          const std::vector<double>& value) {
    std::array<Point3, 3> nodes{{
        {mesh.nodes[tri.node[0]].x, mesh.nodes[tri.node[0]].y, 0.0},
        {mesh.nodes[tri.node[1]].x, mesh.nodes[tri.node[1]].y, 0.0},
        {mesh.nodes[tri.node[2]].x, mesh.nodes[tri.node[2]].y, 0.0}}};
    const auto iso = evaluate_isoparametric(ElementType::tri3, nodes, {1.0/3.0,1.0/3.0,0.0});
    Vector2 g{};
    for (std::size_t i=0;i<3U;++i) {
        g.x += value[tri.node[i]] * iso.gradient_physical[i].x;
        g.y += value[tri.node[i]] * iso.gradient_physical[i].y;
    }
    return g;
}

Node2 centroid(const Mesh2D& mesh, const Tri3& tri) {
    const auto&a=mesh.nodes[tri.node[0]],&b=mesh.nodes[tri.node[1]],&c=mesh.nodes[tri.node[2]];
    return {(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0};
}

} // namespace

Electrostatics2D::Electrostatics2D(Mesh2D mesh, ScalarDiffusion2DConfig config)
    : solver_(std::move(mesh), config) {}

void Electrostatics2D::set_boundary(int patch, ScalarBoundaryCondition2D condition) {
    solver_.set_boundary(patch, std::move(condition));
}

void Electrostatics2D::solve(const std::function<double(Node2)>& relative_permittivity,
                             const std::function<double(Node2)>& charge_over_epsilon0) {
    if (!relative_permittivity || !charge_over_epsilon0) throw std::invalid_argument("electrostatic coefficient callback missing");
    solver_.solve(relative_permittivity, [](Node2){return 0.0;}, charge_over_epsilon0);
}

std::vector<Vector2> Electrostatics2D::element_electric_field() const {
    std::vector<Vector2> field;
    field.reserve(mesh().triangles.size());
    for (const auto& tri : mesh().triangles) {
        const auto g = triangle_gradient(mesh(), tri, potential());
        field.push_back({-g.x,-g.y});
    }
    return field;
}

DCConduction2D::DCConduction2D(Mesh2D mesh, ScalarDiffusion2DConfig config)
    : solver_(std::move(mesh), config) {}

void DCConduction2D::set_boundary(int patch, ScalarBoundaryCondition2D condition) {
    solver_.set_boundary(patch, std::move(condition));
}

void DCConduction2D::solve(const std::function<double(Node2)>& conductivity,
                           const std::function<double(Node2)>& current_source) {
    if (!conductivity) throw std::invalid_argument("DC conductivity callback missing");
    const auto zero=[](Node2){return 0.0;};
    solver_.solve(conductivity, zero, current_source ? current_source : std::function<double(Node2)>(zero));
}

std::vector<Vector2> DCConduction2D::element_current_density(
    const std::function<double(Node2)>& conductivity) const {
    if (!conductivity) throw std::invalid_argument("DC conductivity callback missing");
    std::vector<Vector2> current;
    current.reserve(mesh().triangles.size());
    for (const auto& tri : mesh().triangles) {
        const auto g = triangle_gradient(mesh(), tri, potential());
        const double sigma = conductivity(centroid(mesh(), tri));
        if (!(sigma > 0.0)) throw std::runtime_error("invalid DC conductivity");
        current.push_back({-sigma*g.x,-sigma*g.y});
    }
    return current;
}

} // namespace cfd::fem
