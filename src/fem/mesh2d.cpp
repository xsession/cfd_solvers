#include "cfd/fem/mesh2d.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fem {

void Mesh2D::validate() const {
    if (nodes.empty() || triangles.empty()) throw std::invalid_argument("FEM mesh must contain nodes/elements");
    if (boundary_node.size() != nodes.size()) throw std::invalid_argument("FEM boundary-node mask size mismatch");
    for (const auto& edge : boundary_edges) {
        if (edge.node[0] >= nodes.size() || edge.node[1] >= nodes.size()) {
            throw std::out_of_range("FEM boundary edge node index out of range");
        }
        if (edge.node[0] == edge.node[1]) throw std::invalid_argument("degenerate FEM boundary edge");
    }
    for (const auto& tri : triangles) {
        for (const auto n : tri.node) if (n >= nodes.size()) throw std::out_of_range("FEM triangle node index out of range");
        const auto& a = nodes[tri.node[0]];
        const auto& b = nodes[tri.node[1]];
        const auto& c = nodes[tri.node[2]];
        const double twice_area = (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
        if (!(std::abs(twice_area) > 1.0e-30)) throw std::invalid_argument("degenerate FEM triangle");
    }
}

Mesh2D make_rectangle_tri_mesh(std::size_t nx,
                               std::size_t ny,
                               double width,
                               double height) {
    if (nx == 0U || ny == 0U || !(width > 0.0) || !(height > 0.0)) {
        throw std::invalid_argument("invalid rectangle FEM mesh dimensions");
    }
    Mesh2D mesh;
    mesh.nodes.resize((nx+1U)*(ny+1U));
    mesh.boundary_node.resize(mesh.nodes.size(), 0U);
    const double dx = width/static_cast<double>(nx);
    const double dy = height/static_cast<double>(ny);
    const auto id = [nx](std::size_t i, std::size_t j) { return j*(nx+1U)+i; };
    for (std::size_t j=0; j<=ny; ++j) {
        for (std::size_t i=0; i<=nx; ++i) {
            const std::size_t n=id(i,j);
            mesh.nodes[n] = {static_cast<double>(i)*dx, static_cast<double>(j)*dy};
            mesh.boundary_node[n] = static_cast<unsigned char>(i==0U || i==nx || j==0U || j==ny);
        }
    }
    mesh.triangles.reserve(2U*nx*ny);
    for (std::size_t j=0; j<ny; ++j) {
        for (std::size_t i=0; i<nx; ++i) {
            const auto n00=id(i,j), n10=id(i+1U,j), n01=id(i,j+1U), n11=id(i+1U,j+1U);
            mesh.triangles.push_back({{n00,n10,n11}});
            mesh.triangles.push_back({{n00,n11,n01}});
        }
    }
    // Patch ids: 0=left, 1=right, 2=bottom, 3=top. Edge node order is
    // counter-clockwise around the rectangular boundary.
    mesh.boundary_edges.reserve(2U*(nx+ny));
    for (std::size_t j=0; j<ny; ++j) mesh.boundary_edges.push_back({{id(0U,j+1U),id(0U,j)},0});
    for (std::size_t j=0; j<ny; ++j) mesh.boundary_edges.push_back({{id(nx,j),id(nx,j+1U)},1});
    for (std::size_t i=0; i<nx; ++i) mesh.boundary_edges.push_back({{id(i,0U),id(i+1U,0U)},2});
    for (std::size_t i=0; i<nx; ++i) mesh.boundary_edges.push_back({{id(i+1U,ny),id(i,ny)},3});
    mesh.validate();
    return mesh;
}

} // namespace cfd::fem
