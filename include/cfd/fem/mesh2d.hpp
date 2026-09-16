#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace cfd::fem {

struct Node2 {
    double x{};
    double y{};
};

struct Tri3 {
    std::array<std::size_t, 3> node{};
};

struct BoundaryEdge2 {
    std::array<std::size_t, 2> node{};
    int patch{};
};

struct Mesh2D {
    std::vector<Node2> nodes;
    std::vector<Tri3> triangles;
    std::vector<BoundaryEdge2> boundary_edges;
    std::vector<unsigned char> boundary_node;

    [[nodiscard]] std::size_t node_count() const noexcept { return nodes.size(); }
    [[nodiscard]] std::size_t element_count() const noexcept { return triangles.size(); }
    void validate() const;
};

[[nodiscard]] Mesh2D make_rectangle_tri_mesh(std::size_t nx,
                                             std::size_t ny,
                                             double width = 1.0,
                                             double height = 1.0);

} // namespace cfd::fem
