#pragma once

#include "cfd/fem/reference_element.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace cfd::fem {

struct Tet4 {
    std::array<std::size_t,4> node{};
};

struct Mesh3D {
    std::vector<Point3> nodes;
    std::vector<Tet4> tetrahedra;
    std::vector<unsigned char> boundary_node;

    [[nodiscard]] std::size_t node_count() const noexcept { return nodes.size(); }
    [[nodiscard]] std::size_t element_count() const noexcept { return tetrahedra.size(); }
    void validate() const;
};

[[nodiscard]] Mesh3D make_box_tet_mesh(std::size_t nx,
                                       std::size_t ny,
                                       std::size_t nz,
                                       double width=1.0,
                                       double height=1.0,
                                       double depth=1.0);

} // namespace cfd::fem
