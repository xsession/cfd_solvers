#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <cstddef>
#include <span>
#include <vector>
#include <stdexcept>
namespace cfd::fvm {
[[nodiscard]] std::vector<std::size_t> morton_cell_order(const PolyMesh& mesh);
// Rebuilds the mesh with cells in the supplied permutation while preserving
// owner/neighbour orientation and rebuilding cell-face adjacency.
[[nodiscard]] PolyMesh reorder_cells(const PolyMesh& mesh, std::span<const std::size_t> order);
template <class T>
[[nodiscard]] std::vector<T> reorder_cell_field(std::span<const T> field, std::span<const std::size_t> order) {
    if (field.size() != order.size())
        throw std::invalid_argument("cell reorder size mismatch");
    std::vector<T> out(order.size());
    for (std::size_t i = 0; i < order.size(); ++i)
        out[i] = field[order[i]];
    return out;
}
} // namespace cfd::fvm
