#include "cfd/fvm/reordering.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    try {
        using namespace cfd::fvm;
        const auto mesh = make_cartesian_hexa_mesh(4U, 3U, 1U, 1.0, 0.75, 1.0);
        const auto order = morton_cell_order(mesh);
        const auto reordered = reorder_cells(mesh, order);
        if (reordered.cell_count() != mesh.cell_count() || reordered.face_count() != mesh.face_count())
            throw std::runtime_error("reordering changed mesh cardinality");
        for (const auto& face : reordered.faces()) {
            if (face.owner >= reordered.cell_count() || (!face.boundary() && face.neighbour >= reordered.cell_count()))
                throw std::runtime_error("reordering produced invalid face connectivity");
        }
        std::vector<int> field(mesh.cell_count());
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] = static_cast<int>(i);
        const auto moved = reorder_cell_field<int>(field, order);
        if (moved.front() != static_cast<int>(order.front()))
            throw std::runtime_error("field reorder mismatch");
        std::cout << "connectivity-preserving FVM reorder regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FVM reorder regression failed: " << error.what() << '\n';
        return 1;
    }
}
