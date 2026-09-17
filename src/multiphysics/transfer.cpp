#include "cfd/multiphysics/transfer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::multiphysics {

void CellGrid1D::validate() const {
    if (edge.size() < 2U) throw std::invalid_argument("1-D transfer grid requires at least one cell");
    for (std::size_t i = 1U; i < edge.size(); ++i) {
        if (!(edge[i] > edge[i - 1U])) throw std::invalid_argument("1-D transfer edges must be strictly increasing");
    }
}

std::vector<double> conservative_cell_average_transfer(const CellGrid1D& source,
                                                        std::span<const double> source_average,
                                                        const CellGrid1D& target) {
    source.validate();
    target.validate();
    if (source_average.size() != source.cell_count()) throw std::invalid_argument("source field/grid size mismatch");
    const double scale = std::max({1.0, std::abs(source.edge.front()), std::abs(source.edge.back())});
    if (std::abs(source.edge.front() - target.edge.front()) > 1.0e-12 * scale
        || std::abs(source.edge.back() - target.edge.back()) > 1.0e-12 * scale) {
        throw std::invalid_argument("conservative 1-D transfer requires identical physical extent");
    }

    std::vector<double> target_average(target.cell_count(), 0.0);
    std::size_t source_cell = 0U;
    for (std::size_t t = 0U; t < target.cell_count(); ++t) {
        const double left = target.edge[t];
        const double right = target.edge[t + 1U];
        while (source_cell + 1U < source.edge.size() && source.edge[source_cell + 1U] <= left) ++source_cell;
        std::size_t s = source_cell;
        double integral = 0.0;
        while (s < source.cell_count() && source.edge[s] < right) {
            const double overlap_left = std::max(left, source.edge[s]);
            const double overlap_right = std::min(right, source.edge[s + 1U]);
            if (overlap_right > overlap_left) integral += source_average[s] * (overlap_right - overlap_left);
            if (source.edge[s + 1U] >= right) break;
            ++s;
        }
        target_average[t] = integral / (right - left);
    }
    return target_average;
}

std::vector<double> conservative_face_flux_to_cell_rate(const cfd::fvm::PolyMesh& mesh,
                                                         std::span<const double> integrated_face_flux) {
    if (integrated_face_flux.size() != mesh.face_count()) throw std::invalid_argument("face flux/mesh size mismatch");
    std::vector<double> rate(mesh.cell_count(), 0.0);
    const auto& cells = mesh.cells();
    const auto& faces = mesh.faces();
    for (std::size_t face_index = 0U; face_index < faces.size(); ++face_index) {
        const auto& face = faces[face_index];
        const double flux = integrated_face_flux[face_index];
        rate[face.owner] += flux / cells[face.owner].volume;
        if (!face.boundary()) rate[face.neighbour] -= flux / cells[face.neighbour].volume;
    }
    return rate;
}

} // namespace cfd::multiphysics
