#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <span>
#include <vector>

namespace cfd::multiphysics {

struct CellGrid1D {
    std::vector<double> edge;
    void validate() const;
    [[nodiscard]] std::size_t cell_count() const noexcept { return edge.empty() ? 0U : edge.size() - 1U; }
};

// Conservative remap of cell averages using exact interval overlap. Source and
// target grids must cover the same 1-D physical interval.
[[nodiscard]] std::vector<double> conservative_cell_average_transfer(
    const CellGrid1D& source,
    std::span<const double> source_average,
    const CellGrid1D& target);

// Convert oriented, face-integrated fluxes to a cell-local volumetric rate.
// Internal-face contributions cancel exactly in the volume-weighted global sum.
[[nodiscard]] std::vector<double> conservative_face_flux_to_cell_rate(
    const cfd::fvm::PolyMesh& mesh,
    std::span<const double> integrated_face_flux);

} // namespace cfd::multiphysics
