#pragma once

#include <cstddef>
#include <vector>

namespace cfd::fdtd {

struct CylindricalGridLevel {
    std::size_t radial_cells{};
    std::size_t axial_cells{};
    double radial_spacing{};
    double axial_spacing{};
};

[[nodiscard]] std::vector<CylindricalGridLevel> make_cylindrical_multigrid_levels(std::size_t radial_cells,
                                                                                  std::size_t axial_cells,
                                                                                  double radius_m, double height_m,
                                                                                  std::size_t levels);

} // namespace cfd::fdtd
