#include "cfd/solvers/fdtd/cylindrical_multigrid.hpp"

#include <algorithm>
#include <stdexcept>

namespace cfd::fdtd {

std::vector<CylindricalGridLevel> make_cylindrical_multigrid_levels(std::size_t radial, std::size_t axial,
                                                                    double radius, double height, std::size_t levels) {
    if (radial < 2U || axial < 2U || !(radius > 0.0) || !(height > 0.0) || levels == 0U)
        throw std::invalid_argument("invalid cylindrical multigrid controls");
    std::vector<CylindricalGridLevel> result;
    result.reserve(levels);
    for (std::size_t level = 0; level < levels; ++level) {
        result.push_back({radial, axial, radius / static_cast<double>(radial), height / static_cast<double>(axial)});
        if (radial <= 2U || axial <= 2U)
            break;
        radial = std::max<std::size_t>(2U, radial / 2U);
        axial = std::max<std::size_t>(2U, axial / 2U);
    }
    return result;
}

} // namespace cfd::fdtd
