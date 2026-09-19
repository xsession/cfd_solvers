#include "cfd/solvers/fem/thermoelastic_block.hpp"

#include <stdexcept>

namespace cfd::fem {

ThermoelasticBlockLocalSystem assemble_thermoelastic_block(double area, double stiffness, double conductivity,
                                                           double capacity, double dt, double expansion) {
    if (!(area > 0.0) || !(stiffness > 0.0) || !(conductivity >= 0.0) || !(capacity >= 0.0) || !(dt > 0.0)) {
        throw std::invalid_argument("invalid thermoelastic block controls");
    }
    ThermoelasticBlockLocalSystem result;
    result.matrix[0] = stiffness * area;
    result.matrix[4] = stiffness * area;
    result.matrix[8] = (conductivity + capacity / dt) * area;
    result.matrix[2] = result.matrix[5] = -expansion * stiffness * area;
    result.matrix[6] = result.matrix[7] = result.matrix[2];
    return result;
}

} // namespace cfd::fem
