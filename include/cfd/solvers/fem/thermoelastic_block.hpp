#pragma once

#include <array>

namespace cfd::fem {

struct ThermoelasticBlockLocalSystem {
    std::array<double, 9> matrix{}; // [u_x,u_y,T] block
};

[[nodiscard]] ThermoelasticBlockLocalSystem assemble_thermoelastic_block(double area, double mechanical_stiffness,
                                                                         double thermal_conductivity,
                                                                         double thermal_capacity, double time_step,
                                                                         double thermal_expansion);

} // namespace cfd::fem
