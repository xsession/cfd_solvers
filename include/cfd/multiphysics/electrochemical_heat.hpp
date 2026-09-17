#pragma once
#include <span>
#include <vector>
namespace cfd::multiphysics {
[[nodiscard]] std::vector<double> electrochemical_heat_source(std::span<const double> current_density,std::span<const double> overpotential,std::span<const double> ohmic_heating={});
}
