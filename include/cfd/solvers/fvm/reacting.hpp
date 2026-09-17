#pragma once
#include "cfd/chemistry/kinetics.hpp"
#include <span>
#include <vector>
namespace cfd::fvm {
[[nodiscard]] std::vector<std::vector<double>> reacting_cell_sources(const cfd::chemistry::ReactionNetwork& network,std::span<const std::vector<double>> cell_concentrations,std::span<const double> temperature);
}
