#pragma once
#include "cfd/chemistry/kinetics.hpp"
#include <span>
#include <vector>
namespace cfd::chemistry {
struct ReactionPathEdge {std::size_t from{},to{};double molar_flux{};};
[[nodiscard]] std::vector<ReactionPathEdge> reaction_path(const ReactionNetwork& network,std::span<const double> concentrations,double temperature,double minimum_flux=0.0);
}
