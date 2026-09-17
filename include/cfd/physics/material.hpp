#pragma once
#include <array>
#include <string>
namespace cfd::physics {
struct Material {std::string name;double density{};double heat_capacity{};double thermal_conductivity{};double electrical_conductivity{};std::array<double,3> relative_permittivity{1,1,1};std::array<double,3> relative_permeability{1,1,1};double refractive_index{1.0};};
void validate_material(const Material& material);
}
