#include "cfd/physics/material.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::physics {void validate_material(const Material&m){if(m.name.empty()||m.density<0||m.heat_capacity<0||m.thermal_conductivity<0||m.electrical_conductivity<0||!(m.refractive_index>0))throw std::invalid_argument("invalid shared material");for(double v:m.relative_permittivity)if(!(v>0)||!std::isfinite(v))throw std::invalid_argument("invalid permittivity");for(double v:m.relative_permeability)if(!(v>0)||!std::isfinite(v))throw std::invalid_argument("invalid permeability");}}
