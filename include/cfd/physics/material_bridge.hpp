#pragma once
#include "cfd/physics/material.hpp"
#include "cfd/solvers/fdtd/maxwell3d.hpp"
namespace cfd::physics {
[[nodiscard]] cfd::fdtd::Maxwell3DConfig apply_to_fdtd(cfd::fdtd::Maxwell3DConfig config,const Material& material);
[[nodiscard]] double fem_conductivity(const Material& material);
[[nodiscard]] double optical_refractive_index(const Material& material);
} // namespace cfd::physics
