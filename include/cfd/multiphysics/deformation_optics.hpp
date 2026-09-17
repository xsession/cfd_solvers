#pragma once
#include "cfd/solvers/optics/sequential.hpp"
#include <functional>
namespace cfd::multiphysics {
// Transfers axial structural displacement onto sequential surface vertices.
[[nodiscard]] cfd::optics::SequentialOpticalSystem deform_optical_system_axially(
    const cfd::optics::SequentialOpticalSystem& system,
    const std::function<double(double)>& axial_displacement);
} // namespace cfd::multiphysics
