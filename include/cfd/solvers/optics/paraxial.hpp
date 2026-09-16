#pragma once

#include "cfd/solvers/optics/sequential.hpp"

#include <cstddef>
#include <vector>

namespace cfd::optics {

struct ParaxialRay {
    double height{};
    double angle{}; // radians, paraxial slope angle
    double z{};
};

struct ParaxialTraceResult {
    ParaxialRay ray{};
    bool valid{};
    std::size_t surfaces_traced{};
};

[[nodiscard]] ParaxialTraceResult trace_paraxial(const SequentialOpticalSystem& system,
                                                  ParaxialRay ray,
                                                  double object_space_index=1.0);

// Back focal distance from the last surface for an incoming on-axis parallel
// ray. Returns positive distance when a real paraxial focus lies downstream.
[[nodiscard]] double paraxial_back_focal_distance(const SequentialOpticalSystem& system,
                                                   double object_space_index=1.0);

} // namespace cfd::optics
