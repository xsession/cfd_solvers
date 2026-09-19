#pragma once

#include "cfd/fem/mesh2d.hpp"

#include <array>

namespace cfd::fem {

struct MixedStokesTri3LocalSystem {
    std::array<double, 81> matrix{}; // [u_x0,u_y0,u_x1,u_y1,u_x2,u_y2,p0,p1,p2]
    double area{};
};

// Equal-order P1/P1 mixed Stokes element with pressure-mass stabilization.
// This is a compact local block assembly baseline intended for didactic
// monolithic coupling; callers impose global constraints and boundary data.
[[nodiscard]] MixedStokesTri3LocalSystem assemble_mixed_stokes_tri3(const Mesh2D& mesh, const Tri3& triangle,
                                                                    double viscosity,
                                                                    double pressure_stabilization = 1.0e-8);

} // namespace cfd::fem
