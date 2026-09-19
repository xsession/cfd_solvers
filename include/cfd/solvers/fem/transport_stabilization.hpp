#pragma once

namespace cfd::fem {

struct SupgParameters {
    double tau_s{};
    double element_peclet{};
};

// Element-level SUPG time scale for a scalar advection-diffusion equation.
[[nodiscard]] SupgParameters supg_parameters(double element_length_m, double advection_speed_m_per_s,
                                             double diffusivity_m2_per_s, double timestep_s);

// Face-local DG upwind plus symmetric penalty numerical flux.
[[nodiscard]] double dg_upwind_flux(double left_value, double right_value, double outward_speed, double penalty);

} // namespace cfd::fem
