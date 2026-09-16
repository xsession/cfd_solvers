#pragma once

#include "cfd/solvers/fem/scalar_diffusion2d.hpp"

#include <functional>
#include <vector>

namespace cfd::fem {

struct MagneticFluxDensity2 {
    double x{};
    double y{};
};

// 2-D out-of-plane magnetic vector potential A_z:
//   -div(nu grad A_z) = J_z, nu = 1/mu.
// B = curl(A_z e_z) = (dA_z/dy, -dA_z/dx).
class Magnetostatics2D {
public:
    explicit Magnetostatics2D(Mesh2D mesh, ScalarDiffusion2DConfig config = {});
    void set_boundary(int patch, ScalarBoundaryCondition2D condition);
    void solve(const std::function<double(Node2)>& reluctivity,
               const std::function<double(Node2)>& current_density_z);

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return solver_.mesh(); }
    [[nodiscard]] const std::vector<double>& vector_potential() const noexcept { return solver_.solution(); }
    [[nodiscard]] std::vector<MagneticFluxDensity2> element_flux_density() const;
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return solver_.linear_result(); }

private:
    ScalarDiffusion2D solver_;
};

} // namespace cfd::fem
