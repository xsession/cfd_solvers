#pragma once

#include "cfd/solvers/fem/scalar_diffusion2d.hpp"

#include <functional>
#include <vector>

namespace cfd::fem {

struct Vector2 {
    double x{};
    double y{};
};

class Electrostatics2D {
public:
    explicit Electrostatics2D(Mesh2D mesh, ScalarDiffusion2DConfig config = {});
    void set_boundary(int patch, ScalarBoundaryCondition2D condition);
    // Solves -div(epsilon_r grad(phi)) = rho/epsilon0. Supplying the source in
    // rho/epsilon0 units keeps the linear system well-scaled in SI workflows.
    void solve(const std::function<double(Node2)>& relative_permittivity,
               const std::function<double(Node2)>& charge_over_epsilon0);

    [[nodiscard]] const std::vector<double>& potential() const noexcept { return solver_.solution(); }
    [[nodiscard]] const Mesh2D& mesh() const noexcept { return solver_.mesh(); }
    [[nodiscard]] std::vector<Vector2> element_electric_field() const;
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return solver_.linear_result(); }

private:
    ScalarDiffusion2D solver_;
};

class DCConduction2D {
public:
    explicit DCConduction2D(Mesh2D mesh, ScalarDiffusion2DConfig config = {});
    void set_boundary(int patch, ScalarBoundaryCondition2D condition);
    // Solves -div(sigma grad(phi)) = source.
    void solve(const std::function<double(Node2)>& conductivity,
               const std::function<double(Node2)>& current_source = {});

    [[nodiscard]] const std::vector<double>& potential() const noexcept { return solver_.solution(); }
    [[nodiscard]] const Mesh2D& mesh() const noexcept { return solver_.mesh(); }
    [[nodiscard]] std::vector<Vector2> element_current_density(
        const std::function<double(Node2)>& conductivity) const;
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return solver_.linear_result(); }

private:
    ScalarDiffusion2D solver_;
};

} // namespace cfd::fem
