#pragma once

#include "cfd/solvers/fem/scalar_diffusion2d.hpp"

#include <functional>
#include <vector>

namespace cfd::fem {

struct DarcyVelocity2 {
    double x{};
    double y{};
};

struct Darcy2DConfig {
    double viscosity{1.0e-3};
    double density{1000.0};
    DarcyVelocity2 gravity{0.0,0.0};
    ScalarDiffusion2DConfig linear{};
};

// Saturated isotropic Darcy flow:
//   v = -(K/mu) (grad p - rho g)
//   -div((K/mu) grad p) = q - div((K/mu) rho g)
// The baseline assumes constant gravity and scalar K, so the gravity-divergence
// term vanishes for constant K; gravity still appears in reconstructed velocity.
class Darcy2D {
public:
    explicit Darcy2D(Mesh2D mesh, Darcy2DConfig config = {});

    void set_boundary(int patch, ScalarBoundaryCondition2D condition);
    void solve(const std::function<double(Node2)>& permeability,
               const std::function<double(Node2)>& volumetric_source = {});

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return solver_.mesh(); }
    [[nodiscard]] const std::vector<double>& pressure() const noexcept { return solver_.solution(); }
    [[nodiscard]] std::vector<DarcyVelocity2> element_velocity(const std::function<double(Node2)>& permeability) const;
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return solver_.linear_result(); }

private:
    Darcy2DConfig config_;
    ScalarDiffusion2D solver_;
};

} // namespace cfd::fem
