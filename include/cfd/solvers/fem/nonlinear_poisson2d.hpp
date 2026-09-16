#pragma once

#include "cfd/core/newton.hpp"
#include "cfd/fem/mesh2d.hpp"

#include <functional>
#include <vector>

namespace cfd::fem {

struct NonlinearPoisson2DConfig {
    double diffusivity{1.0};
    double cubic_coefficient{1.0};
    cfd::core::NewtonConfig newton{};
};

// P1 Galerkin solver for
//   -div(k grad u) + beta*u^3 = f
// with Dirichlet data on all boundary nodes. This deliberately small nonlinear
// problem exercises the shared Newton/CSR/ILU-GMRES path used by later material
// and multiphysics nonlinearities.
class NonlinearPoisson2D {
public:
    explicit NonlinearPoisson2D(Mesh2D mesh, NonlinearPoisson2DConfig config = {});

    void solve(const std::function<double(Node2)>& source,
               const std::function<double(Node2)>& dirichlet_value);

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& solution() const noexcept { return solution_; }
    [[nodiscard]] const cfd::core::NewtonResult& nonlinear_result() const noexcept { return nonlinear_result_; }

private:
    Mesh2D mesh_;
    NonlinearPoisson2DConfig config_;
    std::vector<double> solution_;
    cfd::core::NewtonResult nonlinear_result_{};
};

} // namespace cfd::fem
