#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fem/mesh2d.hpp"
#include "cfd/solvers/fem/elasticity2d.hpp"

#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

namespace cfd::fem {

struct AxisymmetricElasticityConfig {
    double youngs_modulus{1.0};
    double poisson_ratio{0.3};
    std::size_t max_iterations{5000U};
    double relative_tolerance{1.0e-10};
};

// Axisymmetric small-strain isotropic elasticity. Mesh x is radial coordinate r
// and mesh y is axial coordinate z. The weak form includes the 2*pi*r measure
// and the hoop strain u_r/r.
class AxisymmetricElasticity {
public:
    explicit AxisymmetricElasticity(Mesh2D mesh, AxisymmetricElasticityConfig config = {});
    void set_dirichlet(const std::function<bool(Node2)>& predicate,
                       const std::function<Displacement2(Node2)>& value,
                       bool fix_radial = true,
                       bool fix_axial = true);
    void solve(const std::function<Displacement2(Node2)>& body_force = {});

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<Displacement2>& displacement() const noexcept { return displacement_; }
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return linear_result_; }

private:
    Mesh2D mesh_;
    AxisymmetricElasticityConfig config_;
    std::vector<Displacement2> displacement_;
    std::vector<double> prescribed_r_;
    std::vector<double> prescribed_z_;
    cfd::core::IterativeSolverResult linear_result_{};
};

} // namespace cfd::fem
