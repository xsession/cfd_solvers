#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fem/mesh2d.hpp"

#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

namespace cfd::fem {

enum class ElasticityMode2D { planeStress, planeStrain };

struct Displacement2 {
    double x{};
    double y{};
};

struct Elasticity2DConfig {
    double youngs_modulus{1.0};
    double poisson_ratio{0.3};
    ElasticityMode2D mode{ElasticityMode2D::planeStress};
    std::size_t max_iterations{4000U};
    double relative_tolerance{1.0e-10};
};

class Elasticity2D {
public:
    explicit Elasticity2D(Mesh2D mesh, Elasticity2DConfig config = {});

    // Constrain selected components at nodes matching predicate.
    void set_dirichlet(const std::function<bool(Node2)>& predicate,
                       const std::function<Displacement2(Node2)>& value,
                       bool fix_x = true,
                       bool fix_y = true);
    void solve(const std::function<Displacement2(Node2)>& body_force = {});
    // Small-strain isotropic thermal expansion on the same nodal Tri3 mesh.
    void solve_thermal(std::span<const double> nodal_temperature,double expansion_coefficient,
                       double reference_temperature,const std::function<Displacement2(Node2)>& body_force = {});

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<Displacement2>& displacement() const noexcept { return displacement_; }
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return linear_result_; }

private:
    Mesh2D mesh_;
    Elasticity2DConfig config_;
    std::vector<Displacement2> displacement_;
    std::vector<double> prescribed_x_;
    std::vector<double> prescribed_y_;
    cfd::core::IterativeSolverResult linear_result_{};
    void solve_impl(const std::function<Displacement2(Node2)>& body_force,
                    std::span<const double> temperature,double expansion,double reference_temperature);
};

} // namespace cfd::fem
