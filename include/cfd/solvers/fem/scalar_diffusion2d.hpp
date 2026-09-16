#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fem/mesh2d.hpp"

#include <cstddef>
#include <functional>
#include <map>
#include <vector>

namespace cfd::fem {

enum class ScalarBoundaryType {
    dirichlet,
    neumann,
    robin
};

struct ScalarBoundaryCondition2D {
    ScalarBoundaryType type{ScalarBoundaryType::dirichlet};
    // Dirichlet: prescribed value.
    // Neumann: outward k*grad(u).n value.
    // Robin: alpha*u + k*grad(u).n = beta, where value returns beta.
    std::function<double(Node2)> value;
    std::function<double(Node2)> alpha;
};

struct ScalarDiffusion2DConfig {
    std::size_t max_iterations{4000U};
    double relative_tolerance{1.0e-11};
};

// General P1 scalar FEM for
//   -div(k grad u) + c u = f
// on a triangle mesh, with per-patch Dirichlet/Neumann/Robin conditions.
class ScalarDiffusion2D {
public:
    explicit ScalarDiffusion2D(Mesh2D mesh, ScalarDiffusion2DConfig config = {});

    void set_boundary(int patch, ScalarBoundaryCondition2D condition);
    void solve(const std::function<double(Node2)>& diffusivity,
               const std::function<double(Node2)>& reaction,
               const std::function<double(Node2)>& source);

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& solution() const noexcept { return solution_; }
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return linear_result_; }

private:
    Mesh2D mesh_;
    ScalarDiffusion2DConfig config_;
    std::map<int, ScalarBoundaryCondition2D> boundary_;
    std::vector<double> solution_;
    cfd::core::IterativeSolverResult linear_result_{};
};

} // namespace cfd::fem
