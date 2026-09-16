#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fem/mesh2d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace cfd::fem {

struct Poisson2DConfig {
    double diffusivity{1.0};
    std::size_t max_iterations{2000U};
    double relative_tolerance{1.0e-10};
};

class Poisson2D {
public:
    explicit Poisson2D(Mesh2D mesh, Poisson2DConfig config = {});

    void solve(const std::function<double(Node2)>& source,
               const std::function<double(Node2)>& dirichlet_value);

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& solution() const noexcept { return solution_; }
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return linear_result_; }

private:
    Mesh2D mesh_;
    Poisson2DConfig config_;
    std::vector<double> solution_;
    cfd::core::IterativeSolverResult linear_result_{};
};

} // namespace cfd::fem
