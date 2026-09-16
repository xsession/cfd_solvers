#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fem/mesh3d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace cfd::fem {

struct Poisson3DConfig {
    double diffusivity{1.0};
    std::size_t max_iterations{5000U};
    double relative_tolerance{1.0e-10};
};

class Poisson3D {
public:
    explicit Poisson3D(Mesh3D mesh,Poisson3DConfig config={});
    void solve(const std::function<double(Point3)>& source,
               const std::function<double(Point3)>& dirichlet_value);
    [[nodiscard]] const Mesh3D& mesh() const noexcept{return mesh_;}
    [[nodiscard]] const std::vector<double>& solution() const noexcept{return solution_;}
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept{return linear_result_;}
private:
    Mesh3D mesh_;
    Poisson3DConfig config_;
    std::vector<double> solution_;
    cfd::core::IterativeSolverResult linear_result_{};
};

} // namespace cfd::fem
