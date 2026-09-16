#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fem/mesh2d.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::fem {

struct Heat2DConfig {
    double conductivity{1.0};
    double volumetric_heat_capacity{1.0};
    double dt{1.0e-3};
    std::size_t max_iterations{2000U};
    double relative_tolerance{1.0e-10};
};

class Heat2D {
public:
    explicit Heat2D(Mesh2D mesh, Heat2DConfig config = {});

    void initialize(const std::function<double(Node2)>& temperature);
    void set_dirichlet(const std::function<double(Node2,double)>& boundary_temperature);
    void step(const std::function<double(Node2,double)>& volumetric_source = {});
    void step_element_source(std::span<const double> element_volumetric_source);
    void run(std::size_t steps, const std::function<double(Node2,double)>& volumetric_source = {});

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& temperature() const noexcept { return temperature_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return linear_result_; }

private:
    Mesh2D mesh_;
    Heat2DConfig config_;
    std::vector<double> temperature_;
    std::function<double(Node2,double)> boundary_temperature_;
    cfd::core::IterativeSolverResult linear_result_{};
    double time_{};

    using ElementSource = std::function<double(std::size_t,Node2,double)>;
    void step_impl(const ElementSource& source);
};

} // namespace cfd::fem
