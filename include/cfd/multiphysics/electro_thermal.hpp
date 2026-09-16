#pragma once

#include "cfd/solvers/fem/electromagnetics2d.hpp"
#include "cfd/solvers/fem/heat2d.hpp"

#include <functional>
#include <vector>

namespace cfd::multiphysics {

// Explicit one-way electro-thermal coupling on a shared Tri3 mesh:
//   DC conduction -> element current density -> Joule heat -> transient heat FEM.
// The two physics solvers retain independent state, while the coupling term is
// transferred conservatively as an element-wise volumetric source.
class JouleHeatingCoupler2D {
public:
    JouleHeatingCoupler2D(cfd::fem::Mesh2D mesh,
                          cfd::fem::Heat2DConfig thermal_config = {},
                          cfd::fem::ScalarDiffusion2DConfig electrical_config = {});

    void set_electrical_boundary(int patch,cfd::fem::ScalarBoundaryCondition2D condition);
    void initialize_temperature(const std::function<double(cfd::fem::Node2)>& temperature);
    void set_thermal_dirichlet(const std::function<double(cfd::fem::Node2,double)>& boundary_temperature);

    void solve_electrical(const std::function<double(cfd::fem::Node2)>& conductivity,
                          const std::function<double(cfd::fem::Node2)>& current_source = {});
    void thermal_step();
    void thermal_run(std::size_t steps);

    [[nodiscard]] const cfd::fem::DCConduction2D& electrical() const noexcept { return electrical_; }
    [[nodiscard]] const cfd::fem::Heat2D& thermal() const noexcept { return thermal_; }
    [[nodiscard]] const std::vector<double>& joule_heating_density() const noexcept { return joule_heating_density_; }

private:
    cfd::fem::DCConduction2D electrical_;
    cfd::fem::Heat2D thermal_;
    std::vector<double> joule_heating_density_;
    bool electrical_solved_{};
};

} // namespace cfd::multiphysics
