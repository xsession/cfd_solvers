#pragma once

#include "cfd/solvers/fem/scalar_diffusion2d.hpp"

#include <functional>
#include <span>
#include <vector>

namespace cfd::fem {

struct MagneticFluxDensity2 {
    double x{};
    double y{};
};

// 2-D out-of-plane magnetic vector potential A_z:
//   -div(nu grad A_z) = J_z, nu = 1/mu.
// B = curl(A_z e_z) = (dA_z/dy, -dA_z/dx).
class Magnetostatics2D {
public:
    explicit Magnetostatics2D(Mesh2D mesh, ScalarDiffusion2DConfig config = {});
    void set_boundary(int patch, ScalarBoundaryCondition2D condition);
    void solve(const std::function<double(Node2)>& reluctivity,
               const std::function<double(Node2)>& current_density_z);

    struct NonlinearResult { std::size_t iterations{}; double maximum_relative_reluctivity_change{}; bool converged{}; };
    // Elementwise Picard nonlinear B-H integration. The law receives |B| [T]
    // and returns reluctivity nu=H/B [1/H/m].
    [[nodiscard]] NonlinearResult solve_nonlinear(
        const std::function<double(double flux_density_t)>& reluctivity_law,
        const std::function<double(Node2)>& current_density_z,
        std::size_t max_iterations=50U,double relative_tolerance=1.0e-6,double relaxation=0.7);

    [[nodiscard]] const Mesh2D& mesh() const noexcept { return solver_.mesh(); }
    [[nodiscard]] const std::vector<double>& vector_potential() const noexcept { return solver_.solution(); }
    [[nodiscard]] std::vector<MagneticFluxDensity2> element_flux_density() const;
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return solver_.linear_result(); }

private:
    ScalarDiffusion2D solver_;
};

class PiecewiseLinearBHCurve {
public:
    // Samples are B [T], H [A/m], strictly increasing from the origin.
    PiecewiseLinearBHCurve(std::vector<double> flux_density_t,std::vector<double> field_strength_a_per_m);
    [[nodiscard]] double field_strength(double flux_density_t) const;
    [[nodiscard]] double reluctivity(double flux_density_t) const;
private:
    std::vector<double> b_,h_;
};

struct RectangularStrandedCoil2D {
    double xmin{},xmax{},ymin{},ymax{};
    double turns{1.0};
    double current_a{1.0};
    [[nodiscard]] double area_m2() const;
    [[nodiscard]] bool contains(Node2 point) const noexcept;
    [[nodiscard]] double current_density_a_per_m2(Node2 point) const;
};
[[nodiscard]] double coil_flux_linkage_wb_turn(
    const Mesh2D& mesh,std::span<const double> vector_potential_a_per_m,
    const RectangularStrandedCoil2D& coil);
[[nodiscard]] double coil_inductance_h(
    const Mesh2D& mesh,std::span<const double> vector_potential_a_per_m,
    const RectangularStrandedCoil2D& coil);

struct MaxwellForceTorque2D {
    double force_x_n_per_m{};
    double force_y_n_per_m{};
    double torque_z_n{}; // per unit out-of-plane depth
};
// Maxwell-stress integral over one boundary patch, assuming the supplied B
// samples are constant per adjacent triangle and permeability is uniform.
[[nodiscard]] MaxwellForceTorque2D maxwell_stress_boundary_force_2d(
    const Mesh2D& mesh,std::span<const MagneticFluxDensity2> element_flux_density,
    int patch,double permeability_h_per_m,Node2 torque_origin={});

} // namespace cfd::fem
