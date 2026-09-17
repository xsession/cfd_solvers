#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/fvm/temporal.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <span>
#include <string_view>
#include <vector>

namespace cfd::fvm {

enum class TurbulenceScalarBoundaryType { zeroGradient, fixedValue };

struct TurbulenceScalarBoundaryCondition {
    TurbulenceScalarBoundaryType type{TurbulenceScalarBoundaryType::zeroGradient};
    double value{};
};

struct RansTransportControls {
    double dt{1.0e-4};
    std::size_t linear_iterations{400};
    std::size_t gmres_restart{30};
    double linear_tolerance{1.0e-10};
    TemporalScheme temporal_scheme{TemporalScheme::euler};
    double crank_nicolson_off_centering{1.0};
};

struct TwoEquationLinearResult {
    cfd::core::IterativeSolverResult first{};
    cfd::core::IterativeSolverResult second{};
};

[[nodiscard]] std::vector<double> turbulence_strain_rate_magnitude(
    const PolyMesh& mesh,
    std::span<const Vec3> velocity);

struct SpalartAllmarasConfig {
    RansTransportControls transport{};
    double molecular_viscosity{1.5e-5};
    double cb1{0.1355};
    double cb2{0.622};
    double cw2{0.3};
    double cw3{2.0};
    double cv1{7.1};
    double cs{0.3};
    double sigma_nu_tilde{2.0 / 3.0};
    double kappa{0.41};
    double minimum_nu_tilde{1.0e-12};
    double minimum_wall_distance{1.0e-9};
};

class SpalartAllmarasTransport {
public:
    explicit SpalartAllmarasTransport(PolyMesh mesh, SpalartAllmarasConfig config = {});
    void initialize(double nu_tilde);
    void initialize(std::span<const double> nu_tilde);
    void set_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_face_flux(std::vector<double> volumetric_flux);
    void set_wall_distance(double distance);
    void set_wall_distance(std::span<const double> distance);
    void set_grid_scale(double scale);
    void set_grid_scale(std::span<const double> scale);
    void set_strain_rate(double strain_rate);
    void set_strain_rate(std::span<const double> strain_rate);
    void set_velocity(std::span<const Vec3> velocity);
    cfd::core::IterativeSolverResult step();
    cfd::core::IterativeSolverResult run(std::size_t steps);
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& nu_tilde() const noexcept { return nu_tilde_; }
    [[nodiscard]] const std::vector<double>& strain_rate() const noexcept { return strain_rate_; }
    [[nodiscard]] const std::vector<double>& wall_distance() const noexcept { return wall_distance_; }
    [[nodiscard]] std::vector<double> kinematic_eddy_viscosity() const;
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return linear_result_; }

private:
    PolyMesh mesh_;
    SpalartAllmarasConfig config_;
    std::vector<TurbulenceScalarBoundaryCondition> boundary_;
    std::vector<double> nu_tilde_, previous_nu_tilde_, face_flux_, wall_distance_, strain_rate_;
    std::size_t steps_{};
    double time_{};
    cfd::core::IterativeSolverResult linear_result_{};
};

struct KEpsilonConfig {
    RansTransportControls transport{};
    double density{1.0};
    double molecular_viscosity{1.5e-5};
    double c_mu{0.09};
    double c1{1.44};
    double c2{1.92};
    double sigma_k{1.0};
    double sigma_epsilon{1.3};
    double production_limiter{20.0};
    double minimum_k{1.0e-12};
    double minimum_epsilon{1.0e-12};
    double maximum_eddy_viscosity_ratio{1.0e6};
};

class KEpsilonTransport {
public:
    explicit KEpsilonTransport(PolyMesh mesh, KEpsilonConfig config = {});
    void initialize(double k, double epsilon);
    void initialize(std::span<const double> k, std::span<const double> epsilon);
    void set_k_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_epsilon_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_face_flux(std::vector<double> volumetric_flux);
    void set_strain_rate(double strain_rate);
    void set_strain_rate(std::span<const double> strain_rate);
    void set_velocity(std::span<const Vec3> velocity);
    TwoEquationLinearResult step();
    TwoEquationLinearResult run(std::size_t steps);
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& k() const noexcept { return k_; }
    [[nodiscard]] const std::vector<double>& epsilon() const noexcept { return epsilon_; }
    [[nodiscard]] std::vector<double> kinematic_eddy_viscosity() const;
    [[nodiscard]] std::vector<double> dynamic_eddy_viscosity() const;
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] const TwoEquationLinearResult& linear_result() const noexcept { return linear_result_; }

private:
    PolyMesh mesh_;
    KEpsilonConfig config_;
    std::vector<TurbulenceScalarBoundaryCondition> k_boundary_, epsilon_boundary_;
    std::vector<double> k_, epsilon_, previous_k_, previous_epsilon_, face_flux_, strain_rate_;
    std::size_t steps_{};
    double time_{};
    TwoEquationLinearResult linear_result_{};
};

enum class SSTDESZonalFilter { none, f1, f2 };

struct KOmegaSSTConfig {
    RansTransportControls transport{};
    double density{1.0};
    double molecular_viscosity{1.5e-5};
    double beta_star{0.09};
    double beta1{0.075};
    double beta2{0.0828};
    double gamma1{5.0 / 9.0};
    double gamma2{0.44};
    double alpha_k1{0.85};
    double alpha_k2{1.0};
    double alpha_omega1{0.5};
    double alpha_omega2{0.856};
    double a1{0.31};
    double production_limiter{10.0};
    double minimum_k{1.0e-12};
    double minimum_omega{1.0e-12};
    double minimum_wall_distance{1.0e-9};
    double maximum_eddy_viscosity_ratio{1.0e6};
    bool des_enabled{false};
    double c_des{0.61};
    SSTDESZonalFilter des_zonal_filter{SSTDESZonalFilter::f2};
};

class KOmegaSSTTransport {
public:
    explicit KOmegaSSTTransport(PolyMesh mesh, KOmegaSSTConfig config = {});
    void initialize(double k, double omega);
    void initialize(std::span<const double> k, std::span<const double> omega);
    void set_k_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_omega_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_face_flux(std::vector<double> volumetric_flux);
    void set_wall_distance(double distance);
    void set_wall_distance(std::span<const double> distance);
    void set_grid_scale(double scale);
    void set_grid_scale(std::span<const double> scale);
    void set_strain_rate(double strain_rate);
    void set_strain_rate(std::span<const double> strain_rate);
    void set_velocity(std::span<const Vec3> velocity);
    TwoEquationLinearResult step();
    TwoEquationLinearResult run(std::size_t steps);
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& k() const noexcept { return k_; }
    [[nodiscard]] const std::vector<double>& omega() const noexcept { return omega_; }
    [[nodiscard]] const std::vector<double>& blending_f1() const noexcept { return blending_f1_; }
    [[nodiscard]] const std::vector<double>& blending_f2() const noexcept { return blending_f2_; }
    [[nodiscard]] const std::vector<double>& hybrid_dissipation_factor() const noexcept { return hybrid_factor_; }
    [[nodiscard]] std::vector<double> kinematic_eddy_viscosity() const;
    [[nodiscard]] std::vector<double> dynamic_eddy_viscosity() const;
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] const TwoEquationLinearResult& linear_result() const noexcept { return linear_result_; }

private:
    PolyMesh mesh_;
    KOmegaSSTConfig config_;
    std::vector<TurbulenceScalarBoundaryCondition> k_boundary_, omega_boundary_;
    std::vector<double> k_, omega_, previous_k_, previous_omega_, face_flux_, wall_distance_, grid_scale_, strain_rate_;
    std::vector<double> blending_f1_, blending_f2_, hybrid_factor_;
    std::size_t steps_{};
    double time_{};
    TwoEquationLinearResult linear_result_{};
    void update_blending();
};

struct SymmetricTensor3 {
    double xx{}, yy{}, zz{}, xy{}, xz{}, yz{};
};

using VelocityGradient3 = std::array<double, 9>; // row-major: dU_i/dx_j

[[nodiscard]] double turbulent_kinetic_energy(const SymmetricTensor3& stress) noexcept;
[[nodiscard]] bool reynolds_stress_is_realizable(const SymmetricTensor3& stress, double tolerance = 1.0e-12) noexcept;

struct ReynoldsStressConfig {
    RansTransportControls transport{};
    double density{1.0};
    double molecular_viscosity{1.5e-5};
    double c_mu{0.09};
    double pressure_strain_slow{1.8};
    double pressure_strain_rapid{0.6};
    double c_epsilon1{1.44};
    double c_epsilon2{1.92};
    double sigma_r{0.82};
    double sigma_epsilon{1.3};
    double minimum_k{1.0e-12};
    double minimum_epsilon{1.0e-12};
    double minimum_normal_stress{1.0e-14};
    double maximum_eddy_viscosity_ratio{1.0e6};
};

struct ReynoldsStressLinearResult {
    std::array<cfd::core::IterativeSolverResult, 6> stress{};
    cfd::core::IterativeSolverResult epsilon{};
};

class ReynoldsStressTransport {
public:
    explicit ReynoldsStressTransport(PolyMesh mesh, ReynoldsStressConfig config = {});
    void initialize_isotropic(double k, double epsilon);
    void initialize(std::span<const SymmetricTensor3> stress, std::span<const double> epsilon);
    void set_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, SymmetricTensor3 value = {});
    void set_epsilon_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_face_flux(std::vector<double> volumetric_flux);
    void set_velocity(std::span<const Vec3> velocity);
    void set_velocity_gradient(std::span<const VelocityGradient3> gradient);
    void set_additional_source(std::function<SymmetricTensor3(Vec3, double)> source);
    ReynoldsStressLinearResult step();
    ReynoldsStressLinearResult run(std::size_t steps);
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] std::vector<SymmetricTensor3> reynolds_stress() const;
    [[nodiscard]] std::vector<double> turbulent_kinetic_energy() const;
    [[nodiscard]] const std::vector<double>& epsilon() const noexcept { return epsilon_; }
    [[nodiscard]] std::vector<double> kinematic_eddy_viscosity() const;
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] const ReynoldsStressLinearResult& linear_result() const noexcept { return linear_result_; }

private:
    PolyMesh mesh_;
    ReynoldsStressConfig config_;
    std::array<std::vector<TurbulenceScalarBoundaryCondition>, 6> stress_boundary_;
    std::vector<TurbulenceScalarBoundaryCondition> epsilon_boundary_;
    std::array<std::vector<double>, 6> stress_, previous_stress_;
    std::vector<double> epsilon_, previous_epsilon_, face_flux_;
    std::vector<VelocityGradient3> velocity_gradient_;
    std::function<SymmetricTensor3(Vec3, double)> additional_source_;
    std::size_t steps_{};
    double time_{};
    ReynoldsStressLinearResult linear_result_{};
    void enforce_realizability();
};

} // namespace cfd::fvm
