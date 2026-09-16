#pragma once

#include "cfd/core/conjugate_gradient.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fvm/pressure_velocity.hpp"
#include "cfd/fvm/schemes.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <string_view>
#include <vector>

namespace cfd::fvm {

struct CollocatedIncompressibleConfig {
    double density{1.0};
    double kinematic_viscosity{1.0e-2};
    double dt{1.0e-2};                 // physical dt for PISO/PIMPLE; pseudo-dt for SIMPLE
    std::size_t momentum_sweeps{3}; // legacy fixed-point fallback
    bool use_krylov_momentum{true};
    std::size_t momentum_iterations{300};
    double momentum_tolerance{1.0e-10};
    std::size_t pressure_iterations{1000};
    double pressure_tolerance{1.0e-10};
    std::size_t nonorthogonal_correctors{2};
    std::size_t pressure_correctors{2}; // PISO correctors
    std::size_t outer_correctors{2};    // PIMPLE outer loops
    double velocity_relaxation{0.7};
    double pressure_relaxation{0.3};
    FaceInterpolationScheme convection_scheme{FaceInterpolationScheme::upwind};
    bool include_convection{true};
};

struct CollocatedIterationInfo {
    double velocity_rms_change{};
    double continuity_l2{};
    double continuity_max{};
    cfd::core::ConjugateGradientResult pressure{};
    std::size_t pressure_solves{};
};

// Cell-centred finite-volume incompressible solver for arbitrary PolyMesh.
// The pressure-free momentum predictor HbyA and V/aP mobility are converted
// to conservative face fluxes with direct face pressure differences
// (Rhie-Chow style), followed by non-orthogonal pressure correction.
class CollocatedIncompressible {
public:
    CollocatedIncompressible(PolyMesh mesh, CollocatedIncompressibleConfig config = {});

    void set_velocity_boundary(std::string_view patch,
                               VelocityBoundaryType type,
                               Vec3 value = {});
    void set_pressure_boundary(std::string_view patch,
                               PressureBoundaryType type,
                               double value = 0.0);

    void initialize_uniform(Vec3 velocity = {}, double pressure = 0.0);
    void initialize_fields(const std::function<Vec3(Vec3)>& velocity,
                           const std::function<double(Vec3)>& pressure);

    // One pseudo-time SIMPLE iteration (one pressure solve sequence, relaxed).
    CollocatedIterationInfo iterate_simple();
    CollocatedIterationInfo solve_simple(std::size_t max_iterations,
                                         double velocity_tolerance = 1.0e-8,
                                         double continuity_tolerance = 1.0e-8);

    // One physical time step. PISO uses one momentum predictor with multiple
    // pressure correctors; PIMPLE repeats the momentum/pressure sequence.
    CollocatedIterationInfo step_piso();
    CollocatedIterationInfo step_pimple();

    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<Vec3>& velocity() const noexcept { return velocity_; }
    [[nodiscard]] const std::vector<double>& pressure() const noexcept { return pressure_; }
    [[nodiscard]] const std::vector<double>& face_flux() const noexcept { return face_flux_; }
    [[nodiscard]] const std::vector<Vec3>& h_by_a() const noexcept { return h_by_a_; }
    [[nodiscard]] const std::vector<double>& pressure_mobility() const noexcept { return pressure_mobility_; }
    [[nodiscard]] const CollocatedIncompressibleConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::array<cfd::core::IterativeSolverResult, 3>& momentum_results() const noexcept { return momentum_results_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }
    [[nodiscard]] double continuity_l2() const;
    [[nodiscard]] double continuity_max() const;
    [[nodiscard]] double kinetic_energy() const;

private:
    PolyMesh mesh_;
    CollocatedIncompressibleConfig config_;
    std::vector<VelocityBoundaryCondition> velocity_boundary_;
    std::vector<PressureBoundaryCondition> pressure_boundary_;
    std::vector<Vec3> velocity_;
    std::vector<Vec3> old_velocity_;
    std::vector<Vec3> h_by_a_;
    std::vector<double> pressure_mobility_;
    std::vector<double> pressure_;
    std::vector<double> face_flux_;
    std::vector<double> pressure_rhs_;
    std::vector<double> pressure_face_coefficient_;
    cfd::core::ConjugateGradientWorkspace pressure_workspace_;
    std::array<cfd::core::KrylovWorkspace, 3> momentum_workspace_;
    std::array<cfd::core::IterativeSolverResult, 3> momentum_results_{};
    cfd::core::ConjugateGradientResult pressure_result_{};
    double time_{};
    std::size_t steps_{};

    [[nodiscard]] bool has_fixed_pressure_boundary() const noexcept;
    void rebuild_flux_from_velocity();
    void momentum_predictor(std::span<const Vec3> time_source);
    cfd::core::ConjugateGradientResult correct_pressure(double pressure_relaxation);
    [[nodiscard]] std::vector<double> predicted_flux() const;
    void update_pressure_coefficients();
    [[nodiscard]] std::vector<double> nonorthogonal_pressure_flux(std::span<const double> p) const;
    void apply_pressure_operator(std::span<const double> x, std::span<double> out) const;
    [[nodiscard]] double velocity_rms_change(std::span<const Vec3> before) const;
    CollocatedIterationInfo coupled_sequence(std::span<const Vec3> time_source,
                                             std::size_t pressure_correctors,
                                             double pressure_relaxation);
};

} // namespace cfd::fvm
