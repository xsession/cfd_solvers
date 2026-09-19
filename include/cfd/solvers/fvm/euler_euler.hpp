#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::fvm {

// Explicit two-fluid baseline. The dispersed volume fraction is transported
// conservatively with its phase velocity; phase momenta use the matching
// upwind phase fluxes and exchange equal-and-opposite interphase drag.
struct EulerEulerConfig {
    double dt{1.0e-3};
    double primary_density{1000.0};
    double dispersed_density{1.0};
    double drag_coefficient{1.0};
    bool quadratic_drag{false};
    Vec3 primary_body_acceleration{};
    Vec3 dispersed_body_acceleration{};
    double minimum_phase_fraction{1.0e-12};
};

class EulerEulerTransport {
public:
    explicit EulerEulerTransport(PolyMesh mesh, EulerEulerConfig config = {});

    void initialize(double dispersed_fraction, Vec3 primary_velocity = {}, Vec3 dispersed_velocity = {});
    void initialize(std::span<const double> dispersed_fraction, std::span<const Vec3> primary_velocity,
                    std::span<const Vec3> dispersed_velocity);

    template <class F>
    void initialize(F&& fraction_function, Vec3 primary_velocity = {}, Vec3 dispersed_velocity = {}) {
        std::vector<double> fraction(mesh_.cell_count());
        for (std::size_t cell = 0; cell < fraction.size(); ++cell)
            fraction[cell] = static_cast<double>(fraction_function(mesh_.cells()[cell].center));
        std::vector<Vec3> primary(mesh_.cell_count(), primary_velocity);
        std::vector<Vec3> dispersed(mesh_.cell_count(), dispersed_velocity);
        initialize(std::span<const double>(fraction), std::span<const Vec3>(primary), std::span<const Vec3>(dispersed));
    }

    // Advance one explicit conservative finite-volume step. Boundary faces
    // are impermeable in this baseline; internal phase fluxes are upwinded.
    void step();

    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const EulerEulerConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::vector<double>& dispersed_fraction() const noexcept { return dispersed_fraction_; }
    [[nodiscard]] const std::vector<Vec3>& primary_velocity() const noexcept { return primary_velocity_; }
    [[nodiscard]] const std::vector<Vec3>& dispersed_velocity() const noexcept { return dispersed_velocity_; }
    [[nodiscard]] const std::vector<double>& primary_face_flux() const noexcept { return primary_face_flux_; }
    [[nodiscard]] const std::vector<double>& dispersed_face_flux() const noexcept { return dispersed_face_flux_; }
    [[nodiscard]] const std::vector<Vec3>& drag_force() const noexcept { return drag_force_; }
    [[nodiscard]] double volume_fraction_integral() const;
    [[nodiscard]] Vec3 total_momentum() const noexcept;
    [[nodiscard]] Vec3 integrated_drag_force() const noexcept;
    [[nodiscard]] double minimum_fraction() const;
    [[nodiscard]] double maximum_fraction() const;
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }

private:
    PolyMesh mesh_;
    EulerEulerConfig config_;
    std::vector<double> dispersed_fraction_;
    std::vector<Vec3> primary_velocity_, dispersed_velocity_;
    std::vector<double> primary_face_flux_, dispersed_face_flux_;
    std::vector<Vec3> drag_force_;
    double time_{};
    std::size_t steps_{};

    void validate_config() const;
    static void validate_velocity(Vec3 velocity);
};

} // namespace cfd::fvm
