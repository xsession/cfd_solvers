#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <algorithm>
#include <span>
#include <vector>

namespace cfd::fvm {

// Explicit lubrication-film controls. The thickness field is coupled to a
// capillary pressure reconstructed from its interface curvature and to a
// body-force potential. The mesh is treated as a unit-depth control volume,
// which makes the baseline useful for 2-D film and coating cases.
struct ThinFilmConfig {
    double dt{1.0e-4};
    double density{1000.0};
    double viscosity{1.0e-3};
    double surface_tension{0.0};
    Vec3 body_acceleration{};
    double minimum_thickness{1.0e-9};
};

class ThinFilmTransport {
public:
    explicit ThinFilmTransport(PolyMesh mesh, ThinFilmConfig config = {});

    void initialize(double thickness);

    template <class F> void initialize(F&& function) {
        for (std::size_t cell = 0; cell < thickness_.size(); ++cell) {
            thickness_[cell] =
                std::max(config_.minimum_thickness, static_cast<double>(function(mesh_.cells()[cell].center)));
        }
        std::fill(pressure_.begin(), pressure_.end(), 0.0);
        std::fill(face_flux_.begin(), face_flux_.end(), 0.0);
    }

    // Optional fixed pressure on boundary faces. With no values supplied the
    // solver uses impermeable boundaries, making the film inventory closed.
    void set_boundary_pressure(std::vector<double> pressure);
    void clear_boundary_pressure() noexcept { boundary_pressure_.clear(); }

    // Rebuild capillary pressure, evaluate the lubrication flux, and advance
    // the thickness by one explicit conservative finite-volume step.
    void step();

    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& thickness() const noexcept { return thickness_; }
    [[nodiscard]] const std::vector<double>& pressure() const noexcept { return pressure_; }
    [[nodiscard]] const std::vector<double>& face_flux() const noexcept { return face_flux_; }
    [[nodiscard]] double inventory() const;

private:
    PolyMesh mesh_;
    ThinFilmConfig config_;
    std::vector<double> thickness_;
    std::vector<double> pressure_;
    std::vector<double> face_flux_;
    std::vector<double> boundary_pressure_;

    void update_pressure();
    void update_face_flux();
};

} // namespace cfd::fvm
