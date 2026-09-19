#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <span>
#include <algorithm>
#include <vector>
namespace cfd::fvm {
struct VofConfig {
    double dt{1e-3};
    // Blend from bounded-linear transport (0) toward the bounded upwind
    // interface flux (1). Keeping the blend bounded is the inexpensive,
    // conservative compression baseline used by VofTransport.
    double compression{1.0};
};
class VofTransport {
public:
    VofTransport(PolyMesh mesh, VofConfig config = {});
    void initialize(double alpha);
    template <class F> void initialize(F f) {
        for (std::size_t c = 0; c < alpha_.size(); ++c)
            alpha_[c] = std::clamp(static_cast<double>(f(mesh_.cells()[c].center)), 0.0, 1.0);
    }
    void set_face_flux(std::vector<double> flux);
    void step();
    [[nodiscard]] const std::vector<double>& alpha() const noexcept { return alpha_; }
    [[nodiscard]] double volume() const;
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }

private:
    PolyMesh mesh_;
    VofConfig config_;
    std::vector<double> alpha_, flux_;
};
[[nodiscard]] std::vector<Vec3> csf_surface_tension(const PolyMesh& mesh, std::span<const double> alpha,
                                                    double surface_tension);
[[nodiscard]] Vec3 contact_angle_normal(Vec3 interface_normal, Vec3 wall_normal, double angle_radians);
} // namespace cfd::fvm
