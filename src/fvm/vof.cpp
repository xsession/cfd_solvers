#include "cfd/solvers/fvm/vof.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/core/parallel.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
VofTransport::VofTransport(PolyMesh m, VofConfig c)
    : mesh_(std::move(m)), config_(c), alpha_(mesh_.cell_count()), flux_(mesh_.face_count()) {
    if (!(c.dt > 0) || c.compression < 0)
        throw std::invalid_argument("invalid VOF controls");
}
void VofTransport::initialize(double a) {
    if (a < 0 || a > 1)
        throw std::invalid_argument("invalid volume fraction");
    std::fill(alpha_.begin(), alpha_.end(), a);
}
void VofTransport::set_face_flux(std::vector<double> f) {
    if (f.size() != mesh_.face_count())
        throw std::invalid_argument("VOF flux size");
    flux_ = std::move(f);
}
double VofTransport::volume() const {
    return cfd::core::parallel_sum(alpha_.size(), [&](std::size_t c) { return alpha_[c] * mesh_.cells()[c].volume; });
}
void VofTransport::step() {
    auto face = interpolate_scalar_to_faces(mesh_, alpha_, FaceInterpolationScheme::bounded_linear, flux_);
    if (config_.compression > 0.0) {
        const auto upwind = interpolate_scalar_to_faces(mesh_, alpha_, FaceInterpolationScheme::upwind, flux_);
        const double blend = std::clamp(config_.compression, 0.0, 1.0);
        for (std::size_t fi = 0; fi < face.size(); ++fi)
            face[fi] = std::clamp((1.0 - blend) * face[fi] + blend * upwind[fi], 0.0, 1.0);
    }
    std::vector<double> next = alpha_;
    for (std::size_t c = 0; c < alpha_.size(); ++c) {
        double sum = 0;
        for (auto fi : mesh_.cell_faces()[c]) {
            const auto& f = mesh_.faces()[fi];
            double oriented = f.owner == c ? flux_[fi] : -flux_[fi];
            sum += oriented * face[fi];
        }
        next[c] = std::clamp(alpha_[c] - config_.dt * sum / mesh_.cells()[c].volume, 0.0, 1.0);
    }
    alpha_.swap(next);
}
std::vector<Vec3> csf_surface_tension(const PolyMesh& m, std::span<const double> a, double sigma) {
    if (a.size() != m.cell_count() || sigma < 0)
        throw std::invalid_argument("invalid CSF inputs");
    auto grad = least_squares_gradient_scalar(m, a);
    std::vector<Vec3> n(a.size());
    for (std::size_t c = 0; c < a.size(); ++c) {
        double mag = magnitude(grad[c]);
        if (mag > 1e-14)
            n[c] = grad[c] / mag;
    }
    auto div = gauss_divergence_vector(m, n);
    std::vector<Vec3> f(a.size());
    for (std::size_t c = 0; c < a.size(); ++c)
        f[c] = grad[c] * (-sigma * div[c]);
    return f;
}
Vec3 contact_angle_normal(Vec3 ni, Vec3 nw, double angle) {
    double wn = magnitude(nw), in = magnitude(ni);
    if (!(wn > 0) || !(in > 0) || angle < 0 || angle > 3.14159265358979323846)
        throw std::invalid_argument("invalid contact angle");
    nw = nw / wn;
    ni = ni / in;
    Vec3 tangent = ni - nw * dot(ni, nw);
    double tm = magnitude(tangent);
    if (tm < 1e-14) {
        tangent = std::abs(nw.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        tangent = tangent - nw * dot(tangent, nw);
        tm = magnitude(tangent);
    }
    tangent = tangent / tm;
    return nw * std::cos(angle) + tangent * std::sin(angle);
}
} // namespace cfd::fvm
