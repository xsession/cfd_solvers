#include "cfd/solvers/optics/ray.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::optics {

double dot(Vec3 a, Vec3 b) noexcept { return a.x*b.x + a.y*b.y + a.z*b.z; }

Vec3 normalized(Vec3 v) {
    const double n = std::sqrt(dot(v, v));
    if (n == 0.0) throw std::invalid_argument("cannot normalize zero vector");
    return {v.x/n, v.y/n, v.z/n};
}

Vec3 reflect(Vec3 incident, Vec3 normal) {
    const Vec3 i = normalized(incident);
    Vec3 n = normalized(normal);
    if (dot(i, n) > 0.0) n = {-n.x, -n.y, -n.z};
    const double d = dot(i, n);
    return normalized({i.x - 2.0*d*n.x, i.y - 2.0*d*n.y, i.z - 2.0*d*n.z});
}

std::optional<Vec3> refract(Vec3 incident, Vec3 normal, double n_from, double n_to) {
    if (n_from <= 0.0 || n_to <= 0.0) throw std::invalid_argument("refractive indices must be positive");
    const Vec3 i = normalized(incident);
    Vec3 n = normalized(normal);
    double cos_i = -dot(i, n);
    if (cos_i < 0.0) {
        n = {-n.x, -n.y, -n.z};
        cos_i = -dot(i, n);
    }
    cos_i = std::clamp(cos_i, 0.0, 1.0);
    const double eta = n_from / n_to;
    const double k = 1.0 - eta*eta*(1.0 - cos_i*cos_i);
    if (k < 0.0) return std::nullopt;
    const double a = eta*cos_i - std::sqrt(k);
    return normalized({eta*i.x + a*n.x, eta*i.y + a*n.y, eta*i.z + a*n.z});
}

} // namespace cfd::optics
