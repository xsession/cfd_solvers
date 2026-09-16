#pragma once

#include <optional>

namespace cfd::optics {

struct Vec3 {
    double x{};
    double y{};
    double z{};
};

struct Ray {
    Vec3 origin{};
    Vec3 direction{0.0, 0.0, 1.0};
    double wavelength_nm{550.0};
};

[[nodiscard]] double dot(Vec3 a, Vec3 b) noexcept;
[[nodiscard]] Vec3 normalized(Vec3 v);
[[nodiscard]] std::optional<Vec3> refract(Vec3 incident,
                                          Vec3 normal,
                                          double n_from,
                                          double n_to);
[[nodiscard]] Vec3 reflect(Vec3 incident, Vec3 normal);

} // namespace cfd::optics
