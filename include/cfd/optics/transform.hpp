#pragma once
#include "cfd/solvers/optics/ray.hpp"
#include <array>
namespace cfd::optics {
class RigidTransform {
public:
    static RigidTransform from_euler_xyz(Vec3 translation,double rx,double ry,double rz);
    [[nodiscard]] Vec3 apply_point(Vec3 p) const noexcept;
    [[nodiscard]] Vec3 apply_vector(Vec3 v) const noexcept;
    [[nodiscard]] Vec3 inverse_point(Vec3 p) const noexcept;
    [[nodiscard]] Vec3 inverse_vector(Vec3 v) const noexcept;
    [[nodiscard]] Ray apply(Ray r) const noexcept;
    [[nodiscard]] Ray inverse(Ray r) const noexcept;
private:Vec3 t_{};std::array<double,9> r_{1,0,0,0,1,0,0,0,1};
};
} // namespace cfd::optics
