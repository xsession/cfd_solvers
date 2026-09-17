#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::multibody {

struct Vec3 {
    double x{};
    double y{};
    double z{};
};

[[nodiscard]] constexpr Vec3 operator+(Vec3 a, Vec3 b) noexcept { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
[[nodiscard]] constexpr Vec3 operator-(Vec3 a, Vec3 b) noexcept { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
[[nodiscard]] constexpr Vec3 operator-(Vec3 a) noexcept { return {-a.x,-a.y,-a.z}; }
[[nodiscard]] constexpr Vec3 operator*(Vec3 a,double s) noexcept { return {a.x*s,a.y*s,a.z*s}; }
[[nodiscard]] constexpr Vec3 operator*(double s,Vec3 a) noexcept { return a*s; }
[[nodiscard]] constexpr Vec3 operator/(Vec3 a,double s) noexcept { return {a.x/s,a.y/s,a.z/s}; }
constexpr Vec3& operator+=(Vec3& a,Vec3 b) noexcept { a=a+b; return a; }
constexpr Vec3& operator-=(Vec3& a,Vec3 b) noexcept { a=a-b; return a; }
constexpr Vec3& operator*=(Vec3& a,double s) noexcept { a=a*s; return a; }

[[nodiscard]] constexpr double dot(Vec3 a,Vec3 b) noexcept { return a.x*b.x+a.y*b.y+a.z*b.z; }
[[nodiscard]] constexpr Vec3 cross(Vec3 a,Vec3 b) noexcept {
    return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}
[[nodiscard]] inline double norm_squared(Vec3 a) noexcept { return dot(a,a); }
[[nodiscard]] inline double norm(Vec3 a) noexcept { return std::sqrt(norm_squared(a)); }
[[nodiscard]] inline Vec3 normalized(Vec3 a,double eps=1e-15) noexcept {
    const double n=norm(a);
    return n>eps ? a/n : Vec3{};
}
[[nodiscard]] inline Vec3 component_mul(Vec3 a,Vec3 b) noexcept { return {a.x*b.x,a.y*b.y,a.z*b.z}; }

struct Quaternion {
    double w{1.0};
    double x{};
    double y{};
    double z{};
};

[[nodiscard]] constexpr Quaternion conjugate(Quaternion q) noexcept { return {q.w,-q.x,-q.y,-q.z}; }
[[nodiscard]] constexpr Quaternion operator*(Quaternion a,Quaternion b) noexcept {
    return {
        a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z,
        a.w*b.x+a.x*b.w+a.y*b.z-a.z*b.y,
        a.w*b.y-a.x*b.z+a.y*b.w+a.z*b.x,
        a.w*b.z+a.x*b.y-a.y*b.x+a.z*b.w
    };
}
[[nodiscard]] inline double norm_squared(Quaternion q) noexcept { return q.w*q.w+q.x*q.x+q.y*q.y+q.z*q.z; }
[[nodiscard]] inline Quaternion normalized(Quaternion q) {
    const double n=std::sqrt(norm_squared(q));
    if(!(n>0.0) || !std::isfinite(n)) throw std::invalid_argument("invalid zero/non-finite quaternion");
    return {q.w/n,q.x/n,q.y/n,q.z/n};
}
[[nodiscard]] inline Vec3 rotate(Quaternion q,Vec3 v) {
    q=normalized(q);
    const Quaternion p{0.0,v.x,v.y,v.z};
    const Quaternion r=q*p*conjugate(q);
    return {r.x,r.y,r.z};
}
[[nodiscard]] inline Vec3 inverse_rotate(Quaternion q,Vec3 v) { return rotate(conjugate(normalized(q)),v); }
[[nodiscard]] inline Quaternion quaternion_from_axis_angle(Vec3 axis,double angle) {
    const Vec3 n=normalized(axis);
    const double h=0.5*angle;
    const double s=std::sin(h);
    return normalized({std::cos(h),n.x*s,n.y*s,n.z*s});
}
inline void integrate_orientation(Quaternion& q,Vec3 angular_velocity_world,double dt) {
    const Quaternion omega{0.0,angular_velocity_world.x,angular_velocity_world.y,angular_velocity_world.z};
    const Quaternion dq=omega*q;
    q=normalized({q.w+0.5*dt*dq.w,q.x+0.5*dt*dq.x,q.y+0.5*dt*dq.y,q.z+0.5*dt*dq.z});
}

[[nodiscard]] inline Vec3 any_perpendicular(Vec3 axis) noexcept {
    const Vec3 n=normalized(axis);
    const Vec3 seed=(std::abs(n.x)<0.7)?Vec3{1.0,0.0,0.0}:Vec3{0.0,1.0,0.0};
    return normalized(cross(n,seed));
}

} // namespace cfd::multibody
