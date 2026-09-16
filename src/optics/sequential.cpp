#include "cfd/solvers/optics/sequential.hpp"

#include "cfd/core/parallel.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::optics {
namespace {

Vec3 add(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 scale(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
Vec3 sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

// Returns sag and ds/dr; a negative radicand is outside the real conic aperture.
bool sag_and_slope(const SequentialSurface& s,double r,double& sag,double& slope) {
    if(s.type==SurfaceType::plane){sag=0.0;slope=0.0;return true;}
    const double c=1.0/s.radius;
    const double k=s.type==SurfaceType::sphere?0.0:s.conic_constant;
    const double radicand=1.0-(1.0+k)*c*c*r*r;
    if(!(radicand>0.0))return false;
    const double root=std::sqrt(radicand);
    sag=c*r*r/(1.0+root);
    slope=c*r/root;
    if(s.type==SurfaceType::even_asphere){
        for(std::size_t i=0;i<s.even_coefficients.size();++i){
            const int exponent=4+2*static_cast<int>(i);
            sag+=s.even_coefficients[i]*std::pow(r,exponent);
            slope+=static_cast<double>(exponent)*s.even_coefficients[i]*std::pow(r,exponent-1);
        }
    }
    return std::isfinite(sag)&&std::isfinite(slope);
}

bool intersect_surface(const SequentialSurface& surface,
                       const Ray& ray,
                       Vec3& hit,
                       Vec3& normal) {
    constexpr double eps = 1.0e-12;
    double t = 0.0;
    if (surface.type == SurfaceType::plane) {
        if (std::abs(ray.direction.z) < eps) return false;
        t = (surface.vertex_z - ray.origin.z) / ray.direction.z;
        if (!(t > eps)) return false;
        hit = add(ray.origin, scale(ray.direction, t));
        normal = {0.0, 0.0, -1.0};
    } else if(surface.type==SurfaceType::conic||surface.type==SurfaceType::even_asphere) {
        if(std::abs(ray.direction.z)<eps)return false;
        t=(surface.vertex_z-ray.origin.z)/ray.direction.z;
        bool converged=false;
        for(std::size_t iteration=0;iteration<50U;++iteration){
            hit=add(ray.origin,scale(ray.direction,t));
            const double r=std::hypot(hit.x,hit.y);
            double sag=0.0,slope=0.0;
            if(!sag_and_slope(surface,r,sag,slope))return false;
            const double residual=hit.z-surface.vertex_z-sag;
            const double sx=r>0.0?slope*hit.x/r:0.0,sy=r>0.0?slope*hit.y/r:0.0;
            normal=normalized({sx,sy,-1.0});
            if(std::abs(residual)<1.0e-12*(1.0+std::abs(sag))){converged=true;break;}
            const double derivative=ray.direction.z-sx*ray.direction.x-sy*ray.direction.y;
            if(std::abs(derivative)<eps)return false;
            t-=residual/derivative;
            if(!std::isfinite(t))return false;
        }
        if(!converged||!(t>eps))return false;
    } else {
        if (std::abs(surface.radius) < eps) return false;
        const Vec3 center{0.0, 0.0, surface.vertex_z + surface.radius};
        const Vec3 oc = sub(ray.origin, center);
        const double b = dot(oc, ray.direction);
        const double c = dot(oc, oc) - surface.radius * surface.radius;
        const double disc = b * b - c;
        if (disc < 0.0) return false;
        const double root = std::sqrt(disc);
        const double t0 = -b - root;
        const double t1 = -b + root;
        if (t0 > eps) t = t0;
        else if (t1 > eps) t = t1;
        else return false;
        hit = add(ray.origin, scale(ray.direction, t));
        normal = normalized(sub(hit, center));
    }
    return hit.x * hit.x + hit.y * hit.y <= surface.aperture_radius * surface.aperture_radius;
}

} // namespace

double surface_sag(const SequentialSurface& surface,double radial_distance) {
    if(!(radial_distance>=0.0)||!std::isfinite(radial_distance))throw std::invalid_argument("invalid sag radius");
    double sag=0.0,slope=0.0;
    if(!sag_and_slope(surface,radial_distance,sag,slope))throw std::domain_error("surface sag outside real conic");
    return sag;
}

SequentialOpticalSystem::SequentialOpticalSystem(double n) : object_space_index_(n) {
    if (!(n > 0.0)) throw std::invalid_argument("object-space index must be positive");
}

void SequentialOpticalSystem::add_surface(SequentialSurface surface) {
    if (!(surface.aperture_radius > 0.0) || !(surface.refractive_index_after > 0.0)
        || !std::isfinite(surface.vertex_z)
        || !std::isfinite(surface.aperture_radius)||!std::isfinite(surface.refractive_index_after)
        || !std::isfinite(surface.conic_constant)
        || (surface.type != SurfaceType::plane
            && (!std::isfinite(surface.radius) || surface.radius == 0.0))) {
        throw std::invalid_argument("invalid sequential optical surface");
    }
    for(double a:surface.even_coefficients)if(!std::isfinite(a))throw std::invalid_argument("invalid asphere coefficient");
    surfaces_.push_back(surface);
}

RayTraceResult SequentialOpticalSystem::trace(Ray ray) const {
    ray.direction = normalized(ray.direction);
    double n = object_space_index_;
    std::size_t count = 0U;
    for (const auto& surface : surfaces_) {
        Vec3 hit{}, normal{};
        if (!intersect_surface(surface, ray, hit, normal)) return {ray, false, count};
        const auto transmitted = refract(ray.direction, normal, n, surface.refractive_index_after);
        if (!transmitted) return {ray, false, count};
        ray.origin = hit;
        ray.direction = *transmitted;
        n = surface.refractive_index_after;
        ++count;
    }
    return {ray, true, count};
}

std::vector<RayTraceResult> SequentialOpticalSystem::trace_many(const std::vector<Ray>& rays) const {
    std::vector<RayTraceResult> out(rays.size());
    cfd::core::parallel_for(rays.size(), [&](std::size_t i) { out[i] = trace(rays[i]); });
    return out;
}

double rms_spot_radius_at_plane(const std::vector<RayTraceResult>& rays, double z) {
    double sx = 0.0;
    double sy = 0.0;
    std::size_t n = 0U;
    for (const auto& result : rays) {
        if (!result.valid || std::abs(result.ray.direction.z) <= 1.0e-14) continue;
        const double t = (z - result.ray.origin.z) / result.ray.direction.z;
        if (t < 0.0) continue;
        sx += result.ray.origin.x + t * result.ray.direction.x;
        sy += result.ray.origin.y + t * result.ray.direction.y;
        ++n;
    }
    if (n == 0U) throw std::runtime_error("no valid rays reach image plane");
    const double cx = sx / static_cast<double>(n);
    const double cy = sy / static_cast<double>(n);
    double sum = 0.0;
    for (const auto& result : rays) {
        if (!result.valid || std::abs(result.ray.direction.z) <= 1.0e-14) continue;
        const double t = (z - result.ray.origin.z) / result.ray.direction.z;
        if (t < 0.0) continue;
        const double x = result.ray.origin.x + t * result.ray.direction.x;
        const double y = result.ray.origin.y + t * result.ray.direction.y;
        sum += (x - cx) * (x - cx) + (y - cy) * (y - cy);
    }
    return std::sqrt(sum / static_cast<double>(n));
}

} // namespace cfd::optics
