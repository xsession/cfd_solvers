#include "cfd/solvers/lbm/particles.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace cfd::lbm {
namespace {
std::size_t wrap(long long i, std::size_t n) {
    const auto m = static_cast<long long>(n);
    i %= m;
    if (i < 0) i += m;
    return static_cast<std::size_t>(i);
}

std::size_t index(std::size_t x, std::size_t y, std::size_t z,
                  std::size_t nx, std::size_t ny) {
    return (z * ny + y) * nx + x;
}

struct KernelPoint {
    std::size_t index{};
    double weight{};
};

std::array<KernelPoint, 8> trilinear_kernel(std::size_t nx, std::size_t ny, std::size_t nz,
                                             ParticleVector3 p) {
    const double fx0 = std::floor(p.x);
    const double fy0 = std::floor(p.y);
    const double fz0 = std::floor(p.z);
    const auto x0 = static_cast<long long>(fx0);
    const auto y0 = static_cast<long long>(fy0);
    const auto z0 = static_cast<long long>(fz0);
    const double tx = p.x - fx0;
    const double ty = p.y - fy0;
    const double tz = p.z - fz0;
    std::array<KernelPoint, 8> points{};
    std::size_t n = 0;
    for (int dz = 0; dz <= 1; ++dz) for (int dy = 0; dy <= 1; ++dy) for (int dx = 0; dx <= 1; ++dx) {
        const double wx = dx ? tx : 1.0 - tx;
        const double wy = dy ? ty : 1.0 - ty;
        const double wz = dz ? tz : 1.0 - tz;
        const auto x = wrap(x0 + dx, nx);
        const auto y = wrap(y0 + dy, ny);
        const auto z = wrap(z0 + dz, nz);
        points[n++] = {index(x, y, z, nx, ny), wx * wy * wz};
    }
    return points;
}

void validate_fields(std::size_t nx, std::size_t ny, std::size_t nz,
                     std::span<const float> ux, std::span<const float> uy, std::span<const float> uz) {
    if (nx < 2 || ny < 2 || nz < 2) throw std::invalid_argument("particle coupling grid must be at least 2x2x2");
    const std::size_t cells = nx * ny * nz;
    if (ux.size() != cells || uy.size() != cells || uz.size() != cells) {
        throw std::invalid_argument("particle coupling velocity field size mismatch");
    }
}

void periodic_wrap(ParticleVector3& p, std::size_t nx, std::size_t ny, std::size_t nz) {
    const auto fold = [](double value, double extent) {
        value = std::fmod(value, extent);
        return value < 0.0 ? value + extent : value;
    };
    p.x = fold(p.x, static_cast<double>(nx));
    p.y = fold(p.y, static_cast<double>(ny));
    p.z = fold(p.z, static_cast<double>(nz));
}
} // namespace

ParticleVector3 interpolate_particle_velocity(
    std::size_t nx, std::size_t ny, std::size_t nz,
    std::span<const float> ux, std::span<const float> uy, std::span<const float> uz,
    ParticleVector3 position) {
    validate_fields(nx, ny, nz, ux, uy, uz);
    periodic_wrap(position, nx, ny, nz);
    ParticleVector3 result{};
    for (const auto point : trilinear_kernel(nx, ny, nz, position)) {
        result.x += point.weight * static_cast<double>(ux[point.index]);
        result.y += point.weight * static_cast<double>(uy[point.index]);
        result.z += point.weight * static_cast<double>(uz[point.index]);
    }
    return result;
}

ParticleCouplingResult advance_immersed_boundary_particles(
    std::size_t nx, std::size_t ny, std::size_t nz,
    std::span<const float> ux, std::span<const float> uy, std::span<const float> uz,
    std::span<ImmersedBoundaryParticle> particles,
    double dt,
    bool two_way_coupling,
    ParticleVector3 gravity) {
    validate_fields(nx, ny, nz, ux, uy, uz);
    if (!(dt > 0.0)) throw std::invalid_argument("particle coupling dt must be positive");
    const std::size_t cells = nx * ny * nz;
    ParticleCouplingResult result;
    if (two_way_coupling) {
        result.acceleration_x.assign(cells, 0.0F);
        result.acceleration_y.assign(cells, 0.0F);
        result.acceleration_z.assign(cells, 0.0F);
    }

    for (auto& particle : particles) {
        if (!(particle.mass > 0.0) || !(particle.response_time > 0.0)) {
            throw std::invalid_argument("particle mass and response time must be positive");
        }
        periodic_wrap(particle.position, nx, ny, nz);
        const auto kernel = trilinear_kernel(nx, ny, nz, particle.position);
        ParticleVector3 fluid{};
        for (const auto point : kernel) {
            fluid.x += point.weight * static_cast<double>(ux[point.index]);
            fluid.y += point.weight * static_cast<double>(uy[point.index]);
            fluid.z += point.weight * static_cast<double>(uz[point.index]);
        }
        const ParticleVector3 drag{
            particle.mass * (fluid.x - particle.velocity.x) / particle.response_time,
            particle.mass * (fluid.y - particle.velocity.y) / particle.response_time,
            particle.mass * (fluid.z - particle.velocity.z) / particle.response_time};

        particle.velocity.x += dt * (drag.x / particle.mass + gravity.x);
        particle.velocity.y += dt * (drag.y / particle.mass + gravity.y);
        particle.velocity.z += dt * (drag.z / particle.mass + gravity.z);
        particle.position.x += dt * particle.velocity.x;
        particle.position.y += dt * particle.velocity.y;
        particle.position.z += dt * particle.velocity.z;
        periodic_wrap(particle.position, nx, ny, nz);

        if (two_way_coupling) {
            result.total_reaction_force.x -= drag.x;
            result.total_reaction_force.y -= drag.y;
            result.total_reaction_force.z -= drag.z;
            for (const auto point : kernel) {
                result.acceleration_x[point.index] += static_cast<float>(-drag.x * point.weight);
                result.acceleration_y[point.index] += static_cast<float>(-drag.y * point.weight);
                result.acceleration_z[point.index] += static_cast<float>(-drag.z * point.weight);
            }
        }
    }
    return result;
}

} // namespace cfd::lbm
