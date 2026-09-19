#pragma once

#include "cfd/solvers/fvm/particles.hpp"

#include <span>
#include <vector>

namespace cfd::fvm {

struct DenseRheologyConfig {
    double fluid_dynamic_viscosity{1.8e-5};
    double particle_density{1000.0};
    double grain_diameter{1.0e-4};
    double maximum_solid_fraction{0.64};
    double static_friction{0.38};
    double dynamic_friction{0.64};
    double inertial_scale{0.30};
    double regularization{1.0e-9};
};

struct DenseRheologyResult {
    std::vector<double> solid_fraction;
    std::vector<double> granular_pressure;
    std::vector<double> effective_viscosity;
    std::vector<Vec3> carrier_momentum_source;
    std::vector<Vec3> particle_momentum_source;
};

// Compact Eulerian-Lagrangian dense-phase closure. It combines a
// Richardson-Zaki/Ergun-like hindered drag with a regularized mu(I) granular
// stress estimate and returns equal-and-opposite carrier/particle sources.
class DenseParticleRheology {
public:
    explicit DenseParticleRheology(PolyMesh mesh, DenseRheologyConfig config = {});

    [[nodiscard]] DenseRheologyResult evaluate(std::span<const Particle> particles,
                                               std::span<const Vec3> carrier_velocity) const;
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const DenseRheologyConfig& config() const noexcept { return config_; }

private:
    PolyMesh mesh_;
    DenseRheologyConfig config_;

    void validate_config() const;
};

} // namespace cfd::fvm
