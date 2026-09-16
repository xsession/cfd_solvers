#include "cfd/solvers/fvm/collocated_incompressible.hpp"

#include "cfd/core/parallel.hpp"
#include "cfd/fvm/operators.hpp"

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>

namespace cfd::fvm {
namespace {

[[nodiscard]] double interpolate_scalar(const PolyMesh& mesh,
                                        const Face& face,
                                        std::span<const double> values) {
    if (face.boundary()) return values[face.owner];
    const double dof = magnitude(face.center - mesh.cells()[face.owner].center);
    const double dnf = magnitude(mesh.cells()[face.neighbour].center - face.center);
    const double d = dof + dnf;
    if (!(d > 0.0)) throw std::runtime_error("degenerate face interpolation distance");
    return (dnf * values[face.owner] + dof * values[face.neighbour]) / d;
}

[[nodiscard]] Vec3 interpolate_vector(const PolyMesh& mesh,
                                      const Face& face,
                                      std::span<const Vec3> values) {
    if (face.boundary()) return values[face.owner];
    const double dof = magnitude(face.center - mesh.cells()[face.owner].center);
    const double dnf = magnitude(mesh.cells()[face.neighbour].center - face.center);
    const double d = dof + dnf;
    if (!(d > 0.0)) throw std::runtime_error("degenerate face interpolation distance");
    return (values[face.owner] * dnf + values[face.neighbour] * dof) / d;
}

[[nodiscard]] double outward_flux(const Face& face, std::size_t cell, double phi) noexcept {
    return face.owner == cell ? phi : -phi;
}

[[nodiscard]] Vec3 boundary_velocity_for_momentum(const Face& face,
                                                   Vec3 owner,
                                                   const VelocityBoundaryCondition& bc) {
    if (bc.type == VelocityBoundaryType::fixedValue) return bc.value;
    if (bc.type == VelocityBoundaryType::zeroGradient) return owner;
    const double area_mag = magnitude(face.area);
    const Vec3 n = face.area / area_mag;
    const Vec3 relative = owner - bc.value;
    return bc.value + relative - n * dot(relative, n);
}

} // namespace

CollocatedIncompressible::CollocatedIncompressible(PolyMesh mesh,
                                                   CollocatedIncompressibleConfig config)
    : mesh_(std::move(mesh)),
      config_(config),
      velocity_boundary_(mesh_.patches().size()),
      pressure_boundary_(mesh_.patches().size()),
      velocity_(mesh_.cell_count()),
      old_velocity_(mesh_.cell_count()),
      h_by_a_(mesh_.cell_count()),
      pressure_mobility_(mesh_.cell_count(), config.dt),
      pressure_(mesh_.cell_count(), 0.0),
      face_flux_(mesh_.face_count(), 0.0),
      pressure_rhs_(mesh_.cell_count(), 0.0),
      pressure_face_coefficient_(mesh_.face_count(), 0.0) {
    if (!(config_.density > 0.0) || !(config_.kinematic_viscosity >= 0.0) || !(config_.dt > 0.0) ||
        config_.momentum_sweeps == 0U || config_.pressure_iterations == 0U ||
        !(config_.pressure_tolerance > 0.0) || config_.pressure_correctors == 0U ||
        config_.outer_correctors == 0U || !(config_.velocity_relaxation > 0.0 && config_.velocity_relaxation <= 1.0) ||
        !(config_.pressure_relaxation > 0.0 && config_.pressure_relaxation <= 1.0)) {
        throw std::invalid_argument("invalid collocated incompressible solver controls");
    }
    rebuild_flux_from_velocity();
}

void CollocatedIncompressible::set_velocity_boundary(std::string_view patch,
                                                      VelocityBoundaryType type,
                                                      Vec3 value) {
    velocity_boundary_[mesh_.patch_index(patch)] = {type, value};
    rebuild_flux_from_velocity();
}

void CollocatedIncompressible::set_pressure_boundary(std::string_view patch,
                                                      PressureBoundaryType type,
                                                      double value) {
    pressure_boundary_[mesh_.patch_index(patch)] = {type, value};
}

void CollocatedIncompressible::initialize_uniform(Vec3 velocity, double pressure) {
    std::fill(velocity_.begin(), velocity_.end(), velocity);
    std::fill(old_velocity_.begin(), old_velocity_.end(), velocity);
    std::fill(h_by_a_.begin(), h_by_a_.end(), velocity);
    std::fill(pressure_.begin(), pressure_.end(), pressure);
    time_ = 0.0;
    steps_ = 0U;
    rebuild_flux_from_velocity();
}

void CollocatedIncompressible::initialize_fields(const std::function<Vec3(Vec3)>& velocity,
                                                  const std::function<double(Vec3)>& pressure) {
    for (std::size_t c = 0; c < mesh_.cell_count(); ++c) {
        velocity_[c] = velocity(mesh_.cells()[c].center);
        pressure_[c] = pressure(mesh_.cells()[c].center);
    }
    old_velocity_ = velocity_;
    h_by_a_ = velocity_;
    time_ = 0.0;
    steps_ = 0U;
    rebuild_flux_from_velocity();
}

bool CollocatedIncompressible::has_fixed_pressure_boundary() const noexcept {
    for (const auto& bc : pressure_boundary_) {
        if (bc.type == PressureBoundaryType::fixedValue) return true;
    }
    return false;
}

void CollocatedIncompressible::rebuild_flux_from_velocity() {
    const auto boundary = boundary_velocity_values(mesh_, velocity_, velocity_boundary_);
    face_flux_ = face_flux_from_velocity(mesh_, velocity_, boundary);
}

void CollocatedIncompressible::momentum_predictor(std::span<const Vec3> time_source) {
    if (time_source.size() != mesh_.cell_count()) throw std::invalid_argument("time source size mismatch");
    const auto p_boundary = boundary_pressure_values(mesh_, pressure_, pressure_boundary_);
    const auto grad_p = gauss_gradient_scalar(mesh_, pressure_, p_boundary);
    std::vector<Vec3> working = velocity_;
    std::vector<Vec3> next(mesh_.cell_count());
    std::vector<double> diagonal(mesh_.cell_count(), 0.0);
    std::vector<Vec3> h(mesh_.cell_count());

    for (std::size_t sweep = 0; sweep < config_.momentum_sweeps; ++sweep) {
        cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
            const double volume = mesh_.cells()[cell].volume;
            double ap = volume / config_.dt;
            Vec3 rhs = time_source[cell] * (volume / config_.dt);

            for (const std::size_t fi : mesh_.cell_faces()[cell]) {
                const auto& face = mesh_.faces()[fi];
                const double phi_out = outward_flux(face, cell, face_flux_[fi]);
                const auto geometry = decompose_face_area(mesh_, fi);
                const double diffusion = config_.kinematic_viscosity * geometry.orthogonal_metric;

                if (!face.boundary()) {
                    const std::size_t other = face.owner == cell ? face.neighbour : face.owner;
                    const double entering = config_.include_convection ? std::max(-phi_out, 0.0) : 0.0;
                    const double leaving = config_.include_convection ? std::max(phi_out, 0.0) : 0.0;
                    ap += diffusion + leaving;
                    rhs += working[other] * (diffusion + entering);
                } else {
                    const auto& bc = velocity_boundary_[face.patch];
                    if (bc.type == VelocityBoundaryType::fixedValue) {
                        const Vec3 ub = boundary_velocity_for_momentum(face, working[cell], bc);
                        const double entering = config_.include_convection ? std::max(-phi_out, 0.0) : 0.0;
                        const double leaving = config_.include_convection ? std::max(phi_out, 0.0) : 0.0;
                        ap += diffusion + leaving;
                        rhs += ub * (diffusion + entering);
                    } else if (config_.include_convection) {
                        // zero-gradient/slip boundaries use the owner value for tangential convection.
                        ap += phi_out;
                    }
                }
            }

            if (!(ap > 0.0) || !std::isfinite(ap)) throw std::runtime_error("non-positive momentum diagonal");
            diagonal[cell] = ap;
            h[cell] = rhs;
            const Vec3 hba = rhs / ap;
            const double mobility = mesh_.cells()[cell].volume / ap;
            const Vec3 candidate = hba - grad_p[cell] * mobility;
            next[cell] = candidate;
        });
        working.swap(next);
    }

    cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
        h_by_a_[cell] = h[cell] / diagonal[cell];
        pressure_mobility_[cell] = mesh_.cells()[cell].volume / diagonal[cell];
        velocity_[cell] = h_by_a_[cell] - grad_p[cell] * pressure_mobility_[cell];
    });
}

std::vector<double> CollocatedIncompressible::predicted_flux() const {
    const auto boundary = boundary_velocity_values(mesh_, h_by_a_, velocity_boundary_);
    return face_flux_from_velocity(mesh_, h_by_a_, boundary);
}

void CollocatedIncompressible::update_pressure_coefficients() {
    cfd::core::parallel_for(mesh_.face_count(), [&](std::size_t fi) {
        const auto& face = mesh_.faces()[fi];
        if (face.boundary() && pressure_boundary_[face.patch].type == PressureBoundaryType::zeroGradient) {
            pressure_face_coefficient_[fi] = 0.0;
            return;
        }
        const double mobility = interpolate_scalar(mesh_, face, pressure_mobility_);
        pressure_face_coefficient_[fi] = mobility * decompose_face_area(mesh_, fi).orthogonal_metric;
    });
}

std::vector<double> CollocatedIncompressible::nonorthogonal_pressure_flux(std::span<const double> p) const {
    const auto p_boundary = boundary_pressure_values(mesh_, p, pressure_boundary_);
    const auto grad = gauss_gradient_scalar(mesh_, p, p_boundary);
    std::vector<double> correction(mesh_.face_count(), 0.0);
    cfd::core::parallel_for(mesh_.face_count(), [&](std::size_t fi) {
        const auto& face = mesh_.faces()[fi];
        if (face.boundary() && pressure_boundary_[face.patch].type == PressureBoundaryType::zeroGradient) return;
        const auto geometry = decompose_face_area(mesh_, fi);
        const double mobility = interpolate_scalar(mesh_, face, pressure_mobility_);
        const Vec3 gf = face.boundary() ? grad[face.owner] : interpolate_vector(mesh_, face, grad);
        correction[fi] = mobility * dot(gf, geometry.nonorthogonal_area);
    });
    return correction;
}

void CollocatedIncompressible::apply_pressure_operator(std::span<const double> x,
                                                        std::span<double> out) const {
    cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
        double sum = 0.0;
        for (const std::size_t fi : mesh_.cell_faces()[cell]) {
            const auto& face = mesh_.faces()[fi];
            const double coeff = pressure_face_coefficient_[fi];
            if (!(coeff > 0.0)) continue;
            if (face.boundary()) {
                sum += coeff * x[cell];
            } else {
                const std::size_t other = face.owner == cell ? face.neighbour : face.owner;
                sum += coeff * (x[cell] - x[other]);
            }
        }
        out[cell] = sum;
    });
}

cfd::core::ConjugateGradientResult CollocatedIncompressible::correct_pressure(double pressure_relaxation) {
    update_pressure_coefficients();
    const auto phi_h = predicted_flux();
    std::vector<double> p_work = pressure_;
    const bool fixed_pressure = has_fixed_pressure_boundary();
    cfd::core::ConjugateGradientResult result{};

    const std::size_t correction_loops = std::max<std::size_t>(1U, config_.nonorthogonal_correctors + 1U);
    for (std::size_t nonorth = 0; nonorth < correction_loops; ++nonorth) {
        const auto correction = nonorthogonal_pressure_flux(p_work);
        const auto p_boundary = boundary_pressure_values(mesh_, p_work, pressure_boundary_);

        cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
            double net_phi = 0.0;
            double net_correction = 0.0;
            double boundary_source = 0.0;
            for (const std::size_t fi : mesh_.cell_faces()[cell]) {
                const auto& face = mesh_.faces()[fi];
                net_phi += outward_flux(face, cell, phi_h[fi]);
                net_correction += outward_flux(face, cell, correction[fi]);
                if (face.boundary() && pressure_boundary_[face.patch].type == PressureBoundaryType::fixedValue) {
                    boundary_source += pressure_face_coefficient_[fi] * p_boundary[fi];
                }
            }
            pressure_rhs_[cell] = -net_phi + net_correction + boundary_source;
        });

        if (!fixed_pressure) {
            double mean = cfd::core::parallel_sum(mesh_.cell_count(), [&](std::size_t cell) {
                return pressure_rhs_[cell];
            });
            mean /= static_cast<double>(mesh_.cell_count());
            cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
                pressure_rhs_[cell] -= mean;
            });
            double pmean = cfd::core::parallel_sum(mesh_.cell_count(), [&](std::size_t cell) {
                return p_work[cell];
            });
            pmean /= static_cast<double>(mesh_.cell_count());
            cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) { p_work[cell] -= pmean; });
        }

        result = cfd::core::conjugate_gradient(
            std::span<const double>(pressure_rhs_.data(), pressure_rhs_.size()),
            std::span<double>(p_work.data(), p_work.size()),
            [&](std::span<const double> x, std::span<double> out) { apply_pressure_operator(x, out); },
            pressure_workspace_, config_.pressure_iterations, config_.pressure_tolerance);
        if (!result.converged) break;
    }

    const std::vector<double> old_pressure = pressure_;
    cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
        pressure_[cell] = old_pressure[cell] + pressure_relaxation * (p_work[cell] - old_pressure[cell]);
    });
    if (!fixed_pressure) {
        double mean = cfd::core::parallel_sum(mesh_.cell_count(), [&](std::size_t cell) { return pressure_[cell]; });
        mean /= static_cast<double>(mesh_.cell_count());
        cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) { pressure_[cell] -= mean; });
    }

    const auto p_boundary = boundary_pressure_values(mesh_, pressure_, pressure_boundary_);
    const auto grad_p = gauss_gradient_scalar(mesh_, pressure_, p_boundary);
    cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
        const Vec3 corrected = h_by_a_[cell] - grad_p[cell] * pressure_mobility_[cell];
        velocity_[cell] = velocity_[cell] + (corrected - velocity_[cell]) *
            (pressure_relaxation < 1.0 ? config_.velocity_relaxation : 1.0);
    });
    face_flux_ = rhie_chow_face_flux(mesh_, h_by_a_, pressure_mobility_, pressure_,
                                     velocity_boundary_, pressure_boundary_, true);
    pressure_result_ = result;
    return result;
}

double CollocatedIncompressible::velocity_rms_change(std::span<const Vec3> before) const {
    const double sum = cfd::core::parallel_sum(mesh_.cell_count(), [&](std::size_t cell) {
        const Vec3 d = velocity_[cell] - before[cell];
        return dot(d, d);
    });
    return std::sqrt(sum / static_cast<double>(mesh_.cell_count()));
}

CollocatedIterationInfo CollocatedIncompressible::coupled_sequence(std::span<const Vec3> time_source,
                                                                   std::size_t pressure_correctors,
                                                                   double pressure_relaxation) {
    const std::vector<Vec3> before = velocity_;
    momentum_predictor(time_source);
    cfd::core::ConjugateGradientResult result{};
    for (std::size_t corr = 0; corr < pressure_correctors; ++corr) {
        result = correct_pressure(pressure_relaxation);
        if (!result.converged) break;
    }
    return {velocity_rms_change(before), continuity_l2(), continuity_max(), result, pressure_correctors};
}

CollocatedIterationInfo CollocatedIncompressible::iterate_simple() {
    old_velocity_ = velocity_;
    return coupled_sequence(old_velocity_, 1U, config_.pressure_relaxation);
}

CollocatedIterationInfo CollocatedIncompressible::solve_simple(std::size_t max_iterations,
                                                               double velocity_tolerance,
                                                               double continuity_tolerance) {
    if (max_iterations == 0U || !(velocity_tolerance > 0.0) || !(continuity_tolerance > 0.0)) {
        throw std::invalid_argument("invalid SIMPLE convergence controls");
    }
    CollocatedIterationInfo info{};
    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
        info = iterate_simple();
        if (!info.pressure.converged) return info;
        if (info.velocity_rms_change <= velocity_tolerance && info.continuity_l2 <= continuity_tolerance) return info;
    }
    return info;
}

CollocatedIterationInfo CollocatedIncompressible::step_piso() {
    old_velocity_ = velocity_;
    auto info = coupled_sequence(old_velocity_, config_.pressure_correctors, 1.0);
    time_ += config_.dt;
    ++steps_;
    return info;
}

CollocatedIterationInfo CollocatedIncompressible::step_pimple() {
    old_velocity_ = velocity_;
    CollocatedIterationInfo info{};
    for (std::size_t outer = 0; outer < config_.outer_correctors; ++outer) {
        info = coupled_sequence(old_velocity_, config_.pressure_correctors, 1.0);
        if (!info.pressure.converged) break;
    }
    time_ += config_.dt;
    ++steps_;
    return info;
}

double CollocatedIncompressible::continuity_l2() const {
    return flux_divergence_l2(mesh_, face_flux_);
}

double CollocatedIncompressible::continuity_max() const {
    return flux_divergence_max(mesh_, face_flux_);
}

double CollocatedIncompressible::kinetic_energy() const {
    const double sum = cfd::core::parallel_sum(mesh_.cell_count(), [&](std::size_t cell) {
        return 0.5 * dot(velocity_[cell], velocity_[cell]) * mesh_.cells()[cell].volume;
    });
    return config_.density * sum;
}

} // namespace cfd::fvm
