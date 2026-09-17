#include "cfd/solvers/fvm/collocated_incompressible.hpp"

#include "cfd/core/parallel.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"
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
        config_.momentum_sweeps == 0U || config_.momentum_iterations == 0U ||
        !(config_.momentum_tolerance > 0.0) || config_.pressure_iterations == 0U ||
        !(config_.pressure_tolerance > 0.0) || config_.pressure_correctors == 0U ||
        config_.outer_correctors == 0U || !(config_.velocity_relaxation > 0.0 && config_.velocity_relaxation <= 1.0) ||
        !(config_.pressure_relaxation > 0.0 && config_.pressure_relaxation <= 1.0)) {
        throw std::invalid_argument("invalid collocated incompressible solver controls");
    }
    if (config_.temporal_scheme == TemporalScheme::crank_nicolson) {
        (void)crank_nicolson_implicit_weight(config_.crank_nicolson_off_centering);
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

void CollocatedIncompressible::set_pressure_amg_cycle(
    std::function<void(std::span<const double>, std::span<double>)> apply) {
    if (!apply) throw std::invalid_argument("pressure AMG cycle callback is empty");
    pressure_amg_cycle_ = std::move(apply);
    config_.pressure_preconditioner = PressurePreconditionerKind::external_amg;
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

void CollocatedIncompressible::momentum_predictor(std::span<const Vec3> time_source,
                                                   bool use_physical_time_scheme) {
    if (time_source.size() != mesh_.cell_count()) throw std::invalid_argument("time source size mismatch");
    const auto p_boundary = boundary_pressure_values(mesh_, pressure_, pressure_boundary_);
    const auto grad_p = gauss_gradient_scalar(mesh_, pressure_, p_boundary);
    const std::size_t n = mesh_.cell_count();

    TemporalScheme scheme = use_physical_time_scheme ? config_.temporal_scheme : TemporalScheme::euler;
    const bool use_bdf2 = scheme == TemporalScheme::backward_bdf2 && steps_ > 0U;
    const double alpha0 = use_bdf2 ? 1.5 : 1.0;
    const double theta = scheme == TemporalScheme::crank_nicolson
        ? crank_nicolson_implicit_weight(config_.crank_nicolson_off_centering) : 1.0;

    cfd::core::CsrBuilder spatial_builder(n, n);
    std::vector<Vec3> boundary_rhs(n);
    for (std::size_t cell = 0; cell < n; ++cell) {
        double ap = 0.0;
        Vec3 row_boundary{};
        for (const std::size_t fi : mesh_.cell_faces()[cell]) {
            const auto& face = mesh_.faces()[fi];
            const double phi_out = outward_flux(face, cell, face_flux_[fi]);
            const auto geometry = decompose_face_area(mesh_, fi);
            const double diffusion = config_.kinematic_viscosity * geometry.orthogonal_metric;
            if (!face.boundary()) {
                const std::size_t other = face.owner == cell ? face.neighbour : face.owner;
                const double entering = config_.include_convection ? std::max(-phi_out, 0.0) : 0.0;
                const double leaving = config_.include_convection ? std::max(phi_out, 0.0) : 0.0;
                const double coupling = diffusion + entering;
                ap += diffusion + leaving;
                if (coupling != 0.0) spatial_builder.add(cell, other, -coupling);
            } else {
                const auto& bc = velocity_boundary_[face.patch];
                if (bc.type == VelocityBoundaryType::fixedValue) {
                    const Vec3 ub = boundary_velocity_for_momentum(face, time_source[cell], bc);
                    const double entering = config_.include_convection ? std::max(-phi_out, 0.0) : 0.0;
                    const double leaving = config_.include_convection ? std::max(phi_out, 0.0) : 0.0;
                    ap += diffusion + leaving;
                    row_boundary += ub * (diffusion + entering);
                } else if (config_.include_convection) {
                    ap += phi_out;
                }
            }
        }
        if (!std::isfinite(ap)) throw std::runtime_error("non-finite momentum spatial diagonal");
        if (ap != 0.0) spatial_builder.add(cell, cell, ap);
        boundary_rhs[cell] = row_boundary;
    }
    const auto spatial_matrix = spatial_builder.build();

    cfd::core::CsrBuilder solve_builder(n, n);
    const auto& ro = spatial_matrix.row_offsets();
    const auto& ci = spatial_matrix.column_indices();
    const auto& av = spatial_matrix.values();
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t k = ro[row]; k < ro[row + 1U]; ++k) {
            solve_builder.add(row, ci[k], theta * av[k]);
        }
        solve_builder.add(row, row, alpha0 * mesh_.cells()[row].volume / config_.dt);
    }
    const auto matrix = solve_builder.build();
    const auto diagonal = matrix.diagonal();

    std::array<std::vector<double>, 3> old_components{
        std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    std::array<std::vector<double>, 3> spatial_old{
        std::vector<double>(n, 0.0), std::vector<double>(n, 0.0), std::vector<double>(n, 0.0)};
    for (std::size_t cell = 0; cell < n; ++cell) {
        old_components[0][cell] = time_source[cell].x;
        old_components[1][cell] = time_source[cell].y;
        old_components[2][cell] = time_source[cell].z;
    }
    if (theta < 1.0) {
        for (std::size_t component = 0; component < 3U; ++component) {
            spatial_matrix.multiply(old_components[component], spatial_old[component]);
        }
    }

    std::array<std::vector<double>, 3> component_rhs{
        std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    std::array<std::vector<double>, 3> component_x{
        std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    for (std::size_t cell = 0; cell < n; ++cell) {
        const double vdt = mesh_.cells()[cell].volume / config_.dt;
        Vec3 temporal = time_source[cell];
        if (use_bdf2) temporal = time_source[cell] * 2.0 - old_velocity_[cell] * 0.5;
        Vec3 full_rhs = temporal * vdt + boundary_rhs[cell] - grad_p[cell] * mesh_.cells()[cell].volume;
        if (theta < 1.0) {
            full_rhs.x -= (1.0 - theta) * spatial_old[0][cell];
            full_rhs.y -= (1.0 - theta) * spatial_old[1][cell];
            full_rhs.z -= (1.0 - theta) * spatial_old[2][cell];
        }
        component_rhs[0][cell] = full_rhs.x;
        component_rhs[1][cell] = full_rhs.y;
        component_rhs[2][cell] = full_rhs.z;
        component_x[0][cell] = velocity_[cell].x;
        component_x[1][cell] = velocity_[cell].y;
        component_x[2][cell] = velocity_[cell].z;
    }

    if (config_.use_krylov_momentum) {
        const cfd::core::Ilu0Preconditioner preconditioner(matrix);
        for (std::size_t component = 0; component < 3U; ++component) {
            momentum_results_[component] = cfd::core::restarted_gmres(
                component_rhs[component], component_x[component],
                [&](std::span<const double> in, std::span<double> out) { matrix.multiply(in, out); },
                [&](std::span<const double> in, std::span<double> out) { preconditioner(in, out); },
                config_.momentum_iterations, 30U, config_.momentum_tolerance);
            if (!momentum_results_[component].converged) {
                throw std::runtime_error("collocated momentum Krylov solve did not converge");
            }
        }
    } else {
        // Legacy deterministic fixed-point fallback, now driven by the same
        // assembled temporal/spatial equation as the production Krylov path.
        for (std::size_t component = 0; component < 3U; ++component) {
            std::vector<double> next(n, 0.0);
            for (std::size_t sweep = 0; sweep < config_.momentum_sweeps; ++sweep) {
                cfd::core::parallel_for(n, [&](std::size_t row) {
                    double sum = component_rhs[component][row];
                    for (std::size_t k = matrix.row_offsets()[row]; k < matrix.row_offsets()[row + 1U]; ++k) {
                        const std::size_t col = matrix.column_indices()[k];
                        if (col != row) sum -= matrix.values()[k] * component_x[component][col];
                    }
                    if (!(diagonal[row] > 0.0) || !std::isfinite(diagonal[row]))
                        throw std::runtime_error("non-positive momentum diagonal");
                    next[row] = sum / diagonal[row];
                });
                component_x[component].swap(next);
            }
        }
        momentum_results_ = {};
    }

    cfd::core::parallel_for(n, [&](std::size_t cell) {
        velocity_[cell] = {component_x[0][cell], component_x[1][cell], component_x[2][cell]};
        pressure_mobility_[cell] = mesh_.cells()[cell].volume / diagonal[cell];
        h_by_a_[cell] = velocity_[cell] + grad_p[cell] * pressure_mobility_[cell];
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

        if (config_.pressure_preconditioner == PressurePreconditionerKind::none) {
            result = cfd::core::conjugate_gradient(
                std::span<const double>(pressure_rhs_.data(), pressure_rhs_.size()),
                std::span<double>(p_work.data(), p_work.size()),
                [&](std::span<const double> x, std::span<double> out) { apply_pressure_operator(x, out); },
                pressure_workspace_, config_.pressure_iterations, config_.pressure_tolerance);
        } else {
            std::vector<double> pressure_diagonal(mesh_.cell_count(), 0.0);
            cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
                double diagonal = 0.0;
                for (const std::size_t fi : mesh_.cell_faces()[cell]) diagonal += pressure_face_coefficient_[fi];
                pressure_diagonal[cell] = diagonal;
            });
            cfd::core::IterativeSolverResult iterative{};
            if (config_.pressure_preconditioner == PressurePreconditionerKind::jacobi) {
                const cfd::core::JacobiPreconditioner jacobi(pressure_diagonal);
                iterative = cfd::core::preconditioned_conjugate_gradient(
                    pressure_rhs_, p_work,
                    [&](std::span<const double> x, std::span<double> out) { apply_pressure_operator(x, out); },
                    [&](std::span<const double> r, std::span<double> z) { jacobi(r, z); },
                    pressure_krylov_workspace_, config_.pressure_iterations, config_.pressure_tolerance);
            } else {
                if (!pressure_amg_cycle_) throw std::runtime_error("external pressure AMG cycle was not configured");
                iterative = cfd::core::preconditioned_conjugate_gradient(
                    pressure_rhs_, p_work,
                    [&](std::span<const double> x, std::span<double> out) { apply_pressure_operator(x, out); },
                    [&](std::span<const double> r, std::span<double> z) { pressure_amg_cycle_(r, z); },
                    pressure_krylov_workspace_, config_.pressure_iterations, config_.pressure_tolerance);
            }
            result = {iterative.iterations, iterative.residual_rms, iterative.converged};
        }
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
                                                                   double pressure_relaxation,
                                                                   bool use_physical_time_scheme) {
    const std::vector<Vec3> before = velocity_;
    momentum_predictor(time_source, use_physical_time_scheme);
    cfd::core::ConjugateGradientResult result{};
    for (std::size_t corr = 0; corr < pressure_correctors; ++corr) {
        result = correct_pressure(pressure_relaxation);
        if (!result.converged) break;
    }
    return {velocity_rms_change(before), continuity_l2(), continuity_max(), result, pressure_correctors};
}

CollocatedIterationInfo CollocatedIncompressible::iterate_simple() {
    const std::vector<Vec3> pseudo_time_source = velocity_;
    return coupled_sequence(pseudo_time_source, 1U, config_.pressure_relaxation, false);
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
    const std::vector<Vec3> current_time_level = velocity_;
    auto info = coupled_sequence(current_time_level, config_.pressure_correctors, 1.0, true);
    old_velocity_ = current_time_level;
    time_ += config_.dt;
    ++steps_;
    return info;
}

CollocatedIterationInfo CollocatedIncompressible::step_pimple() {
    const std::vector<Vec3> current_time_level = velocity_;
    CollocatedIterationInfo info{};
    for (std::size_t outer = 0; outer < config_.outer_correctors; ++outer) {
        info = coupled_sequence(current_time_level, config_.pressure_correctors, 1.0, true);
        if (!info.pressure.converged) break;
    }
    old_velocity_ = current_time_level;
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
