#include "cfd/solvers/fvm/rans_transport.hpp"

#include "cfd/core/csr_matrix.hpp"
#include "cfd/fvm/advanced_models.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/fvm/pressure_velocity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace cfd::fvm {
namespace {

[[nodiscard]] double outward_flux(const Face& face, std::size_t cell, double flux) noexcept {
    return face.owner == cell ? flux : -flux;
}

void validate_controls(const RansTransportControls& c) {
    if (!(c.dt > 0.0) || c.linear_iterations == 0U || c.gmres_restart == 0U || !(c.linear_tolerance > 0.0)) {
        throw std::invalid_argument("invalid RANS transport controls");
    }
    if (c.temporal_scheme == TemporalScheme::crank_nicolson) {
        (void)crank_nicolson_implicit_weight(c.crank_nicolson_off_centering);
    }
}

void validate_nonnegative_field(std::span<const double> values, std::size_t size, const char* what, bool strictly_positive = false) {
    if (values.size() != size) throw std::invalid_argument(std::string(what) + " size mismatch");
    for (const double value : values) {
        if (!std::isfinite(value) || (strictly_positive ? !(value > 0.0) : value < 0.0)) {
            throw std::invalid_argument(std::string("invalid ") + what + " value");
        }
    }
}

void validate_flux(std::span<const double> values, std::size_t size) {
    if (values.size() != size) throw std::invalid_argument("RANS face-flux size mismatch");
    for (const double value : values) if (!std::isfinite(value)) throw std::invalid_argument("non-finite RANS face flux");
}

[[nodiscard]] double harmonic_mean(double a, double b) noexcept {
    if (!(a > 0.0) || !(b > 0.0)) return std::max(0.0, 0.5 * (a + b));
    return 2.0 * a * b / (a + b);
}

[[nodiscard]] cfd::core::IterativeSolverResult advance_scalar(
    const PolyMesh& mesh,
    const RansTransportControls& controls,
    std::span<const TurbulenceScalarBoundaryCondition> boundary,
    std::span<const double> face_flux,
    std::span<const double> diffusivity,
    std::span<const double> source,
    std::span<const double> sink,
    std::vector<double>& values,
    std::vector<double>& previous_values,
    std::size_t steps,
    double minimum_value) {

    const std::size_t n = mesh.cell_count();
    if (boundary.size() != mesh.patches().size() || face_flux.size() != mesh.face_count()
        || diffusivity.size() != n || source.size() != n || sink.size() != n || values.size() != n || previous_values.size() != n) {
        throw std::invalid_argument("RANS scalar transport size mismatch");
    }

    cfd::core::CsrBuilder spatial_builder(n, n);
    std::vector<double> boundary_rhs(n, 0.0);
    for (std::size_t cell = 0; cell < n; ++cell) {
        if (!std::isfinite(diffusivity[cell]) || diffusivity[cell] < 0.0 || !std::isfinite(source[cell])
            || !std::isfinite(sink[cell]) || sink[cell] < 0.0) {
            throw std::runtime_error("invalid RANS equation coefficient");
        }
        double diagonal = sink[cell] * mesh.cells()[cell].volume;
        double rhs = 0.0;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            const auto& face = mesh.faces()[fi];
            const double phi = outward_flux(face, cell, face_flux[fi]);
            if (!face.boundary()) {
                const std::size_t other = face.owner == cell ? face.neighbour : face.owner;
                const double gamma_face = harmonic_mean(diffusivity[cell], diffusivity[other]);
                const double diffusion = gamma_face * decompose_face_area(mesh, fi).orthogonal_metric;
                diagonal += diffusion + std::max(phi, 0.0);
                const double off_diagonal = -diffusion + std::min(phi, 0.0);
                if (off_diagonal != 0.0) spatial_builder.add(cell, other, off_diagonal);
            } else {
                const auto& bc = boundary[face.patch];
                if (bc.type == TurbulenceScalarBoundaryType::fixedValue) {
                    if (!std::isfinite(bc.value) || bc.value < minimum_value) {
                        throw std::runtime_error("invalid fixed turbulence boundary value");
                    }
                    const double distance = magnitude(face.center - mesh.cells()[cell].center);
                    if (!(distance > 0.0)) throw std::runtime_error("degenerate RANS boundary distance");
                    const double diffusion = diffusivity[cell] * magnitude(face.area) / distance;
                    diagonal += diffusion + std::max(phi, 0.0);
                    rhs += diffusion * bc.value - std::min(phi, 0.0) * bc.value;
                } else {
                    // zeroGradient: the boundary value equals the owner value.
                    diagonal += phi;
                }
            }
        }
        if (!std::isfinite(diagonal)) throw std::runtime_error("non-finite RANS transport diagonal");
        if (diagonal != 0.0) spatial_builder.add(cell, cell, diagonal);
        boundary_rhs[cell] = rhs;
    }
    const auto spatial = spatial_builder.build();

    double alpha0 = 1.0;
    double theta = 1.0;
    const bool use_bdf2 = controls.temporal_scheme == TemporalScheme::backward_bdf2 && steps > 0U;
    if (use_bdf2) alpha0 = 1.5;
    if (controls.temporal_scheme == TemporalScheme::crank_nicolson) {
        theta = crank_nicolson_implicit_weight(controls.crank_nicolson_off_centering);
    }

    cfd::core::CsrBuilder matrix_builder(n, n);
    const auto& row_offsets = spatial.row_offsets();
    const auto& columns = spatial.column_indices();
    const auto& coefficients = spatial.values();
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t p = row_offsets[row]; p < row_offsets[row + 1U]; ++p) {
            matrix_builder.add(row, columns[p], theta * coefficients[p]);
        }
        matrix_builder.add(row, row, alpha0 * mesh.cells()[row].volume / controls.dt);
    }
    const auto matrix = matrix_builder.build();
    cfd::core::Ilu0Preconditioner preconditioner(matrix);

    std::vector<double> explicit_spatial(n, 0.0), rhs(n, 0.0), candidate(values);
    if (theta < 1.0) spatial.multiply(values, explicit_spatial);
    for (std::size_t cell = 0; cell < n; ++cell) {
        const double volume_dt = mesh.cells()[cell].volume / controls.dt;
        const double temporal = use_bdf2
            ? (2.0 * values[cell] - 0.5 * previous_values[cell]) * volume_dt
            : values[cell] * volume_dt;
        rhs[cell] = temporal + source[cell] * mesh.cells()[cell].volume + boundary_rhs[cell]
            - (1.0 - theta) * explicit_spatial[cell];
    }

    const auto result = cfd::core::restarted_gmres(
        rhs,
        candidate,
        [&](std::span<const double> input, std::span<double> output) { matrix.multiply(input, output); },
        [&](std::span<const double> input, std::span<double> output) { preconditioner(input, output); },
        controls.linear_iterations,
        controls.gmres_restart,
        controls.linear_tolerance);
    if (!result.converged) throw std::runtime_error("RANS scalar transport linear solve did not converge");

    previous_values = values;
    for (double& value : candidate) {
        if (!std::isfinite(value)) throw std::runtime_error("RANS transport produced a non-finite field value");
        value = std::max(value, minimum_value);
    }
    values.swap(candidate);
    return result;
}

[[nodiscard]] std::vector<double> bounded_kepsilon_nut(
    std::span<const double> k,
    std::span<const double> epsilon,
    const KEpsilonConfig& config) {
    std::vector<double> result(k.size(), 0.0);
    const double maximum = config.maximum_eddy_viscosity_ratio * config.molecular_viscosity;
    for (std::size_t i = 0; i < k.size(); ++i) {
        const double dynamic = k_epsilon_eddy_viscosity(config.density, k[i], epsilon[i], config.c_mu);
        result[i] = std::min(dynamic / config.density, maximum);
    }
    return result;
}

[[nodiscard]] double sst_f2(double k, double omega, double wall_distance, const KOmegaSSTConfig& c) noexcept {
    const double d = std::max(wall_distance, c.minimum_wall_distance);
    const double w = std::max(omega, c.minimum_omega);
    const double arg_a = 2.0 * std::sqrt(std::max(k, c.minimum_k)) / (c.beta_star * w * d);
    const double arg_b = 500.0 * c.molecular_viscosity / (d * d * w);
    const double arg = std::clamp(std::max(arg_a, arg_b), 0.0, 100.0);
    return std::tanh(arg * arg);
}

} // namespace

std::vector<double> turbulence_strain_rate_magnitude(const PolyMesh& mesh, std::span<const Vec3> velocity) {
    if (velocity.size() != mesh.cell_count()) throw std::invalid_argument("turbulence velocity size mismatch");
    std::vector<double> x(velocity.size()), y(velocity.size()), z(velocity.size());
    for (std::size_t i = 0; i < velocity.size(); ++i) {
        x[i] = velocity[i].x;
        y[i] = velocity[i].y;
        z[i] = velocity[i].z;
    }
    const auto gx = least_squares_gradient_scalar(mesh, x);
    const auto gy = least_squares_gradient_scalar(mesh, y);
    const auto gz = least_squares_gradient_scalar(mesh, z);
    std::vector<double> strain(velocity.size(), 0.0);
    for (std::size_t n = 0; n < velocity.size(); ++n) {
        const double gradient[3][3] = {
            {gx[n].x, gx[n].y, gx[n].z},
            {gy[n].x, gy[n].y, gy[n].z},
            {gz[n].x, gz[n].y, gz[n].z}};
        double ss = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                const double sij = 0.5 * (gradient[i][j] + gradient[j][i]);
                ss += sij * sij;
            }
        }
        strain[n] = std::sqrt(2.0 * ss);
    }
    return strain;
}

SpalartAllmarasTransport::SpalartAllmarasTransport(PolyMesh mesh, SpalartAllmarasConfig config)
    : mesh_(std::move(mesh)), config_(config), boundary_(mesh_.patches().size()),
      nu_tilde_(mesh_.cell_count(), config.minimum_nu_tilde), previous_nu_tilde_(nu_tilde_),
      face_flux_(mesh_.face_count(), 0.0), wall_distance_(mesh_.cell_count(), 1.0),
      strain_rate_(mesh_.cell_count(), 0.0) {
    validate_controls(config_.transport);
    if (!(config_.molecular_viscosity > 0.0) || !(config_.cb1 > 0.0) || !(config_.cb2 >= 0.0)
        || !(config_.cw2 >= 0.0) || !(config_.cw3 > 0.0) || !(config_.cv1 > 0.0) || !(config_.cs >= 0.0)
        || !(config_.sigma_nu_tilde > 0.0) || !(config_.kappa > 0.0) || !(config_.minimum_nu_tilde >= 0.0)
        || !(config_.minimum_wall_distance > 0.0)) {
        throw std::invalid_argument("invalid Spalart-Allmaras configuration");
    }
}

void SpalartAllmarasTransport::initialize(double value) {
    if (!std::isfinite(value) || value < config_.minimum_nu_tilde) throw std::invalid_argument("invalid nu_tilde initialization");
    std::fill(nu_tilde_.begin(), nu_tilde_.end(), value);
    previous_nu_tilde_ = nu_tilde_;
    steps_ = 0U;
    time_ = 0.0;
}

void SpalartAllmarasTransport::initialize(std::span<const double> values) {
    validate_nonnegative_field(values, mesh_.cell_count(), "nu_tilde");
    nu_tilde_.assign(values.begin(), values.end());
    for (double& value : nu_tilde_) value = std::max(value, config_.minimum_nu_tilde);
    previous_nu_tilde_ = nu_tilde_;
    steps_ = 0U;
    time_ = 0.0;
}

void SpalartAllmarasTransport::set_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value) {
    if (!std::isfinite(value) || (type == TurbulenceScalarBoundaryType::fixedValue && value < config_.minimum_nu_tilde)) {
        throw std::invalid_argument("invalid Spalart-Allmaras boundary value");
    }
    boundary_[mesh_.patch_index(patch)] = {type, value};
}

void SpalartAllmarasTransport::set_face_flux(std::vector<double> flux) {
    validate_flux(flux, mesh_.face_count());
    face_flux_ = std::move(flux);
}

void SpalartAllmarasTransport::set_wall_distance(double distance) {
    if (!std::isfinite(distance) || distance < config_.minimum_wall_distance) throw std::invalid_argument("invalid SA wall distance");
    std::fill(wall_distance_.begin(), wall_distance_.end(), distance);
}

void SpalartAllmarasTransport::set_wall_distance(std::span<const double> distance) {
    validate_nonnegative_field(distance, mesh_.cell_count(), "SA wall distance", true);
    wall_distance_.assign(distance.begin(), distance.end());
    for (double& value : wall_distance_) value = std::max(value, config_.minimum_wall_distance);
}

void SpalartAllmarasTransport::set_strain_rate(double strain) {
    if (!std::isfinite(strain) || strain < 0.0) throw std::invalid_argument("invalid SA strain rate");
    std::fill(strain_rate_.begin(), strain_rate_.end(), strain);
}

void SpalartAllmarasTransport::set_strain_rate(std::span<const double> strain) {
    validate_nonnegative_field(strain, mesh_.cell_count(), "SA strain rate");
    strain_rate_.assign(strain.begin(), strain.end());
}

void SpalartAllmarasTransport::set_velocity(std::span<const Vec3> velocity) {
    strain_rate_ = turbulence_strain_rate_magnitude(mesh_, velocity);
}

std::vector<double> SpalartAllmarasTransport::kinematic_eddy_viscosity() const {
    std::vector<double> result(mesh_.cell_count());
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = spalart_allmaras_eddy_viscosity(config_.molecular_viscosity, nu_tilde_[i], config_.cv1);
    }
    return result;
}

cfd::core::IterativeSolverResult SpalartAllmarasTransport::step() {
    const auto gradient = least_squares_gradient_scalar(mesh_, nu_tilde_);
    std::vector<double> diffusivity(mesh_.cell_count()), source(mesh_.cell_count()), sink(mesh_.cell_count());
    const double cw1 = config_.cb1 / (config_.kappa * config_.kappa)
        + (1.0 + config_.cb2) / config_.sigma_nu_tilde;
    const double cw3_6 = std::pow(config_.cw3, 6.0);
    for (std::size_t i = 0; i < mesh_.cell_count(); ++i) {
        const double nt = std::max(nu_tilde_[i], config_.minimum_nu_tilde);
        const double d = std::max(wall_distance_[i], config_.minimum_wall_distance);
        const double chi = nt / config_.molecular_viscosity;
        const double chi3 = chi * chi * chi;
        const double cv3 = config_.cv1 * config_.cv1 * config_.cv1;
        const double fv1 = chi3 / (chi3 + cv3);
        const double fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const double kappa2d2 = config_.kappa * config_.kappa * d * d;
        const double raw_stilda = strain_rate_[i] + nt * fv2 / kappa2d2;
        const double stilda = std::max(raw_stilda, config_.cs * strain_rate_[i]);
        const double denominator = std::max(stilda * kappa2d2, 1.0e-30);
        const double r = std::clamp(nt / denominator, 0.0, 10.0);
        const double r6 = std::pow(r, 6.0);
        const double g = r + config_.cw2 * (r6 - r);
        const double g6 = std::pow(g, 6.0);
        const double fw = g * std::pow((1.0 + cw3_6) / std::max(g6 + cw3_6, 1.0e-30), 1.0 / 6.0);
        const double grad2 = dot(gradient[i], gradient[i]);
        diffusivity[i] = (config_.molecular_viscosity + nt) / config_.sigma_nu_tilde;
        source[i] = config_.cb1 * stilda * nt + (config_.cb2 / config_.sigma_nu_tilde) * grad2;
        sink[i] = cw1 * fw * nt / (d * d);
    }
    linear_result_ = advance_scalar(mesh_, config_.transport, boundary_, face_flux_, diffusivity, source, sink,
                                    nu_tilde_, previous_nu_tilde_, steps_, config_.minimum_nu_tilde);
    ++steps_;
    time_ += config_.transport.dt;
    return linear_result_;
}

cfd::core::IterativeSolverResult SpalartAllmarasTransport::run(std::size_t steps) {
    for (std::size_t i = 0; i < steps; ++i) step();
    return linear_result_;
}

KEpsilonTransport::KEpsilonTransport(PolyMesh mesh, KEpsilonConfig config)
    : mesh_(std::move(mesh)), config_(config), k_boundary_(mesh_.patches().size()), epsilon_boundary_(mesh_.patches().size()),
      k_(mesh_.cell_count(), config.minimum_k), epsilon_(mesh_.cell_count(), config.minimum_epsilon),
      previous_k_(k_), previous_epsilon_(epsilon_), face_flux_(mesh_.face_count(), 0.0), strain_rate_(mesh_.cell_count(), 0.0) {
    validate_controls(config_.transport);
    if (!(config_.density > 0.0) || !(config_.molecular_viscosity > 0.0) || !(config_.c_mu > 0.0)
        || !(config_.c1 > 0.0) || !(config_.c2 > 0.0) || !(config_.sigma_k > 0.0) || !(config_.sigma_epsilon > 0.0)
        || !(config_.production_limiter > 0.0) || !(config_.minimum_k > 0.0) || !(config_.minimum_epsilon > 0.0)
        || !(config_.maximum_eddy_viscosity_ratio > 0.0)) {
        throw std::invalid_argument("invalid k-epsilon configuration");
    }
}

void KEpsilonTransport::initialize(double k, double epsilon) {
    if (!std::isfinite(k) || !std::isfinite(epsilon) || k < config_.minimum_k || epsilon < config_.minimum_epsilon) {
        throw std::invalid_argument("invalid k-epsilon initialization");
    }
    std::fill(k_.begin(), k_.end(), k);
    std::fill(epsilon_.begin(), epsilon_.end(), epsilon);
    previous_k_ = k_;
    previous_epsilon_ = epsilon_;
    steps_ = 0U;
    time_ = 0.0;
}

void KEpsilonTransport::initialize(std::span<const double> k, std::span<const double> epsilon) {
    validate_nonnegative_field(k, mesh_.cell_count(), "k", true);
    validate_nonnegative_field(epsilon, mesh_.cell_count(), "epsilon", true);
    k_.assign(k.begin(), k.end());
    epsilon_.assign(epsilon.begin(), epsilon.end());
    for (double& value : k_) value = std::max(value, config_.minimum_k);
    for (double& value : epsilon_) value = std::max(value, config_.minimum_epsilon);
    previous_k_ = k_;
    previous_epsilon_ = epsilon_;
    steps_ = 0U;
    time_ = 0.0;
}

void KEpsilonTransport::set_k_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value) {
    if (!std::isfinite(value) || (type == TurbulenceScalarBoundaryType::fixedValue && value < config_.minimum_k)) {
        throw std::invalid_argument("invalid k boundary value");
    }
    k_boundary_[mesh_.patch_index(patch)] = {type, value};
}

void KEpsilonTransport::set_epsilon_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value) {
    if (!std::isfinite(value) || (type == TurbulenceScalarBoundaryType::fixedValue && value < config_.minimum_epsilon)) {
        throw std::invalid_argument("invalid epsilon boundary value");
    }
    epsilon_boundary_[mesh_.patch_index(patch)] = {type, value};
}

void KEpsilonTransport::set_face_flux(std::vector<double> flux) {
    validate_flux(flux, mesh_.face_count());
    face_flux_ = std::move(flux);
}

void KEpsilonTransport::set_strain_rate(double strain) {
    if (!std::isfinite(strain) || strain < 0.0) throw std::invalid_argument("invalid k-epsilon strain rate");
    std::fill(strain_rate_.begin(), strain_rate_.end(), strain);
}

void KEpsilonTransport::set_strain_rate(std::span<const double> strain) {
    validate_nonnegative_field(strain, mesh_.cell_count(), "k-epsilon strain rate");
    strain_rate_.assign(strain.begin(), strain.end());
}

void KEpsilonTransport::set_velocity(std::span<const Vec3> velocity) {
    strain_rate_ = turbulence_strain_rate_magnitude(mesh_, velocity);
}

std::vector<double> KEpsilonTransport::kinematic_eddy_viscosity() const {
    return bounded_kepsilon_nut(k_, epsilon_, config_);
}

std::vector<double> KEpsilonTransport::dynamic_eddy_viscosity() const {
    auto result = kinematic_eddy_viscosity();
    for (double& value : result) value *= config_.density;
    return result;
}

TwoEquationLinearResult KEpsilonTransport::step() {
    const auto nut = bounded_kepsilon_nut(k_, epsilon_, config_);
    const std::size_t n = mesh_.cell_count();
    std::vector<double> gamma_k(n), gamma_epsilon(n), source_k(n), source_epsilon(n), sink_k(n), sink_epsilon(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double k = std::max(k_[i], config_.minimum_k);
        const double epsilon = std::max(epsilon_[i], config_.minimum_epsilon);
        const double production = std::min(nut[i] * strain_rate_[i] * strain_rate_[i], config_.production_limiter * epsilon);
        gamma_k[i] = config_.molecular_viscosity + nut[i] / config_.sigma_k;
        gamma_epsilon[i] = config_.molecular_viscosity + nut[i] / config_.sigma_epsilon;
        source_k[i] = production;
        sink_k[i] = epsilon / k;
        source_epsilon[i] = config_.c1 * (epsilon / k) * production;
        sink_epsilon[i] = config_.c2 * epsilon / k;
    }
    linear_result_.first = advance_scalar(mesh_, config_.transport, k_boundary_, face_flux_, gamma_k, source_k, sink_k,
                                          k_, previous_k_, steps_, config_.minimum_k);
    linear_result_.second = advance_scalar(mesh_, config_.transport, epsilon_boundary_, face_flux_, gamma_epsilon,
                                           source_epsilon, sink_epsilon, epsilon_, previous_epsilon_, steps_, config_.minimum_epsilon);
    ++steps_;
    time_ += config_.transport.dt;
    return linear_result_;
}

TwoEquationLinearResult KEpsilonTransport::run(std::size_t steps) {
    for (std::size_t i = 0; i < steps; ++i) step();
    return linear_result_;
}

KOmegaSSTTransport::KOmegaSSTTransport(PolyMesh mesh, KOmegaSSTConfig config)
    : mesh_(std::move(mesh)), config_(config), k_boundary_(mesh_.patches().size()), omega_boundary_(mesh_.patches().size()),
      k_(mesh_.cell_count(), config.minimum_k), omega_(mesh_.cell_count(), config.minimum_omega),
      previous_k_(k_), previous_omega_(omega_), face_flux_(mesh_.face_count(), 0.0), wall_distance_(mesh_.cell_count(), 1.0),
      grid_scale_(mesh_.cell_count(), 0.0), strain_rate_(mesh_.cell_count(), 0.0), blending_f1_(mesh_.cell_count(), 0.0),
      blending_f2_(mesh_.cell_count(), 0.0), hybrid_factor_(mesh_.cell_count(), 1.0) {
    for (std::size_t i = 0; i < mesh_.cell_count(); ++i) grid_scale_[i] = std::cbrt(mesh_.cells()[i].volume);
    validate_controls(config_.transport);
    if (!(config_.density > 0.0) || !(config_.molecular_viscosity > 0.0) || !(config_.beta_star > 0.0)
        || !(config_.beta1 > 0.0) || !(config_.beta2 > 0.0) || !(config_.gamma1 > 0.0) || !(config_.gamma2 > 0.0)
        || !(config_.alpha_k1 > 0.0) || !(config_.alpha_k2 > 0.0) || !(config_.alpha_omega1 > 0.0)
        || !(config_.alpha_omega2 > 0.0) || !(config_.a1 > 0.0) || !(config_.production_limiter > 0.0)
        || !(config_.minimum_k > 0.0) || !(config_.minimum_omega > 0.0) || !(config_.minimum_wall_distance > 0.0)
        || !(config_.maximum_eddy_viscosity_ratio > 0.0) || !(config_.c_des > 0.0)) {
        throw std::invalid_argument("invalid k-omega SST configuration");
    }
    update_blending();
}

void KOmegaSSTTransport::initialize(double k, double omega) {
    if (!std::isfinite(k) || !std::isfinite(omega) || k < config_.minimum_k || omega < config_.minimum_omega) {
        throw std::invalid_argument("invalid k-omega SST initialization");
    }
    std::fill(k_.begin(), k_.end(), k);
    std::fill(omega_.begin(), omega_.end(), omega);
    previous_k_ = k_;
    previous_omega_ = omega_;
    steps_ = 0U;
    time_ = 0.0;
    update_blending();
}

void KOmegaSSTTransport::initialize(std::span<const double> k, std::span<const double> omega) {
    validate_nonnegative_field(k, mesh_.cell_count(), "SST k", true);
    validate_nonnegative_field(omega, mesh_.cell_count(), "SST omega", true);
    k_.assign(k.begin(), k.end());
    omega_.assign(omega.begin(), omega.end());
    for (double& value : k_) value = std::max(value, config_.minimum_k);
    for (double& value : omega_) value = std::max(value, config_.minimum_omega);
    previous_k_ = k_;
    previous_omega_ = omega_;
    steps_ = 0U;
    time_ = 0.0;
    update_blending();
}

void KOmegaSSTTransport::set_k_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value) {
    if (!std::isfinite(value) || (type == TurbulenceScalarBoundaryType::fixedValue && value < config_.minimum_k)) {
        throw std::invalid_argument("invalid SST k boundary value");
    }
    k_boundary_[mesh_.patch_index(patch)] = {type, value};
}

void KOmegaSSTTransport::set_omega_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value) {
    if (!std::isfinite(value) || (type == TurbulenceScalarBoundaryType::fixedValue && value < config_.minimum_omega)) {
        throw std::invalid_argument("invalid SST omega boundary value");
    }
    omega_boundary_[mesh_.patch_index(patch)] = {type, value};
}

void KOmegaSSTTransport::set_face_flux(std::vector<double> flux) {
    validate_flux(flux, mesh_.face_count());
    face_flux_ = std::move(flux);
}

void KOmegaSSTTransport::set_wall_distance(double distance) {
    if (!std::isfinite(distance) || distance < config_.minimum_wall_distance) throw std::invalid_argument("invalid SST wall distance");
    std::fill(wall_distance_.begin(), wall_distance_.end(), distance);
    update_blending();
}

void KOmegaSSTTransport::set_wall_distance(std::span<const double> distance) {
    validate_nonnegative_field(distance, mesh_.cell_count(), "SST wall distance", true);
    wall_distance_.assign(distance.begin(), distance.end());
    for (double& value : wall_distance_) value = std::max(value, config_.minimum_wall_distance);
    update_blending();
}

void KOmegaSSTTransport::set_grid_scale(double scale) {
    if (!std::isfinite(scale) || !(scale > 0.0)) throw std::invalid_argument("invalid SST-DES grid scale");
    std::fill(grid_scale_.begin(), grid_scale_.end(), scale);
    update_blending();
}

void KOmegaSSTTransport::set_grid_scale(std::span<const double> scale) {
    validate_nonnegative_field(scale, mesh_.cell_count(), "SST-DES grid scale", true);
    grid_scale_.assign(scale.begin(), scale.end());
    update_blending();
}

void KOmegaSSTTransport::set_strain_rate(double strain) {
    if (!std::isfinite(strain) || strain < 0.0) throw std::invalid_argument("invalid SST strain rate");
    std::fill(strain_rate_.begin(), strain_rate_.end(), strain);
    update_blending();
}

void KOmegaSSTTransport::set_strain_rate(std::span<const double> strain) {
    validate_nonnegative_field(strain, mesh_.cell_count(), "SST strain rate");
    strain_rate_.assign(strain.begin(), strain.end());
    update_blending();
}

void KOmegaSSTTransport::set_velocity(std::span<const Vec3> velocity) {
    strain_rate_ = turbulence_strain_rate_magnitude(mesh_, velocity);
    update_blending();
}

void KOmegaSSTTransport::update_blending() {
    if (mesh_.cell_count() == 0U) return;
    const auto grad_k = least_squares_gradient_scalar(mesh_, k_);
    const auto grad_omega = least_squares_gradient_scalar(mesh_, omega_);
    for (std::size_t i = 0; i < mesh_.cell_count(); ++i) {
        const double k = std::max(k_[i], config_.minimum_k);
        const double omega = std::max(omega_[i], config_.minimum_omega);
        const double d = std::max(wall_distance_[i], config_.minimum_wall_distance);
        const double cross = dot(grad_k[i], grad_omega[i]);
        const double cd_kw = std::max(2.0 * config_.alpha_omega2 * cross / omega, 1.0e-20);
        const double arg_a = std::sqrt(k) / (config_.beta_star * omega * d);
        const double arg_b = 500.0 * config_.molecular_viscosity / (d * d * omega);
        const double arg_c = 4.0 * config_.alpha_omega2 * k / (cd_kw * d * d);
        const double arg1 = std::clamp(std::min(std::max(arg_a, arg_b), arg_c), 0.0, 100.0);
        const double arg1_2 = arg1 * arg1;
        blending_f1_[i] = std::tanh(arg1_2 * arg1_2);
        blending_f2_[i] = sst_f2(k, omega, d, config_);
        hybrid_factor_[i] = 1.0;
        if (config_.des_enabled) {
            const double turbulent_length = std::sqrt(k) / (config_.beta_star * omega);
            const double les_length = config_.c_des * std::max(grid_scale_[i], 1.0e-30);
            const double raw = std::max(turbulent_length / les_length, 1.0);
            double shield = 1.0;
            if (config_.des_zonal_filter == SSTDESZonalFilter::f1) shield = 1.0 - blending_f1_[i];
            if (config_.des_zonal_filter == SSTDESZonalFilter::f2) shield = 1.0 - blending_f2_[i];
            hybrid_factor_[i] = 1.0 + (raw - 1.0) * std::clamp(shield, 0.0, 1.0);
        }
    }
}

std::vector<double> KOmegaSSTTransport::kinematic_eddy_viscosity() const {
    std::vector<double> result(mesh_.cell_count());
    const double maximum = config_.maximum_eddy_viscosity_ratio * config_.molecular_viscosity;
    for (std::size_t i = 0; i < result.size(); ++i) {
        const double dynamic = k_omega_sst_eddy_viscosity(config_.density, k_[i], omega_[i], strain_rate_[i], blending_f2_[i], config_.a1);
        result[i] = std::min(dynamic / config_.density, maximum);
    }
    return result;
}

std::vector<double> KOmegaSSTTransport::dynamic_eddy_viscosity() const {
    auto result = kinematic_eddy_viscosity();
    for (double& value : result) value *= config_.density;
    return result;
}

TwoEquationLinearResult KOmegaSSTTransport::step() {
    update_blending();
    const auto nut = kinematic_eddy_viscosity();
    const auto grad_k = least_squares_gradient_scalar(mesh_, k_);
    const auto grad_omega = least_squares_gradient_scalar(mesh_, omega_);
    const std::size_t n = mesh_.cell_count();
    std::vector<double> gamma_k(n), gamma_omega(n), source_k(n), source_omega(n), sink_k(n), sink_omega(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double k = std::max(k_[i], config_.minimum_k);
        const double omega = std::max(omega_[i], config_.minimum_omega);
        const double f1 = blending_f1_[i];
        const double alpha_k = f1 * config_.alpha_k1 + (1.0 - f1) * config_.alpha_k2;
        const double alpha_omega = f1 * config_.alpha_omega1 + (1.0 - f1) * config_.alpha_omega2;
        const double beta = f1 * config_.beta1 + (1.0 - f1) * config_.beta2;
        const double gamma = f1 * config_.gamma1 + (1.0 - f1) * config_.gamma2;
        const double production = std::min(
            nut[i] * strain_rate_[i] * strain_rate_[i],
            config_.production_limiter * config_.beta_star * k * omega);
        const double production_ratio = nut[i] > 1.0e-30 ? production / nut[i] : 0.0;
        const double cross = 2.0 * (1.0 - f1) * config_.alpha_omega2
            * dot(grad_k[i], grad_omega[i]) / omega;
        gamma_k[i] = config_.molecular_viscosity + alpha_k * nut[i];
        gamma_omega[i] = config_.molecular_viscosity + alpha_omega * nut[i];
        source_k[i] = production;
        sink_k[i] = config_.beta_star * omega * hybrid_factor_[i];
        source_omega[i] = gamma * production_ratio + std::max(cross, 0.0);
        sink_omega[i] = beta * omega + (cross < 0.0 ? -cross / omega : 0.0);
    }
    linear_result_.first = advance_scalar(mesh_, config_.transport, k_boundary_, face_flux_, gamma_k, source_k, sink_k,
                                          k_, previous_k_, steps_, config_.minimum_k);
    linear_result_.second = advance_scalar(mesh_, config_.transport, omega_boundary_, face_flux_, gamma_omega,
                                           source_omega, sink_omega, omega_, previous_omega_, steps_, config_.minimum_omega);
    ++steps_;
    time_ += config_.transport.dt;
    update_blending();
    return linear_result_;
}

TwoEquationLinearResult KOmegaSSTTransport::run(std::size_t steps) {
    for (std::size_t i = 0; i < steps; ++i) step();
    return linear_result_;
}

double turbulent_kinetic_energy(const SymmetricTensor3& stress) noexcept {
    return 0.5 * (stress.xx + stress.yy + stress.zz);
}

bool reynolds_stress_is_realizable(const SymmetricTensor3& r, double tolerance) noexcept {
    if (!std::isfinite(r.xx) || !std::isfinite(r.yy) || !std::isfinite(r.zz) || !std::isfinite(r.xy)
        || !std::isfinite(r.xz) || !std::isfinite(r.yz) || tolerance < 0.0) return false;
    if (r.xx < -tolerance || r.yy < -tolerance || r.zz < -tolerance) return false;
    if (r.xx * r.yy - r.xy * r.xy < -tolerance) return false;
    if (r.xx * r.zz - r.xz * r.xz < -tolerance) return false;
    if (r.yy * r.zz - r.yz * r.yz < -tolerance) return false;
    const double determinant = r.xx * r.yy * r.zz + 2.0 * r.xy * r.xz * r.yz
        - r.xx * r.yz * r.yz - r.yy * r.xz * r.xz - r.zz * r.xy * r.xy;
    return determinant >= -tolerance;
}

namespace {

[[nodiscard]] SymmetricTensor3 stress_from_components(const std::array<std::vector<double>, 6>& fields, std::size_t i) {
    return {fields[0][i], fields[1][i], fields[2][i], fields[3][i], fields[4][i], fields[5][i]};
}

void set_stress_components(std::array<std::vector<double>, 6>& fields, std::size_t i, const SymmetricTensor3& r) {
    fields[0][i] = r.xx; fields[1][i] = r.yy; fields[2][i] = r.zz;
    fields[3][i] = r.xy; fields[4][i] = r.xz; fields[5][i] = r.yz;
}

[[nodiscard]] SymmetricTensor3 production_tensor(const SymmetricTensor3& r, const VelocityGradient3& g) noexcept {
    const double rm[3][3] = {{r.xx, r.xy, r.xz}, {r.xy, r.yy, r.yz}, {r.xz, r.yz, r.zz}};
    double p[3][3]{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double value = 0.0;
            for (int k = 0; k < 3; ++k) {
                value -= rm[i][k] * g[static_cast<std::size_t>(j * 3 + k)]
                    + rm[j][k] * g[static_cast<std::size_t>(i * 3 + k)];
            }
            p[i][j] = value;
        }
    }
    return {p[0][0], p[1][1], p[2][2], 0.5 * (p[0][1] + p[1][0]),
            0.5 * (p[0][2] + p[2][0]), 0.5 * (p[1][2] + p[2][1])};
}

[[nodiscard]] SymmetricTensor3 lrr_source(const SymmetricTensor3& r, const SymmetricTensor3& p,
                                           double epsilon, const ReynoldsStressConfig& c) noexcept {
    const double k = std::max(turbulent_kinetic_energy(r), c.minimum_k);
    const double two_thirds_k = 2.0 * k / 3.0;
    const double p_trace_third = (p.xx + p.yy + p.zz) / 3.0;
    const double slow = c.pressure_strain_slow * epsilon / k;
    SymmetricTensor3 out{};
    out.xx = p.xx - slow * (r.xx - two_thirds_k) - c.pressure_strain_rapid * (p.xx - p_trace_third) - 2.0 * epsilon / 3.0;
    out.yy = p.yy - slow * (r.yy - two_thirds_k) - c.pressure_strain_rapid * (p.yy - p_trace_third) - 2.0 * epsilon / 3.0;
    out.zz = p.zz - slow * (r.zz - two_thirds_k) - c.pressure_strain_rapid * (p.zz - p_trace_third) - 2.0 * epsilon / 3.0;
    out.xy = p.xy - slow * r.xy - c.pressure_strain_rapid * p.xy;
    out.xz = p.xz - slow * r.xz - c.pressure_strain_rapid * p.xz;
    out.yz = p.yz - slow * r.yz - c.pressure_strain_rapid * p.yz;
    return out;
}

[[nodiscard]] SymmetricTensor3 project_realizable(SymmetricTensor3 r, double normal_floor) noexcept {
    r.xx = std::max(r.xx, normal_floor);
    r.yy = std::max(r.yy, normal_floor);
    r.zz = std::max(r.zz, normal_floor);
    const double margin = 1.0 - 1.0e-12;
    r.xy = std::clamp(r.xy, -margin * std::sqrt(r.xx * r.yy), margin * std::sqrt(r.xx * r.yy));
    r.xz = std::clamp(r.xz, -margin * std::sqrt(r.xx * r.zz), margin * std::sqrt(r.xx * r.zz));
    r.yz = std::clamp(r.yz, -margin * std::sqrt(r.yy * r.zz), margin * std::sqrt(r.yy * r.zz));
    if (reynolds_stress_is_realizable(r, 1.0e-14)) return r;
    const double xy = r.xy, xz = r.xz, yz = r.yz;
    double low = 0.0, high = 1.0;
    for (int iteration = 0; iteration < 64; ++iteration) {
        const double mid = 0.5 * (low + high);
        SymmetricTensor3 candidate = r;
        candidate.xy = mid * xy; candidate.xz = mid * xz; candidate.yz = mid * yz;
        if (reynolds_stress_is_realizable(candidate, 1.0e-14)) low = mid; else high = mid;
    }
    r.xy = low * xy; r.xz = low * xz; r.yz = low * yz;
    return r;
}

} // namespace

ReynoldsStressTransport::ReynoldsStressTransport(PolyMesh mesh, ReynoldsStressConfig config)
    : mesh_(std::move(mesh)), config_(config), epsilon_boundary_(mesh_.patches().size()),
      epsilon_(mesh_.cell_count(), config.minimum_epsilon), previous_epsilon_(epsilon_),
      face_flux_(mesh_.face_count(), 0.0), velocity_gradient_(mesh_.cell_count()) {
    validate_controls(config_.transport);
    if (!(config_.density > 0.0) || !(config_.molecular_viscosity > 0.0) || !(config_.c_mu > 0.0)
        || !(config_.pressure_strain_slow >= 0.0) || !(config_.pressure_strain_rapid >= 0.0)
        || !(config_.c_epsilon1 > 0.0) || !(config_.c_epsilon2 > 0.0) || !(config_.sigma_r > 0.0)
        || !(config_.sigma_epsilon > 0.0) || !(config_.minimum_k > 0.0) || !(config_.minimum_epsilon > 0.0)
        || !(config_.minimum_normal_stress > 0.0) || !(config_.maximum_eddy_viscosity_ratio > 0.0)) {
        throw std::invalid_argument("invalid Reynolds-stress configuration");
    }
    for (std::size_t component = 0; component < 6U; ++component) {
        stress_boundary_[component].assign(mesh_.patches().size(), {});
        stress_[component].assign(mesh_.cell_count(), component < 3U ? 2.0 * config_.minimum_k / 3.0 : 0.0);
        previous_stress_[component] = stress_[component];
    }
}

void ReynoldsStressTransport::initialize_isotropic(double k, double epsilon) {
    if (!std::isfinite(k) || !std::isfinite(epsilon) || k < config_.minimum_k || epsilon < config_.minimum_epsilon) {
        throw std::invalid_argument("invalid isotropic Reynolds-stress initialization");
    }
    const double normal = 2.0 * k / 3.0;
    for (std::size_t component = 0; component < 6U; ++component) {
        std::fill(stress_[component].begin(), stress_[component].end(), component < 3U ? normal : 0.0);
        previous_stress_[component] = stress_[component];
    }
    std::fill(epsilon_.begin(), epsilon_.end(), epsilon);
    previous_epsilon_ = epsilon_;
    steps_ = 0U;
    time_ = 0.0;
}

void ReynoldsStressTransport::initialize(std::span<const SymmetricTensor3> stress, std::span<const double> epsilon) {
    if (stress.size() != mesh_.cell_count()) throw std::invalid_argument("Reynolds-stress field size mismatch");
    validate_nonnegative_field(epsilon, mesh_.cell_count(), "RSM epsilon", true);
    for (std::size_t i = 0; i < stress.size(); ++i) {
        if (!reynolds_stress_is_realizable(stress[i]) || cfd::fvm::turbulent_kinetic_energy(stress[i]) < config_.minimum_k) {
            throw std::invalid_argument("non-realizable Reynolds-stress initialization");
        }
        set_stress_components(stress_, i, stress[i]);
    }
    for (std::size_t component = 0; component < 6U; ++component) previous_stress_[component] = stress_[component];
    epsilon_.assign(epsilon.begin(), epsilon.end());
    for (double& value : epsilon_) value = std::max(value, config_.minimum_epsilon);
    previous_epsilon_ = epsilon_;
    steps_ = 0U;
    time_ = 0.0;
}

void ReynoldsStressTransport::set_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, SymmetricTensor3 value) {
    if (type == TurbulenceScalarBoundaryType::fixedValue && !reynolds_stress_is_realizable(value)) {
        throw std::invalid_argument("non-realizable Reynolds-stress boundary value");
    }
    const std::size_t index = mesh_.patch_index(patch);
    const double values[6] = {value.xx, value.yy, value.zz, value.xy, value.xz, value.yz};
    for (std::size_t component = 0; component < 6U; ++component) {
        stress_boundary_[component][index] = {type, values[component]};
    }
}

void ReynoldsStressTransport::set_epsilon_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value) {
    if (!std::isfinite(value) || (type == TurbulenceScalarBoundaryType::fixedValue && value < config_.minimum_epsilon)) {
        throw std::invalid_argument("invalid RSM epsilon boundary value");
    }
    epsilon_boundary_[mesh_.patch_index(patch)] = {type, value};
}

void ReynoldsStressTransport::set_face_flux(std::vector<double> flux) {
    validate_flux(flux, mesh_.face_count());
    face_flux_ = std::move(flux);
}

void ReynoldsStressTransport::set_velocity(std::span<const Vec3> velocity) {
    if (velocity.size() != mesh_.cell_count()) throw std::invalid_argument("RSM velocity size mismatch");
    std::vector<double> u(velocity.size()), v(velocity.size()), w(velocity.size());
    for (std::size_t i = 0; i < velocity.size(); ++i) { u[i] = velocity[i].x; v[i] = velocity[i].y; w[i] = velocity[i].z; }
    const auto gu = least_squares_gradient_scalar(mesh_, u);
    const auto gv = least_squares_gradient_scalar(mesh_, v);
    const auto gw = least_squares_gradient_scalar(mesh_, w);
    for (std::size_t i = 0; i < velocity.size(); ++i) {
        velocity_gradient_[i] = {gu[i].x, gu[i].y, gu[i].z, gv[i].x, gv[i].y, gv[i].z, gw[i].x, gw[i].y, gw[i].z};
    }
}

void ReynoldsStressTransport::set_velocity_gradient(std::span<const VelocityGradient3> gradient) {
    if (gradient.size() != mesh_.cell_count()) throw std::invalid_argument("RSM velocity-gradient size mismatch");
    for (const auto& tensor : gradient) for (const double value : tensor) if (!std::isfinite(value)) throw std::invalid_argument("invalid RSM velocity gradient");
    velocity_gradient_.assign(gradient.begin(), gradient.end());
}

void ReynoldsStressTransport::set_additional_source(std::function<SymmetricTensor3(Vec3, double)> source) {
    additional_source_ = std::move(source);
}

std::vector<SymmetricTensor3> ReynoldsStressTransport::reynolds_stress() const {
    std::vector<SymmetricTensor3> result(mesh_.cell_count());
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = stress_from_components(stress_, i);
    return result;
}

std::vector<double> ReynoldsStressTransport::turbulent_kinetic_energy() const {
    std::vector<double> result(mesh_.cell_count());
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = cfd::fvm::turbulent_kinetic_energy(stress_from_components(stress_, i));
    return result;
}

std::vector<double> ReynoldsStressTransport::kinematic_eddy_viscosity() const {
    const auto k = turbulent_kinetic_energy();
    std::vector<double> result(k.size());
    const double maximum = config_.maximum_eddy_viscosity_ratio * config_.molecular_viscosity;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = std::min(config_.c_mu * k[i] * k[i] / std::max(epsilon_[i], config_.minimum_epsilon), maximum);
    }
    return result;
}

void ReynoldsStressTransport::enforce_realizability() {
    for (std::size_t i = 0; i < mesh_.cell_count(); ++i) {
        const auto projected = project_realizable(stress_from_components(stress_, i), config_.minimum_normal_stress);
        set_stress_components(stress_, i, projected);
    }
}

ReynoldsStressLinearResult ReynoldsStressTransport::step() {
    const std::size_t n = mesh_.cell_count();
    const auto k = turbulent_kinetic_energy();
    const auto nut = kinematic_eddy_viscosity();
    std::vector<double> gamma_r(n), gamma_epsilon(n), epsilon_source(n), epsilon_sink(n);
    std::array<std::vector<double>, 6> source;
    for (auto& component : source) component.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto r = stress_from_components(stress_, i);
        const auto p = production_tensor(r, velocity_gradient_[i]);
        auto total_source = lrr_source(r, p, epsilon_[i], config_);
        if (additional_source_) {
            const auto extra = additional_source_(mesh_.cells()[i].center, time_);
            total_source.xx += extra.xx; total_source.yy += extra.yy; total_source.zz += extra.zz;
            total_source.xy += extra.xy; total_source.xz += extra.xz; total_source.yz += extra.yz;
        }
        const double values[6] = {total_source.xx, total_source.yy, total_source.zz, total_source.xy, total_source.xz, total_source.yz};
        for (std::size_t component = 0; component < 6U; ++component) source[component][i] = values[component];
        gamma_r[i] = config_.molecular_viscosity + nut[i] / config_.sigma_r;
        gamma_epsilon[i] = config_.molecular_viscosity + nut[i] / config_.sigma_epsilon;
        const double ki = std::max(k[i], config_.minimum_k);
        const double eps = std::max(epsilon_[i], config_.minimum_epsilon);
        const double pk = std::max(0.0, 0.5 * (p.xx + p.yy + p.zz));
        epsilon_source[i] = config_.c_epsilon1 * (eps / ki) * pk;
        epsilon_sink[i] = config_.c_epsilon2 * eps / ki;
    }
    std::vector<double> zero_sink(n, 0.0);
    for (std::size_t component = 0; component < 6U; ++component) {
        const double floor = component < 3U ? config_.minimum_normal_stress : -std::numeric_limits<double>::max();
        linear_result_.stress[component] = advance_scalar(mesh_, config_.transport, stress_boundary_[component], face_flux_,
            gamma_r, source[component], zero_sink, stress_[component], previous_stress_[component], steps_, floor);
    }
    linear_result_.epsilon = advance_scalar(mesh_, config_.transport, epsilon_boundary_, face_flux_, gamma_epsilon,
        epsilon_source, epsilon_sink, epsilon_, previous_epsilon_, steps_, config_.minimum_epsilon);
    enforce_realizability();
    ++steps_;
    time_ += config_.transport.dt;
    return linear_result_;
}

ReynoldsStressLinearResult ReynoldsStressTransport::run(std::size_t steps) {
    for (std::size_t i = 0; i < steps; ++i) step();
    return linear_result_;
}

} // namespace cfd::fvm
