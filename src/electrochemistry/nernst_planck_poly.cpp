#include "cfd/solvers/electrochemistry/nernst_planck_poly.hpp"

#include "cfd/chemistry/kinetics.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"
#include "cfd/fvm/pressure_velocity.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfd::electrochemistry {
namespace {

[[nodiscard]] double bernoulli_exponential_fit(double x) noexcept {
    const double ax = std::abs(x);
    if (ax < 1.0e-6) {
        // B(x)=x/(exp(x)-1), series avoids cancellation around zero.
        return 1.0 - 0.5 * x + x * x / 12.0;
    }
    if (x > 50.0) return x * std::exp(-x);
    if (x < -50.0) return -x;
    return x / std::expm1(x);
}

} // namespace

NernstPlanckPolyMesh::NernstPlanckPolyMesh(cfd::fvm::PolyMesh mesh, NernstPlanckPolyConfig config)
    : mesh_(std::move(mesh)),
      config_(config),
      potential_(mesh_.cell_count(), 0.0),
      potential_boundary_(mesh_.patches().size()),
      face_flux_(mesh_.face_count(), 0.0) {
    if (!(config_.temperature > 0.0) || !(config_.dt > 0.0) || !(config_.relative_permittivity > 0.0)) {
        throw std::invalid_argument("invalid PolyMesh Nernst-Planck controls");
    }
}

std::size_t NernstPlanckPolyMesh::add_species(std::string name,
                                               int charge,
                                               double diffusivity,
                                               double initial_concentration) {
    return add_species(std::move(name), charge, diffusivity,
                       [initial_concentration](cfd::fvm::Vec3) { return initial_concentration; });
}

std::size_t NernstPlanckPolyMesh::add_species(
    std::string name,
    int charge,
    double diffusivity,
    const std::function<double(cfd::fvm::Vec3)>& initial_concentration) {
    if (name.empty() || !(diffusivity >= 0.0)) throw std::invalid_argument("invalid transport species");
    PolyTransportSpecies species;
    species.name = std::move(name);
    species.charge = charge;
    species.diffusivity = diffusivity;
    species.concentration.resize(mesh_.cell_count());
    species.boundary.resize(mesh_.patches().size());
    for (std::size_t c = 0; c < mesh_.cell_count(); ++c) {
        species.concentration[c] = initial_concentration(mesh_.cells()[c].center);
        if (!(species.concentration[c] >= 0.0) || !std::isfinite(species.concentration[c])) {
            throw std::invalid_argument("initial species concentration must be finite and non-negative");
        }
    }
    species_.push_back(std::move(species));
    return species_.size() - 1U;
}

void NernstPlanckPolyMesh::set_species_boundary(std::size_t species,
                                                 std::string_view patch,
                                                 SpeciesBoundaryType type,
                                                 double value) {
    if (species >= species_.size()) throw std::out_of_range("species index out of range");
    if (type == SpeciesBoundaryType::fixedConcentration && (!(value >= 0.0) || !std::isfinite(value))) {
        throw std::invalid_argument("fixed concentration must be finite and non-negative");
    }
    if (!std::isfinite(value)) throw std::invalid_argument("species boundary value must be finite");
    auto& bc = species_[species].boundary[mesh_.patch_index(patch)];
    bc = {};
    bc.type = type;
    bc.value = value;
}

void NernstPlanckPolyMesh::set_species_butler_volmer_boundary(
    std::size_t species,
    std::string_view patch,
    double electrode_potential,
    double equilibrium_potential,
    double exchange_current_density,
    int electrons,
    double stoichiometric,
    double alpha_anodic,
    double alpha_cathodic) {
    if (species >= species_.size()) throw std::out_of_range("species index out of range");
    if (!(exchange_current_density >= 0.0) || !std::isfinite(exchange_current_density) || electrons <= 0 ||
        !std::isfinite(electrode_potential) || !std::isfinite(equilibrium_potential) ||
        !std::isfinite(stoichiometric) || !(alpha_anodic > 0.0) || !(alpha_cathodic > 0.0)) {
        throw std::invalid_argument("invalid Butler-Volmer boundary parameters");
    }
    auto& bc = species_[species].boundary[mesh_.patch_index(patch)];
    bc = {};
    bc.type = SpeciesBoundaryType::butlerVolmerFlux;
    bc.electrode_potential = electrode_potential;
    bc.equilibrium_potential = equilibrium_potential;
    bc.exchange_current_density = exchange_current_density;
    bc.alpha_anodic = alpha_anodic;
    bc.alpha_cathodic = alpha_cathodic;
    bc.electrons = electrons;
    bc.stoichiometric = stoichiometric;
}

void NernstPlanckPolyMesh::set_potential(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("potential must be finite");
    std::fill(potential_.begin(), potential_.end(), value);
}

void NernstPlanckPolyMesh::set_potential(const std::function<double(cfd::fvm::Vec3)>& value) {
    for (std::size_t c = 0; c < mesh_.cell_count(); ++c) {
        potential_[c] = value(mesh_.cells()[c].center);
        if (!std::isfinite(potential_[c])) throw std::invalid_argument("potential must be finite");
    }
}


void NernstPlanckPolyMesh::set_potential_boundary(std::string_view patch,
                                                   PotentialBoundaryType type,
                                                   double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("potential boundary value must be finite");
    potential_boundary_[mesh_.patch_index(patch)] = {type, value};
}

void NernstPlanckPolyMesh::set_face_flux(std::vector<double> volumetric_flux) {
    if (volumetric_flux.size() != mesh_.face_count()) throw std::invalid_argument("Nernst-Planck face-flux size mismatch");
    face_flux_ = std::move(volumetric_flux);
}

void NernstPlanckPolyMesh::step(std::size_t steps) {
    for (std::size_t i = 0; i < steps; ++i) step_once();
}


void NernstPlanckPolyMesh::step_poisson_nernst_planck(
    std::size_t steps,
    std::size_t poisson_max_iterations,
    double poisson_relative_tolerance) {
    for (std::size_t i = 0; i < steps; ++i) {
        solve_poisson_potential(poisson_max_iterations, poisson_relative_tolerance);
        step_once();
    }
}

void NernstPlanckPolyMesh::step_once() {
    const double migration_scale = faraday_constant / (cfd::chemistry::gas_constant * config_.temperature);
    for (auto& species : species_) {
        std::vector<double> flux(mesh_.face_count(), 0.0); // oriented owner -> neighbour/outside
        cfd::core::parallel_for(mesh_.face_count(), [&](std::size_t fi) {
            const auto& face = mesh_.faces()[fi];
            if (face.boundary()) {
                const auto& bc = species.boundary[face.patch];
                const double area = cfd::fvm::magnitude(face.area);
                if (bc.type == SpeciesBoundaryType::noFlux) {
                    flux[fi] = 0.0;
                    return;
                }
                if (bc.type == SpeciesBoundaryType::fixedMolarFlux) {
                    flux[fi] = bc.value * area;
                    return;
                }
                if (bc.type == SpeciesBoundaryType::butlerVolmerFlux) {
                    const auto& pbc = potential_boundary_[face.patch];
                    const double electrolyte_potential = pbc.type == PotentialBoundaryType::fixedPotential
                        ? pbc.value : potential_[face.owner];
                    const double eta = bc.electrode_potential - electrolyte_potential - bc.equilibrium_potential;
                    const double current = butler_volmer_current_density(
                        bc.exchange_current_density, eta, config_.temperature,
                        static_cast<double>(bc.electrons), bc.alpha_anodic, bc.alpha_cathodic);
                    // Positive anodic current creates a positive-stoichiometry
                    // product inside the electrolyte, i.e. negative outward flux.
                    const double molar_flux_out = -bc.stoichiometric * current
                        / (static_cast<double>(bc.electrons) * faraday_constant);
                    flux[fi] = molar_flux_out * area;
                    return;
                }

                const double co = species.concentration[face.owner];
                const double cb = bc.value;
                const auto geometry = cfd::fvm::decompose_face_area(mesh_, fi);
                const double advective = face_flux_[fi] >= 0.0 ? face_flux_[fi] * co : face_flux_[fi] * cb;
                const auto& pbc = potential_boundary_[face.patch];
                const double dphi = pbc.type == PotentialBoundaryType::fixedPotential
                    ? pbc.value - potential_[face.owner] : 0.0;
                if (config_.electromigration_scheme == ElectromigrationFluxScheme::scharfetterGummel
                    && species.charge != 0) {
                    const double psi = static_cast<double>(species.charge) * migration_scale * dphi;
                    const double drift_diffusion = species.diffusivity * geometry.orthogonal_metric
                        * (co * bernoulli_exponential_fit(psi) - cb * bernoulli_exponential_fit(-psi));
                    flux[fi] = advective + drift_diffusion;
                } else {
                    const double diffusive = -species.diffusivity * geometry.orthogonal_metric * (cb - co);
                    const double migration = -static_cast<double>(species.charge) * species.diffusivity * migration_scale
                        * 0.5 * (co + cb) * geometry.orthogonal_metric * dphi;
                    flux[fi] = advective + diffusive + migration;
                }
                return;
            }

            const std::size_t owner = face.owner;
            const std::size_t neighbour = face.neighbour;
            const double co = species.concentration[owner];
            const double cn = species.concentration[neighbour];
            const double phi = face_flux_[fi];
            const double cf = phi >= 0.0 ? co : cn;
            const auto geometry = cfd::fvm::decompose_face_area(mesh_, fi);
            const double dc = cn - co;
            const double dphi = potential_[neighbour] - potential_[owner];
            if (config_.electromigration_scheme == ElectromigrationFluxScheme::scharfetterGummel
                && species.charge != 0) {
                const double psi = static_cast<double>(species.charge) * migration_scale * dphi;
                const double drift_diffusion = species.diffusivity * geometry.orthogonal_metric
                    * (co * bernoulli_exponential_fit(psi) - cn * bernoulli_exponential_fit(-psi));
                flux[fi] = phi * cf + drift_diffusion;
            } else {
                const double diffusion = -species.diffusivity * geometry.orthogonal_metric * dc;
                const double migration = -static_cast<double>(species.charge) * species.diffusivity * migration_scale
                    * (0.5 * (co + cn)) * geometry.orthogonal_metric * dphi;
                flux[fi] = phi * cf + diffusion + migration;
            }
        });

        std::vector<double> next(species.concentration.size(), 0.0);
        cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
            double net_out = 0.0;
            for (const std::size_t fi : mesh_.cell_faces()[cell]) {
                const auto& face = mesh_.faces()[fi];
                net_out += face.owner == cell ? flux[fi] : -flux[fi];
            }
            const double updated = species.concentration[cell]
                - config_.dt * net_out / mesh_.cells()[cell].volume;
            if (!std::isfinite(updated) || updated < -1.0e-12) {
                throw std::runtime_error("Nernst-Planck explicit step produced a negative/invalid concentration; reduce dt");
            }
            next[cell] = std::max(0.0, updated);
        });
        species.concentration.swap(next);
    }
    time_ += config_.dt;
}

cfd::core::ConjugateGradientResult NernstPlanckPolyMesh::solve_electroneutral_potential(
    std::size_t max_iterations,
    double relative_tolerance) {
    if (species_.empty()) throw std::runtime_error("electroneutral potential requires at least one ionic species");
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) throw std::invalid_argument("invalid potential solve controls");

    const double conductivity_scale = faraday_constant * faraday_constant
        / (cfd::chemistry::gas_constant * config_.temperature);
    std::vector<double> coefficient(mesh_.face_count(), 0.0);
    std::vector<double> diffusion_current(mesh_.face_count(), 0.0);
    cfd::core::parallel_for(mesh_.face_count(), [&](std::size_t fi) {
        const auto& face = mesh_.faces()[fi];
        if (face.boundary()) return; // insulating/current-free baseline boundary
        const auto geometry = cfd::fvm::decompose_face_area(mesh_, fi);
        double conductivity = 0.0;
        double diffusive = 0.0;
        for (const auto& species : species_) {
            if (species.charge == 0 || species.diffusivity == 0.0) continue;
            const double co = species.concentration[face.owner];
            const double cn = species.concentration[face.neighbour];
            const double cf = 0.5 * (co + cn);
            const double z = static_cast<double>(species.charge);
            conductivity += z * z * species.diffusivity * cf;
            diffusive += z * species.diffusivity * (cn - co);
        }
        coefficient[fi] = conductivity_scale * conductivity * geometry.orthogonal_metric;
        diffusion_current[fi] = -faraday_constant * diffusive * geometry.orthogonal_metric;
    });

    std::vector<double> rhs(mesh_.cell_count(), 0.0);
    cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
        double value = 0.0;
        for (const std::size_t fi : mesh_.cell_faces()[cell]) {
            const auto& face = mesh_.faces()[fi];
            if (face.boundary()) continue;
            const double outward_diffusion = face.owner == cell ? diffusion_current[fi] : -diffusion_current[fi];
            value -= outward_diffusion;
        }
        rhs[cell] = value;
    });
    double rhs_mean = cfd::core::parallel_sum(rhs.size(), [&](std::size_t i) { return rhs[i]; });
    rhs_mean /= static_cast<double>(rhs.size());
    cfd::core::parallel_for(rhs.size(), [&](std::size_t i) { rhs[i] -= rhs_mean; });
    double phi_mean = cfd::core::parallel_sum(potential_.size(), [&](std::size_t i) { return potential_[i]; });
    phi_mean /= static_cast<double>(potential_.size());
    cfd::core::parallel_for(potential_.size(), [&](std::size_t i) { potential_[i] -= phi_mean; });

    potential_result_ = cfd::core::conjugate_gradient(
        rhs, potential_,
        [&](std::span<const double> x, std::span<double> out) {
            cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
                double sum = 0.0;
                for (const std::size_t fi : mesh_.cell_faces()[cell]) {
                    const auto& face = mesh_.faces()[fi];
                    if (face.boundary()) continue;
                    const std::size_t other = face.owner == cell ? face.neighbour : face.owner;
                    sum += coefficient[fi] * (x[cell] - x[other]);
                }
                out[cell] = sum;
            });
        },
        potential_workspace_, max_iterations, relative_tolerance);

    if (!potential_result_.converged) throw std::runtime_error("electroneutral potential solve did not converge");
    phi_mean = cfd::core::parallel_sum(potential_.size(), [&](std::size_t i) { return potential_[i]; });
    phi_mean /= static_cast<double>(potential_.size());
    cfd::core::parallel_for(potential_.size(), [&](std::size_t i) { potential_[i] -= phi_mean; });
    return potential_result_;
}

cfd::core::ConjugateGradientResult NernstPlanckPolyMesh::solve_poisson_potential(
    std::size_t max_iterations,
    double relative_tolerance) {
    if (species_.empty()) throw std::runtime_error("Poisson potential requires at least one charged species");
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) throw std::invalid_argument("invalid Poisson solve controls");

    constexpr double vacuum_permittivity = 8.8541878128e-12;
    const double epsilon = vacuum_permittivity * config_.relative_permittivity;
    std::vector<double> rhs(mesh_.cell_count(), 0.0);
    cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
        double charge_density = 0.0;
        for (const auto& sp : species_) {
            charge_density += static_cast<double>(sp.charge) * sp.concentration[cell];
        }
        // Scale the Poisson equation by 1/epsilon. This leaves the solution
        // unchanged while keeping the matrix/RHS near engineering magnitudes
        // so the shared dimensional CG tolerance does not see ~1e-12 SI
        // coefficients.
        rhs[cell] = (faraday_constant / epsilon) * charge_density * mesh_.cells()[cell].volume;
    });

    bool has_dirichlet = false;
    for (const auto& bc : potential_boundary_) {
        has_dirichlet = has_dirichlet || bc.type == PotentialBoundaryType::fixedPotential;
    }

    if (!has_dirichlet) {
        const double net_charge = cfd::core::parallel_sum(rhs.size(), [&](std::size_t i) { return rhs[i]; });
        const double scale = cfd::core::parallel_sum(rhs.size(), [&](std::size_t i) { return std::abs(rhs[i]); });
        if (std::abs(net_charge) > 1.0e-10 * std::max(1.0, scale)) {
            throw std::runtime_error("pure-insulating Poisson problem requires net-zero charge");
        }
        const double mean = net_charge / static_cast<double>(rhs.size());
        cfd::core::parallel_for(rhs.size(), [&](std::size_t i) { rhs[i] -= mean; });
        double phi_mean = cfd::core::parallel_sum(potential_.size(), [&](std::size_t i) { return potential_[i]; });
        phi_mean /= static_cast<double>(potential_.size());
        cfd::core::parallel_for(potential_.size(), [&](std::size_t i) { potential_[i] -= phi_mean; });
    } else {
        for (std::size_t fi = 0; fi < mesh_.face_count(); ++fi) {
            const auto& face = mesh_.faces()[fi];
            if (!face.boundary()) continue;
            const auto& bc = potential_boundary_[face.patch];
            if (bc.type != PotentialBoundaryType::fixedPotential) continue;
            const double coeff = cfd::fvm::decompose_face_area(mesh_, fi).orthogonal_metric;
            rhs[face.owner] += coeff * bc.value;
        }
    }

    potential_result_ = cfd::core::conjugate_gradient(
        rhs, potential_,
        [&](std::span<const double> x, std::span<double> out) {
            cfd::core::parallel_for(mesh_.cell_count(), [&](std::size_t cell) {
                double sum = 0.0;
                for (const std::size_t fi : mesh_.cell_faces()[cell]) {
                    const auto& face = mesh_.faces()[fi];
                    const double coeff = cfd::fvm::decompose_face_area(mesh_, fi).orthogonal_metric;
                    if (!face.boundary()) {
                        const std::size_t other = face.owner == cell ? face.neighbour : face.owner;
                        sum += coeff * (x[cell] - x[other]);
                    } else if (potential_boundary_[face.patch].type == PotentialBoundaryType::fixedPotential) {
                        sum += coeff * x[cell];
                    }
                }
                out[cell] = sum;
            });
        },
        potential_workspace_, max_iterations, relative_tolerance);

    if (!potential_result_.converged) throw std::runtime_error("Poisson potential solve did not converge");
    if (!has_dirichlet) {
        double phi_mean = cfd::core::parallel_sum(potential_.size(), [&](std::size_t i) { return potential_[i]; });
        phi_mean /= static_cast<double>(potential_.size());
        cfd::core::parallel_for(potential_.size(), [&](std::size_t i) { potential_[i] -= phi_mean; });
    }
    return potential_result_;
}

double NernstPlanckPolyMesh::total_charge() const {
    return faraday_constant * cfd::core::parallel_sum(mesh_.cell_count(), [&](std::size_t cell) {
        double molar_charge = 0.0;
        for (const auto& sp : species_) molar_charge += static_cast<double>(sp.charge) * sp.concentration[cell];
        return molar_charge * mesh_.cells()[cell].volume;
    });
}

std::vector<double> NernstPlanckPolyMesh::current_flux_faces() const {
    const double migration_scale = faraday_constant * faraday_constant
        / (cfd::chemistry::gas_constant * config_.temperature);
    std::vector<double> current(mesh_.face_count(), 0.0);
    cfd::core::parallel_for(mesh_.face_count(), [&](std::size_t fi) {
        const auto& face = mesh_.faces()[fi];
        if (face.boundary()) return;
        const auto geometry = cfd::fvm::decompose_face_area(mesh_, fi);
        double diffusion = 0.0;
        double conductivity = 0.0;
        for (const auto& species : species_) {
            if (species.charge == 0 || species.diffusivity == 0.0) continue;
            const double co = species.concentration[face.owner];
            const double cn = species.concentration[face.neighbour];
            const double z = static_cast<double>(species.charge);
            diffusion += z * species.diffusivity * (cn - co);
            conductivity += z * z * species.diffusivity * 0.5 * (co + cn);
        }
        current[fi] = -faraday_constant * diffusion * geometry.orthogonal_metric
            - migration_scale * conductivity * geometry.orthogonal_metric
                * (potential_[face.neighbour] - potential_[face.owner]);
    });
    return current;
}

double NernstPlanckPolyMesh::total_amount(std::size_t species) const {
    if (species >= species_.size()) throw std::out_of_range("species index out of range");
    return cfd::core::parallel_sum(mesh_.cell_count(), [&](std::size_t c) {
        return species_[species].concentration[c] * mesh_.cells()[c].volume;
    });
}

double NernstPlanckPolyMesh::minimum_concentration(std::size_t species) const {
    if (species >= species_.size()) throw std::out_of_range("species index out of range");
    return *std::min_element(species_[species].concentration.begin(), species_[species].concentration.end());
}

} // namespace cfd::electrochemistry
