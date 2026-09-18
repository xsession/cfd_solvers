#pragma once

#include "cfd/battery/lithium_ion.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace cfd::battery {

// Full Doyle-Fuller-Newman / pseudo-two-dimensional (P2D) porous-electrode cell
// model. This is the distributed limit of the single-particle models: the solid
// particle concentration, the electrolyte concentration and the solid/electrolyte
// potentials are resolved over the electrode and separator thickness instead of
// being collapsed to electrode-averaged states.
//
// Geometry is one-dimensional along the through-cell axis x in [0, L] with
// L = L_neg + L_sep + L_pos. The domain is divided into one electrolyte node per
// electrode node plus the separator. Each electrode node owns one
// SphericalDiffusionParticle whose surface reaction flux is set from a local
// Butler-Volmer-consistent partition of the imposed terminal current.
struct DfnConfig {
    double area_m2{0.010};
    double negative_thickness_m{80.0e-6};
    double separator_thickness_m{25.0e-6};
    double positive_thickness_m{75.0e-6};
    ParticleConfig negative_particle{};
    ParticleConfig positive_particle{};
    double negative_solid_conductivity_s_per_m{100.0};
    double positive_solid_conductivity_s_per_m{10.0};
    double negative_exchange_current_density_a_per_m2{4.0};
    double positive_exchange_current_density_a_per_m2{4.0};
    double negative_porosity{0.30};
    double positive_porosity{0.30};
    double separator_porosity{0.50};
    double bruggeman_exponent{1.5};
    double negative_active_volume_fraction{0.62};
    double positive_active_volume_fraction{0.58};
    double electrolyte_initial_concentration_mol_per_m3{1000.0};
    double electrolyte_diffusivity_m2_per_s{2.0e-10};
    double electrolyte_conductivity_s_per_m{1.0};
    double transference_number{0.38};
    std::function<double(double)> negative_ocv_v{};
    std::function<double(double)> positive_ocv_v{};
    std::size_t negative_nodes{24U};
    std::size_t separator_nodes{9U};
    std::size_t positive_nodes{22U};
    double contact_resistance_ohm{0.012};
    double initial_temperature_k{298.15};
};

[[nodiscard]] DfnConfig default_graphite_nmc_dfn_config();

struct DfnStepResult {
    double time_s{};
    double current_a{}; // positive = discharge
    double terminal_voltage_v{};
    double open_circuit_voltage_v{};
    double reaction_overpotential_v{};
    double ohmic_drop_v{}; // solid + electrolyte + contact
    double negative_surface_stoichiometry{}; // area-averaged over negative nodes
    double positive_surface_stoichiometry{};
    double electrolyte_min_concentration_mol_per_m3{};
    double electrolyte_max_concentration_mol_per_m3{};
    double temperature_k{};
    double heat_generation_w{};
};

class DoyleFullerNewmanModel {
public:
    explicit DoyleFullerNewmanModel(DfnConfig config = {});
    void reset();
    [[nodiscard]] DfnStepResult step(double current_a, double dt_s);
    [[nodiscard]] DfnStepResult state(double current_a = 0.0) const;
    [[nodiscard]] double time_s() const noexcept { return time_s_; }
    [[nodiscard]] double temperature_k() const noexcept { return temperature_k_; }
    [[nodiscard]] const DfnConfig& config() const noexcept { return config_; }

    // Distributed state. Electrolyte fields span the whole cell
    // (negative_nodes + separator_nodes + positive_nodes). Solid potentials and
    // local reaction rates are zero in the separator.
    [[nodiscard]] const std::vector<double>& electrolyte_concentration() const noexcept { return electrolyte_concentration_; }
    [[nodiscard]] const std::vector<double>& electrolyte_potential() const noexcept { return electrolyte_potential_v_; }
    [[nodiscard]] const std::vector<double>& solid_potential() const noexcept { return solid_potential_v_; }
    [[nodiscard]] const std::vector<double>& local_reaction_current() const noexcept { return local_reaction_current_; }
    [[nodiscard]] const std::vector<SphericalDiffusionParticle>& negative_particles() const noexcept { return negative_particles_; }
    [[nodiscard]] const std::vector<SphericalDiffusionParticle>& positive_particles() const noexcept { return positive_particles_; }

    // Total reversible lithium inventory (particles + electrolyte) in mol, used
    // to verify global conservation across the distributed state.
    [[nodiscard]] double total_lithium_mol() const noexcept;

private:
    // Electrode-specific helpers.
    [[nodiscard]] std::size_t total_nodes() const noexcept;
    [[nodiscard]] double cell_length_m() const noexcept;
    [[nodiscard]] double node_spacing_m() const noexcept;
    [[nodiscard]] double region_porosity(double x_m) const;
    [[nodiscard]] double solid_conductivity(double x_m) const;

    void update_reaction_partition(double current_a);
    void advance_particles(double current_a, double dt_s);
    void advance_electrolyte_concentration(double current_a, double dt_s);
    void update_potentials(double current_a);
    [[nodiscard]] DfnStepResult evaluate(double current_a) const;

    DfnConfig config_{};
    std::vector<SphericalDiffusionParticle> negative_particles_{};
    std::vector<SphericalDiffusionParticle> positive_particles_{};
    std::vector<double> electrolyte_concentration_{};
    std::vector<double> electrolyte_potential_v_{};
    std::vector<double> solid_potential_v_{};
    std::vector<double> local_reaction_current_{}; // A/m^2 cross-section, zero in separator
    double temperature_k_{};
    double time_s_{};
};

} // namespace cfd::battery
