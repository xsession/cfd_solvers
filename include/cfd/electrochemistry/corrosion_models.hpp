#pragma once
#include "cfd/chemistry/aqueous_equilibrium.hpp"
#include <complex>
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::electrochemistry {
struct RedoxActivityTerm {
    double activity{};
    double stoichiometric_power{};
};
[[nodiscard]] double reaction_quotient(std::span<const RedoxActivityTerm> products,
                                       std::span<const RedoxActivityTerm> reactants);
[[nodiscard]] double concentration_dependent_nernst(double standard_potential, double temperature, double electrons,
                                                    std::span<const RedoxActivityTerm> products,
                                                    std::span<const RedoxActivityTerm> reactants);
struct CathodicReactionResult {
    double equilibrium_potential{};
    double activation_current_density{};
    double limited_current_density{};
};
[[nodiscard]] CathodicReactionResult oxygen_reduction_cathodic(double electrode_potential, double standard_potential,
                                                               double oxygen_activity, double proton_activity,
                                                               double exchange_current_density,
                                                               double limiting_current_density,
                                                               double temperature = 298.15, double electrons = 4.0);
[[nodiscard]] CathodicReactionResult
hydrogen_evolution_cathodic(double electrode_potential, double hydrogen_pressure_activity, double proton_activity,
                            double exchange_current_density, double limiting_current_density,
                            double temperature = 298.15, double electrons = 2.0);
[[nodiscard]] double charge_transfer_resistance(double exchange_current_density, double temperature = 298.15,
                                                double electrons = 1.0, double anodic_transfer = 0.5,
                                                double cathodic_transfer = 0.5);
[[nodiscard]] std::complex<double> randles_impedance(double frequency_hz, double solution_resistance,
                                                     double charge_transfer_resistance_ohm,
                                                     double double_layer_capacitance_f, double warburg_sigma = 0.0);
struct DoubleLayerState {
    double capacitance_per_area{};
    double potential{};
    [[nodiscard]] double charge_density() const noexcept { return capacitance_per_area * potential; }
    void advance_current(double current_density, double dt);
};

struct PassivationModel {
    double active_equilibrium_potential{};
    double exchange_current_density{};
    double passive_onset{};
    double transpassive_onset{};
    double passive_current_density{};
    double transpassive_slope{}; // A/(m2 V)
    double electrons{1.0};
};
[[nodiscard]] double passivation_current_density(double electrode_potential, const PassivationModel& model,
                                                 double temperature = 298.15);

struct ProductLayerState {
    double thickness{};   // m
    double resistivity{}; // ohm m
    double growth_fraction{1.0};
    [[nodiscard]] double area_specific_resistance() const noexcept { return thickness * resistivity; }
    void advance_from_anodic_current(double current_density, double dt, double molar_mass, double density,
                                     double electrons = 1.0);
};

[[nodiscard]] double faradaic_recession_distance(double current_density, double dt, double molar_mass,
                                                 double metal_density, double electrons = 1.0);
[[nodiscard]] double cathodic_protection_current_density(double protected_potential, double free_corrosion_potential,
                                                         double polarization_resistance);

[[nodiscard]] double bruggeman_effective_transport(double bulk_value, double porosity, double exponent = 1.5);
[[nodiscard]] double solid_phase_ohmic_drop(double current_density, double length, double conductivity);
[[nodiscard]] double marcus_current_density(double exchange_current_density, double overpotential,
                                            double reorganization_energy_j_per_mol, double temperature = 298.15,
                                            double electrons = 1.0);
[[nodiscard]] double bikerman_activity_correction(double concentration, double occupied_volume_per_mole,
                                                  double total_concentration);
struct ModifiedPnpFlux {
    double molar_flux{};
    double activity{};
};
[[nodiscard]] ModifiedPnpFlux modified_nernst_planck_flux(double concentration_left, double concentration_right,
                                                          double diffusivity, double potential_gradient, int charge,
                                                          double temperature, double cell_width,
                                                          double occupied_volume_per_mole, double total_concentration,
                                                          double advective_velocity = 0.0);
struct ImmersedInterfaceWeight {
    double electrolyte_fraction{}, delta{};
};
[[nodiscard]] ImmersedInterfaceWeight smooth_electrode_interface(double signed_distance, double thickness);
struct CutCellInterface {
    double electrolyte_fraction{};
    double electrode_fraction{};
    double interface_area_fraction{};
};
// Conservative 1-D cut-cell geometry contract. Signed distance is negative in
// the electrode; the returned fractions sum to one and the interface-area
// fraction is one only when the cell is actually cut.
[[nodiscard]] CutCellInterface cut_cell_electrode_interface(double left_signed_distance, double right_signed_distance,
                                                            double cell_width);
struct CoupledElectrodeResult {
    double surface_concentration{}, current_density{}, molar_flux{};
    std::size_t iterations{};
};
[[nodiscard]] CoupledElectrodeResult
solve_anodic_bv_mass_transfer(double bulk_concentration, double mass_transfer_coefficient, double standard_potential,
                              double electrode_potential, double exchange_current_density, double temperature = 298.15,
                              double electrons = 2.0, double reference_concentration = 1.0);
[[nodiscard]] double acid_speciation_redox_potential(double standard_potential, double temperature, double electrons,
                                                     double proton_stoichiometry,
                                                     std::span<const cfd::chemistry::AcidFamily> families,
                                                     double strong_charge_mol_per_litre = 0.0);
struct PittingCriterion {
    double critical_chloride{}, critical_potential{};
};
[[nodiscard]] bool pitting_initiates(double chloride_activity, double electrode_potential,
                                     const PittingCriterion& criterion);
[[nodiscard]] double stress_assisted_exchange_current(double base_exchange_current, double hydrostatic_stress,
                                                      double activation_volume, double temperature = 298.15);
[[nodiscard]] double phase_field_corrosion_step(double order_parameter, double laplacian,
                                                double electrochemical_driving, double mobility, double gradient_energy,
                                                double dt);
[[nodiscard]] double recessed_level_set(double signed_distance, double recession_distance);
[[nodiscard]] double marcus_hush_chidsey_current_density(double exchange_current_density, double overpotential,
                                                         double reorganization_energy_j_per_mol,
                                                         double temperature = 298.15, double electrons = 1.0,
                                                         std::size_t quadrature_points = 1200);
[[nodiscard]] std::vector<double> recession_remesh_1d(std::span<const double> coordinates, double recession_distance);
struct StressCorrosionUpdate {
    double stress_amplified_current_density{};
    double recession_distance{};
};
[[nodiscard]] StressCorrosionUpdate coupled_stress_corrosion_step(double base_exchange_current,
                                                                  double hydrostatic_stress, double activation_volume,
                                                                  double temperature, double dt, double molar_mass,
                                                                  double metal_density, double electrons = 1.0);
class PhaseFieldCorrosion2D {
public:
    PhaseFieldCorrosion2D(std::size_t nx, std::size_t ny, double dx, double dy, double mobility,
                          double gradient_energy);
    void initialize_interface(double x0, double y0, double width);
    void step(double dt, std::span<const double> driving);
    [[nodiscard]] const std::vector<double>& order_parameter() const noexcept { return phi_; }
    [[nodiscard]] std::size_t nx() const noexcept { return nx_; }
    [[nodiscard]] std::size_t ny() const noexcept { return ny_; }

private:
    std::size_t nx_{}, ny_{};
    double dx_{}, dy_{}, mobility_{}, gradient_{};
    std::vector<double> phi_, next_;
};
class PhaseFieldCorrosion1D {
public:
    PhaseFieldCorrosion1D(std::size_t cells, double dx, double mobility, double gradient_energy);
    void initialize_interface(double location, double width);
    void step(double dt, std::span<const double> driving);
    [[nodiscard]] const std::vector<double>& order_parameter() const noexcept { return phi_; }

private:
    std::size_t cells_{};
    double dx_{}, mobility_{}, gradient_{};
    std::vector<double> phi_, next_;
};

} // namespace cfd::electrochemistry
