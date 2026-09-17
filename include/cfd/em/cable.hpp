#pragma once

#include "cfd/rf/nport.hpp"

#include <complex>
#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::em {

using Complex = std::complex<double>;

struct MulticonductorRlcg {
    std::size_t conductors{};
    // Row-major per-unit-length matrices [ohm/m, H/m, S/m, F/m].
    std::vector<double> resistance;
    std::vector<double> inductance;
    std::vector<double> conductance;
    std::vector<double> capacitance;
    void validate() const;
};

struct CableFrequencyResult {
    double frequency_hz{};
    std::vector<Complex> input_voltage_v;
    std::vector<Complex> input_current_a;
    std::vector<Complex> load_voltage_v;
    std::vector<Complex> load_current_a;
    cfd::rf::ComplexMatrix transfer_matrix;
};

// Frequency-domain multiconductor telegrapher solve. A [1/1] Pade step of the
// first-order RLCG system is cascaded along the cable, then the unknown input
// currents are solved from the full matrix load termination V_L=Z_L I_L.
[[nodiscard]] CableFrequencyResult solve_multiconductor_cable(
    double length_m,std::size_t segments,double frequency_hz,
    const MulticonductorRlcg& rlcg,
    std::span<const Complex> imposed_input_voltage_v,
    const cfd::rf::ComplexMatrix& load_impedance_ohm);

using RlcgSampler = std::function<MulticonductorRlcg(double frequency_hz)>;
[[nodiscard]] std::vector<CableFrequencyResult> sweep_multiconductor_cable(
    double length_m,std::size_t segments,std::span<const double> frequencies_hz,
    const RlcgSampler& rlcg_sampler,std::span<const Complex> imposed_input_voltage_v,
    const cfd::rf::ComplexMatrix& load_impedance_ohm);

struct CableShieldTransferModel {
    double dc_transfer_resistance_ohm_per_m{};
    double transfer_inductance_h_per_m{};
    double skin_corner_hz{1.0e6};
};
[[nodiscard]] Complex cable_shield_transfer_impedance_ohm_per_m(
    double frequency_hz,const CableShieldTransferModel& model);
[[nodiscard]] Complex cable_shield_induced_voltage_v(
    double frequency_hz,Complex shield_current_a,double exposed_length_m,
    const CableShieldTransferModel& model);
[[nodiscard]] Complex field_to_cable_open_circuit_voltage_v(
    Complex parallel_electric_field_v_per_m,double effective_height_m,
    Complex coupling_factor={1.0,0.0});

} // namespace cfd::em
