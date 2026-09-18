#pragma once

#include <cstddef>
#include <vector>

namespace cfd::fem {

struct PiezoelectricBar1DConfig {
    std::size_t elements{16};
    double length_m{1.0e-3};
    double area_m2{1.0e-4};
    double elastic_modulus_pa{6.0e10};
    double piezoelectric_stress_c_per_m2{15.0};
    double permittivity_f_per_m{1.5e-8};
};

// Small-strain 1-D coupled piezoelectric finite-element baseline. Node 0 is
// mechanically fixed and electrically grounded. The far electrode is open
// circuit unless a voltage is prescribed.
class PiezoelectricBar1D {
public:
    explicit PiezoelectricBar1D(PiezoelectricBar1DConfig config = {});
    void solve_receive(double tip_force_n);
    void solve_transmit(double applied_voltage_v, double tip_force_n = 0.0);

    [[nodiscard]] const std::vector<double>& displacement_m() const noexcept { return displacement_; }
    [[nodiscard]] const std::vector<double>& potential_v() const noexcept { return potential_; }
    [[nodiscard]] double open_circuit_voltage_v() const noexcept { return potential_.back(); }
    [[nodiscard]] double tip_displacement_m() const noexcept { return displacement_.back(); }
    [[nodiscard]] double electrode_charge_c() const noexcept { return electrode_charge_c_; }
    [[nodiscard]] const PiezoelectricBar1DConfig& config() const noexcept { return config_; }
private:
    PiezoelectricBar1DConfig config_;
    std::vector<double> displacement_, potential_;
    double electrode_charge_c_{};
    void solve(double tip_force_n, bool prescribe_tip_voltage, double tip_voltage_v);
};

} // namespace cfd::fem
