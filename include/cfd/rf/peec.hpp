#pragma once

#include <complex>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace cfd::rf {

struct Point3 {
    double x{}, y{}, z{};
};

struct FilamentSegment {
    Point3 start{}, end{};
    double radius_m{};
    double conductivity_s_per_m{5.8e7};
};

struct PeecMatrices {
    std::size_t size{};
    std::vector<double> resistance_ohm;      // diagonal branch resistance
    std::vector<double> partial_inductance_h; // row-major Lp matrix
};

struct PeecCapacitanceMatrices {
    std::size_t size{};
    // Row-major coefficient-of-potential matrix P where V = P q.
    std::vector<double> coefficient_of_potential_v_per_c;
    // Row-major Maxwell capacitance matrix C = P^-1.
    std::vector<double> capacitance_f;
};

class PeecFilamentSystem {
public:
    void add_segment(FilamentSegment segment);
    [[nodiscard]] const std::vector<FilamentSegment>& segments() const noexcept { return segments_; }
    [[nodiscard]] PeecMatrices extract(std::size_t quadrature_order = 6U) const;
    [[nodiscard]] PeecCapacitanceMatrices extract_capacitance(double relative_permittivity = 1.0,
                                                               std::size_t quadrature_order = 6U) const;
    [[nodiscard]] std::vector<std::complex<double>> impedance_matrix(double frequency_hz,
                                                                     std::size_t quadrature_order = 6U) const;
    [[nodiscard]] std::vector<std::complex<double>> impedance_matrix_skin_effect(double frequency_hz,
                                                                                 std::size_t quadrature_order = 6U) const;
    // Emits a coupled-inductor SPICE subcircuit with one terminal pair per filament.
    // Mutual coupling is represented by K elements derived from the PEEC L matrix.
    [[nodiscard]] std::string to_spice_subcircuit(std::string_view name = "PEEC",
                                                  std::size_t quadrature_order = 6U) const;
    // Solve Z I = V for branch currents. Intended as a quasi-static conductor
    // extraction baseline, not a radiation-capable full-wave antenna solver.
    [[nodiscard]] std::vector<std::complex<double>> solve_currents(
        double frequency_hz,std::span<const std::complex<double>> applied_branch_voltage,
        std::size_t quadrature_order = 6U) const;

private:
    std::vector<FilamentSegment> segments_;
};

[[nodiscard]] double filament_length(const FilamentSegment& segment);
[[nodiscard]] double dc_resistance(const FilamentSegment& segment);
[[nodiscard]] double self_partial_inductance(const FilamentSegment& segment);
[[nodiscard]] double skin_depth_m(double conductivity_s_per_m,double frequency_hz,double relative_permeability = 1.0);
[[nodiscard]] double round_wire_ac_resistance(const FilamentSegment& segment,double frequency_hz,
                                              double relative_permeability = 1.0);
[[nodiscard]] double mutual_partial_inductance(const FilamentSegment& first,
                                                const FilamentSegment& second,
                                                std::size_t quadrature_order = 6U);

} // namespace cfd::rf
