#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::fdtd {

[[nodiscard]] std::vector<double> make_graded_axis(double length_m, std::size_t cells, double first_cell_width_m,
                                                   double growth_ratio);

struct SubcellWeight {
    std::size_t i{}, j{}, k{};
    double weight{};
};
[[nodiscard]] std::vector<SubcellWeight> deposit_thin_wire_subcell(double x_m, double y_m, double z_m, std::size_t nx,
                                                                   std::size_t ny, std::size_t nz, double dx_m,
                                                                   double dy_m, double dz_m);

struct HuygensEquivalentSource {
    std::array<double, 3> electric_surface_current{}; // n x H
    std::array<double, 3> magnetic_surface_current{}; // -n x E
};
[[nodiscard]] HuygensEquivalentSource huygens_equivalent_source(std::array<double, 3> normal,
                                                                std::array<double, 3> electric_field,
                                                                std::array<double, 3> magnetic_field);

[[nodiscard]] std::complex<double> impedance_sheet_transmission(std::complex<double> sheet_impedance_ohm,
                                                                std::complex<double> medium1_impedance_ohm,
                                                                std::complex<double> medium2_impedance_ohm);

struct EmcTransferSample {
    double frequency_hz{};
    std::complex<double> gain{};
};
[[nodiscard]] std::complex<double> interpolate_emc_transfer(std::span<const EmcTransferSample> samples,
                                                            double frequency_hz);

} // namespace cfd::fdtd
