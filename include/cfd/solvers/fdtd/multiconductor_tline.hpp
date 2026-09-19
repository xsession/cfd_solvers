#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::fdtd {

struct MulticonductorTlineConfig {
    std::size_t conductors{};
    std::vector<double> inductance_h;  // diagonal L baseline
    std::vector<double> capacitance_f; // diagonal C baseline
    double dt_s{};
};

struct MulticonductorTlineState {
    std::vector<double> voltage;
    std::vector<double> current;
};

void advance_multiconductor_tline(const MulticonductorTlineConfig& config, MulticonductorTlineState& state,
                                  std::span<const double> source_voltage);

} // namespace cfd::fdtd
