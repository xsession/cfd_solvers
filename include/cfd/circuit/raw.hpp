#pragma once

#include "cfd/circuit/spice.hpp"

#include <complex>
#include <iosfwd>
#include <span>
#include <string>
#include <vector>

namespace cfd::circuit {

struct RawTrace {
    std::string name;
    std::string type{"voltage"};
    std::vector<Complex> values;
};

struct SpiceRawDataset {
    std::string title{"cfd_solvers"};
    std::string plot_name;
    std::string scale_name;
    std::string scale_type;
    bool complex_values{};
    std::vector<double> scale;
    std::vector<RawTrace> traces;
};

// ASCII SPICE RAW interchange baseline. The writer intentionally uses the
// documented text representation so files remain inspectable and portable.
void write_spice_raw(std::ostream& output,const SpiceRawDataset& dataset);
[[nodiscard]] SpiceRawDataset read_spice_raw(std::istream& input);

[[nodiscard]] SpiceRawDataset raw_from_transient(std::span<const TransientPoint> points,
                                                 std::span<const std::string> node_names,
                                                 std::string title = "cfd_solvers transient");
[[nodiscard]] SpiceRawDataset raw_from_ac(std::span<const AcSolution> points,
                                          std::span<const std::string> node_names,
                                          std::string title = "cfd_solvers AC analysis");

} // namespace cfd::circuit
