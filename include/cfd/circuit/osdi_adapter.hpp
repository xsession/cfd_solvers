#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

namespace cfd::circuit {

// Stable clean-core stamping seam around an external OSDI descriptor. An
// adapter owns the vendor ABI object and exposes only terminal currents and a
// dense Jacobian to Circuit's Newton/MNA layer.
class OsdiDescriptorAdapter {
public:
    using StampEvaluator = std::function<void(std::span<const double>, std::span<double>, std::span<double>)>;

    OsdiDescriptorAdapter(std::size_t terminal_count, StampEvaluator evaluator);
    void stamp(std::span<const double> terminal_voltage, std::span<double> terminal_current,
               std::span<double> jacobian_row_major) const;
    [[nodiscard]] std::size_t terminal_count() const noexcept { return terminal_count_; }

private:
    std::size_t terminal_count_{};
    StampEvaluator evaluator_;
};

struct OpenVafCompilePlan {
    std::vector<std::string> argv;
    std::string source;
    std::string output_library;
};

// Generates a reproducible OpenVAF command without making the optional tool a
// dependency of the solver build.
[[nodiscard]] OpenVafCompilePlan make_openvaf_compile_plan(std::string source, std::string output_library,
                                                           std::string compiler = "openvaf");

} // namespace cfd::circuit
