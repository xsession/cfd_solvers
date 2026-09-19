#include "cfd/circuit/osdi_adapter.hpp"

#include <stdexcept>
#include <utility>

namespace cfd::circuit {

OsdiDescriptorAdapter::OsdiDescriptorAdapter(std::size_t terminal_count, StampEvaluator evaluator)
    : terminal_count_(terminal_count), evaluator_(std::move(evaluator)) {
    if (terminal_count_ == 0U || !evaluator_)
        throw std::invalid_argument("invalid OSDI descriptor adapter");
}

void OsdiDescriptorAdapter::stamp(std::span<const double> terminal_voltage, std::span<double> terminal_current,
                                  std::span<double> jacobian) const {
    if (terminal_voltage.size() != terminal_count_ || terminal_current.size() != terminal_count_ ||
        jacobian.size() != terminal_count_ * terminal_count_) {
        throw std::invalid_argument("OSDI descriptor stamp buffer size mismatch");
    }
    evaluator_(terminal_voltage, terminal_current, jacobian);
}

OpenVafCompilePlan make_openvaf_compile_plan(std::string source, std::string output_library, std::string compiler) {
    if (source.empty() || output_library.empty() || compiler.empty())
        throw std::invalid_argument("OpenVAF compile plan paths must not be empty");
    return {{std::move(compiler), "--output", output_library, source}, std::move(source), std::move(output_library)};
}

} // namespace cfd::circuit
