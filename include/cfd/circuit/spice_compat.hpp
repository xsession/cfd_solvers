#pragma once

#include <string_view>
#include <vector>

namespace cfd::circuit {

enum class SpiceCompatibilityMode { native, pspice, ltspice, hspice };

struct SpiceCompatibilityProfile {
    SpiceCompatibilityMode mode{};
    bool engineering_suffixes{};
    bool common_model_cards{};
    bool subcircuits_and_parameters{};
    std::vector<std::string_view> supported_directives;
};

// Describes the clean parser's deliberately portable interchange subset. The
// profile is explicit so applications can reject vendor extensions instead of
// silently treating them as native syntax.
[[nodiscard]] inline SpiceCompatibilityProfile spice_compatibility_profile(SpiceCompatibilityMode mode) {
    return {
        mode,
        true,
        true,
        true,
        {".MODEL", ".SUBCKT", ".PARAM", ".FUNC", ".INCLUDE", ".LIB", ".IC", ".NODESET", ".TRAN", ".AC", ".DC", ".OP"}};
}

} // namespace cfd::circuit
