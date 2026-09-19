#pragma once

namespace cfd::fvm {

struct SootRadiationConfig {
    double soot_yield{0.02};
    double oxidation_rate{1.0};
    double absorption_per_mass{2.0e4};
    double emissivity_floor{0.0};
};

struct SootRadiationState {
    double soot_mass_fraction{};
    double temperature{300.0};
    double heat_release{};
};

struct SootRadiationResult {
    double soot_source{};
    double absorption_coefficient{};
    double radiative_source{};
};

// Compact soot source/radiation hook. It is intentionally a model object,
// not a hidden chemistry side effect: callers can feed soot_source into a
// transported scalar and radiative_source into ThermalTransport.
class SootRadiationModel {
public:
    explicit SootRadiationModel(SootRadiationConfig config = {});
    [[nodiscard]] SootRadiationResult evaluate(SootRadiationState state, double environment_temperature) const;

private:
    SootRadiationConfig config_;
};

} // namespace cfd::fvm
