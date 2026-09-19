#pragma once

#include <cstddef>

namespace cfd::fvm {

struct FlameModelState {
    double progress{};         // 0 = unburned, 1 = burned
    double mixture_fraction{}; // 0 = oxidizer, 1 = fuel
    double temperature{300.0};
};

struct FlameModelResult {
    double progress_source{};
    double mixture_fraction_source{};
    double heat_release{};
};

class FlameModel {
public:
    virtual ~FlameModel() = default;
    [[nodiscard]] virtual FlameModelResult evaluate(FlameModelState state) const = 0;
};

struct PremixedFlameModelConfig {
    double reaction_rate{1.0};
    double heat_release{1.0e6};
    double ignition_temperature{300.0};
};

class PremixedFlameModel final : public FlameModel {
public:
    explicit PremixedFlameModel(PremixedFlameModelConfig config = {});
    [[nodiscard]] FlameModelResult evaluate(FlameModelState state) const override;

private:
    PremixedFlameModelConfig config_;
};

struct NonPremixedFlameModelConfig {
    double stoichiometric_mixture_fraction{0.5};
    double mixing_width{0.08};
    double mixing_rate{10.0};
    double heat_release{1.0e6};
    double ignition_temperature{300.0};
};

// A mixture-fraction/progress interface for diffusion flames. The Gaussian
// stoichiometric window is a compact, bounded hook that can later be replaced
// by a flamelet table without changing the FVM transport contract.
class NonPremixedFlameModel final : public FlameModel {
public:
    explicit NonPremixedFlameModel(NonPremixedFlameModelConfig config = {});
    [[nodiscard]] FlameModelResult evaluate(FlameModelState state) const override;

private:
    NonPremixedFlameModelConfig config_;
};

} // namespace cfd::fvm
