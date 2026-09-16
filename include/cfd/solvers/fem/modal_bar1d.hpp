#pragma once

#include "cfd/core/eigen_solvers.hpp"

#include <cstddef>
#include <vector>

namespace cfd::fem {

struct BarModal1DConfig {
    std::size_t elements{80U};
    double length{1.0};
    double area{1.0};
    double young_modulus{200.0e9};
    double density{7800.0};
    cfd::core::GeneralizedEigenConfig eigen{};
};

struct BarMode1D {
    double angular_frequency{};
    double frequency_hz{};
    std::vector<double> displacement; // includes fixed node 0
    cfd::core::GeneralizedEigenpair eigenpair;
};

// Axial vibration of a uniform fixed-free bar using Line2 FEM with consistent
// mass. Provides a small, deterministic generalized-eigenvalue regression for
// the structural dynamics/eigenfrequency stack.
class BarModal1D {
public:
    explicit BarModal1D(BarModal1DConfig config = {});
    [[nodiscard]] std::vector<BarMode1D> solve(std::size_t mode_count) const;
    [[nodiscard]] const BarModal1DConfig& config() const noexcept { return config_; }

private:
    BarModal1DConfig config_;
};

} // namespace cfd::fem
