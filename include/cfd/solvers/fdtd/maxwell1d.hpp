#pragma once

#include "cfd/core/aligned_allocator.hpp"

#include <cstddef>

namespace cfd::fdtd {

struct Maxwell1DConfig {
    std::size_t cells{1024};
    double dx{1.0e-3};
    double courant{0.99};
    double epsilon_r{1.0};
    double mu_r{1.0};
};

class Maxwell1D {
public:
    explicit Maxwell1D(Maxwell1DConfig config);

    void initialize_gaussian(double center_fraction = 0.3, double width_fraction = 0.04);
    void step(std::size_t count = 1);

    [[nodiscard]] double energy() const;
    [[nodiscard]] const cfd::core::AlignedVector<double>& electric() const noexcept { return ez_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& magnetic() const noexcept { return hy_; }
    [[nodiscard]] double dt() const noexcept { return dt_; }

private:
    Maxwell1DConfig config_;
    double dt_{};
    double epsilon_{};
    double mu_{};
    cfd::core::AlignedVector<double> ez_;
    cfd::core::AlignedVector<double> hy_;
};

} // namespace cfd::fdtd
