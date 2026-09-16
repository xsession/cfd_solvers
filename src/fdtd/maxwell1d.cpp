#include "cfd/solvers/fdtd/maxwell1d.hpp"

#include "cfd/core/parallel.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::fdtd {

namespace {
constexpr double epsilon0 = 8.8541878128e-12;
constexpr double mu0 = 1.25663706212e-6;
}

Maxwell1D::Maxwell1D(Maxwell1DConfig config)
    : config_(config),
      epsilon_(epsilon0 * config.epsilon_r),
      mu_(mu0 * config.mu_r),
      ez_(config.cells, 0.0),
      hy_(config.cells - 1, 0.0) {
    if (config.cells < 4 || config.dx <= 0.0 || config.courant <= 0.0 || config.courant >= 1.0 ||
        config.epsilon_r <= 0.0 || config.mu_r <= 0.0) {
        throw std::invalid_argument("invalid Maxwell1D configuration");
    }
    const double wave_speed = 1.0 / std::sqrt(epsilon_ * mu_);
    dt_ = config.courant * config.dx / wave_speed;
}

void Maxwell1D::initialize_gaussian(double center_fraction, double width_fraction) {
    const double center = center_fraction * static_cast<double>(config_.cells - 1);
    const double width = width_fraction * static_cast<double>(config_.cells);
    for (std::size_t i = 0; i < config_.cells; ++i) {
        const double z = (static_cast<double>(i) - center) / width;
        ez_[i] = std::exp(-0.5 * z * z);
    }
}

void Maxwell1D::step(std::size_t count) {
    const double h_coeff = dt_ / (mu_ * config_.dx);
    const double e_coeff = dt_ / (epsilon_ * config_.dx);
    for (std::size_t step_i = 0; step_i < count; ++step_i) {
        cfd::core::parallel_for(hy_.size(), [&](std::size_t i) {
            hy_[i] += h_coeff * (ez_[i] - ez_[i + 1]);
        });
        cfd::core::parallel_for(config_.cells - 2, [&](std::size_t j) {
            const std::size_t i = j + 1;
            ez_[i] += e_coeff * (hy_[i - 1] - hy_[i]);
        });
        // Simple PEC boundaries. Higher-order ABC/PML belongs to the EM roadmap.
        ez_.front() = 0.0;
        ez_.back() = 0.0;
    }
}

double Maxwell1D::energy() const {
    double e = 0.0;
    for (double v : ez_) e += 0.5 * epsilon_ * v * v * config_.dx;
    for (double v : hy_) e += 0.5 * mu_ * v * v * config_.dx;
    return e;
}

} // namespace cfd::fdtd
