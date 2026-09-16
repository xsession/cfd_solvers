#include "cfd/solvers/fvm/diffusion2d.hpp"

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fvm {

Diffusion2D::Diffusion2D(Diffusion2DConfig config)
    : config_(config),
      grid_(config.nx, config.ny,
            config.lx / static_cast<double>(config.nx - 1),
            config.ly / static_cast<double>(config.ny - 1)),
      phi_(grid_.cells(), 0.0),
      next_(grid_.cells(), 0.0),
      source_(grid_.cells(), 0.0) {
    if (config_.diffusivity <= 0.0) throw std::invalid_argument("diffusivity must be positive");
}

void Diffusion2D::set_source(double value) {
    std::fill(source_.begin(), source_.end(), value);
}

void Diffusion2D::set_dirichlet(double left, double right, double bottom, double top) {
    left_ = left;
    right_ = right;
    bottom_ = bottom;
    top_ = top;
    apply_boundaries(phi_);
    apply_boundaries(next_);
}

void Diffusion2D::apply_boundaries(cfd::core::AlignedVector<double>& field) {
    for (std::size_t y = 0; y < grid_.ny; ++y) {
        field[grid_.index(0, y)] = left_;
        field[grid_.index(grid_.nx - 1, y)] = right_;
    }
    for (std::size_t x = 0; x < grid_.nx; ++x) {
        field[grid_.index(x, 0)] = bottom_;
        field[grid_.index(x, grid_.ny - 1)] = top_;
    }
}

double Diffusion2D::iterate(std::size_t iterations) {
    double max_delta = 0.0;
    const double idx2 = 1.0 / (grid_.dx * grid_.dx);
    const double idy2 = 1.0 / (grid_.dy * grid_.dy);
    const double denom = 2.0 * config_.diffusivity * (idx2 + idy2);

    for (std::size_t iter = 0; iter < iterations; ++iter) {
        cfd::core::parallel_for(grid_.cells(), [&](std::size_t i) {
            const std::size_t x = i % grid_.nx;
            const std::size_t y = i / grid_.nx;
            if (x == 0 || y == 0 || x + 1 == grid_.nx || y + 1 == grid_.ny) {
                next_[i] = phi_[i];
                return;
            }
            const double east = phi_[grid_.index(x + 1, y)];
            const double west = phi_[grid_.index(x - 1, y)];
            const double north = phi_[grid_.index(x, y + 1)];
            const double south = phi_[grid_.index(x, y - 1)];
            next_[i] = (config_.diffusivity * (idx2 * (east + west) + idy2 * (north + south))
                       + source_[i]) / denom;
        });

        apply_boundaries(next_);
        max_delta = 0.0;
        for (std::size_t i = 0; i < grid_.cells(); ++i) {
            max_delta = std::max(max_delta, std::abs(next_[i] - phi_[i]));
        }
        phi_.swap(next_);
    }
    return max_delta;
}

} // namespace cfd::fvm
