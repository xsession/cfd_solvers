#pragma once

#include "cfd/core/parallel.hpp"

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>

namespace cfd::core {

// A regular-grid seven-point operator for cell-centered scalar problems.
// The operator intentionally has no mesh/index indirection: it is the
// structured counterpart to CsrMatrix for regular Cartesian stencils.
class StructuredSevenPointOperator {
public:
    StructuredSevenPointOperator(std::size_t nx, std::size_t ny, std::size_t nz, double diagonal = 6.0,
                                 double neighbor = -1.0)
        : nx_(nx), ny_(ny), nz_(nz), diagonal_(diagonal), neighbor_(neighbor) {
        if (nx_ == 0U || ny_ == 0U || nz_ == 0U || !std::isfinite(diagonal_) || !std::isfinite(neighbor_)) {
            throw std::invalid_argument("invalid structured seven-point operator dimensions or coefficients");
        }
    }

    [[nodiscard]] std::size_t nx() const noexcept { return nx_; }
    [[nodiscard]] std::size_t ny() const noexcept { return ny_; }
    [[nodiscard]] std::size_t nz() const noexcept { return nz_; }
    [[nodiscard]] std::size_t cells() const noexcept { return nx_ * ny_ * nz_; }
    [[nodiscard]] double diagonal() const noexcept { return diagonal_; }
    [[nodiscard]] double neighbor() const noexcept { return neighbor_; }

    // These are algorithmic traffic estimates, not hardware counters. They
    // make results comparable across implementations and architectures.
    [[nodiscard]] std::size_t estimated_apply_bytes() const noexcept { return cells() * 8U * sizeof(double); }

    [[nodiscard]] std::size_t estimated_fused_residual_bytes() const noexcept { return cells() * 9U * sizeof(double); }

    void apply(std::span<const double> x, std::span<double> y) const {
        validate_vectors(x, y);
        parallel_for_weighted(ny_ * nz_, nx_, [&](std::size_t line) { apply_line(x, y, line); });
    }

    // Computes rhs - A*x and its L2 norm in one field traversal. This is the
    // structured equivalent of fusing operator application with convergence
    // monitoring, avoiding a second full-grid residual pass.
    [[nodiscard]] double apply_residual_l2(std::span<const double> rhs, std::span<const double> x,
                                           std::span<double> residual) const {
        if (rhs.size() != cells() || x.size() != cells() || residual.size() != cells()) {
            throw std::invalid_argument("structured stencil residual buffer mismatch");
        }
        const double squared = parallel_sum_weighted(
            ny_ * nz_, nx_, [&](std::size_t line) { return residual_line(rhs, x, residual, line); });
        return std::sqrt(squared);
    }

private:
    std::size_t nx_{};
    std::size_t ny_{};
    std::size_t nz_{};
    double diagonal_{};
    double neighbor_{};

    void apply_line(std::span<const double> x, std::span<double> y, std::size_t line) const noexcept {
        const std::size_t row_start = line * nx_;
        const std::size_t row_end = row_start + nx_;
        const std::size_t plane = nx_ * ny_;
        const bool first_row = line % ny_ == 0U;
        const bool last_row = (line + 1U) % ny_ == 0U;
        const bool first_plane = line < ny_;
        const bool last_plane = line + ny_ >= ny_ * nz_;
        for (std::size_t flat = row_start; flat < row_end; ++flat) {
            double value = diagonal_ * x[flat];
            if (flat > row_start)
                value += neighbor_ * x[flat - 1U];
            if (flat + 1U < row_end)
                value += neighbor_ * x[flat + 1U];
            if (!first_row)
                value += neighbor_ * x[flat - nx_];
            if (!last_row)
                value += neighbor_ * x[flat + nx_];
            if (!first_plane)
                value += neighbor_ * x[flat - plane];
            if (!last_plane)
                value += neighbor_ * x[flat + plane];
            y[flat] = value;
        }
    }

    [[nodiscard]] double residual_line(std::span<const double> rhs, std::span<const double> x,
                                       std::span<double> residual, std::size_t line) const noexcept {
        const std::size_t row_start = line * nx_;
        const std::size_t row_end = row_start + nx_;
        const std::size_t plane = nx_ * ny_;
        const bool first_row = line % ny_ == 0U;
        const bool last_row = (line + 1U) % ny_ == 0U;
        const bool first_plane = line < ny_;
        const bool last_plane = line + ny_ >= ny_ * nz_;
        double squared = 0.0;
        for (std::size_t flat = row_start; flat < row_end; ++flat) {
            double value = diagonal_ * x[flat];
            if (flat > row_start)
                value += neighbor_ * x[flat - 1U];
            if (flat + 1U < row_end)
                value += neighbor_ * x[flat + 1U];
            if (!first_row)
                value += neighbor_ * x[flat - nx_];
            if (!last_row)
                value += neighbor_ * x[flat + nx_];
            if (!first_plane)
                value += neighbor_ * x[flat - plane];
            if (!last_plane)
                value += neighbor_ * x[flat + plane];
            const double value_residual = rhs[flat] - value;
            residual[flat] = value_residual;
            squared += value_residual * value_residual;
        }
        return squared;
    }

    void validate_vectors(std::span<const double> x, std::span<double> y) const {
        if (x.size() != cells() || y.size() != cells()) {
            throw std::invalid_argument("structured stencil buffer mismatch");
        }
    }
};

} // namespace cfd::core
