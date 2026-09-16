#include "cfd/solvers/fem/poisson1d.hpp"

#include <stdexcept>

namespace cfd::fem {

Poisson1D::Poisson1D(std::size_t elements, double length)
    : elements_(elements), length_(length), x_(elements + 1), u_(elements + 1, 0.0) {
    if (elements < 2 || length <= 0.0) throw std::invalid_argument("Poisson1D needs >=2 elements and positive length");
    const double h = length_ / static_cast<double>(elements_);
    for (std::size_t i = 0; i <= elements_; ++i) x_[i] = static_cast<double>(i) * h;
}

void Poisson1D::solve(const std::function<double(double)>& source,
                      double left_value,
                      double right_value) {
    const std::size_t n = elements_ - 1;
    const double h = length_ / static_cast<double>(elements_);

    cfd::core::AlignedVector<double> lower(n, -1.0 / h);
    cfd::core::AlignedVector<double> diag(n, 2.0 / h);
    cfd::core::AlignedVector<double> upper(n, -1.0 / h);
    cfd::core::AlignedVector<double> rhs(n, 0.0);

    for (std::size_t j = 0; j < n; ++j) {
        const std::size_t node = j + 1;
        // Two adjacent linear elements; midpoint quadrature is exact for constant source.
        const double xl = x_[node] - 0.5 * h;
        const double xr = x_[node] + 0.5 * h;
        rhs[j] = 0.5 * h * (source(xl) + source(xr));
    }

    rhs.front() += left_value / h;
    rhs.back() += right_value / h;

    for (std::size_t i = 1; i < n; ++i) {
        const double m = lower[i] / diag[i - 1];
        diag[i] -= m * upper[i - 1];
        rhs[i] -= m * rhs[i - 1];
    }

    cfd::core::AlignedVector<double> interior(n, 0.0);
    interior[n - 1] = rhs[n - 1] / diag[n - 1];
    for (std::size_t k = n - 1; k-- > 0;) {
        interior[k] = (rhs[k] - upper[k] * interior[k + 1]) / diag[k];
    }

    u_.front() = left_value;
    u_.back() = right_value;
    for (std::size_t j = 0; j < n; ++j) u_[j + 1] = interior[j];
}

} // namespace cfd::fem
