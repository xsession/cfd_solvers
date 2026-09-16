#pragma once

#include "cfd/core/aligned_allocator.hpp"

#include <cstddef>
#include <functional>

namespace cfd::fem {

class Poisson1D {
public:
    explicit Poisson1D(std::size_t elements, double length = 1.0);

    void solve(const std::function<double(double)>& source,
               double left_value = 0.0,
               double right_value = 0.0);

    [[nodiscard]] const cfd::core::AlignedVector<double>& x() const noexcept { return x_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& u() const noexcept { return u_; }

private:
    std::size_t elements_{};
    double length_{};
    cfd::core::AlignedVector<double> x_;
    cfd::core::AlignedVector<double> u_;
};

} // namespace cfd::fem
