#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/grid.hpp"

#include <cstddef>

namespace cfd::fvm {

struct Diffusion2DConfig {
    std::size_t nx{128};
    std::size_t ny{128};
    double lx{1.0};
    double ly{1.0};
    double diffusivity{1.0};
};

class Diffusion2D {
public:
    explicit Diffusion2D(Diffusion2DConfig config);

    void set_source(double value);
    void set_dirichlet(double left, double right, double bottom, double top);
    [[nodiscard]] double iterate(std::size_t iterations = 1);

    [[nodiscard]] const cfd::core::AlignedVector<double>& field() const noexcept { return phi_; }
    [[nodiscard]] const cfd::core::Grid2D& grid() const noexcept { return grid_; }

private:
    Diffusion2DConfig config_;
    cfd::core::Grid2D grid_;
    cfd::core::AlignedVector<double> phi_;
    cfd::core::AlignedVector<double> next_;
    cfd::core::AlignedVector<double> source_;
    double left_{0.0};
    double right_{0.0};
    double bottom_{0.0};
    double top_{0.0};

    void apply_boundaries(cfd::core::AlignedVector<double>& field);
};

} // namespace cfd::fvm
