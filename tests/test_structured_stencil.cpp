#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/structured_stencil.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

std::size_t index(std::size_t nx, std::size_t ny, std::size_t i, std::size_t j, std::size_t k) {
    return (k * ny + j) * nx + i;
}

double max_difference(std::span<const double> a, std::span<const double> b) {
    double error = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        error = std::max(error, std::abs(a[i] - b[i]));
    return error;
}

} // namespace

int main() {
    constexpr std::size_t nx = 9U;
    constexpr std::size_t ny = 7U;
    constexpr std::size_t nz = 5U;
    const std::size_t cells = nx * ny * nz;

    cfd::core::CsrBuilder builder(cells, cells);
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const auto row = index(nx, ny, i, j, k);
                builder.add(row, row, 6.0);
                if (i > 0U)
                    builder.add(row, index(nx, ny, i - 1U, j, k), -1.0);
                if (i + 1U < nx)
                    builder.add(row, index(nx, ny, i + 1U, j, k), -1.0);
                if (j > 0U)
                    builder.add(row, index(nx, ny, i, j - 1U, k), -1.0);
                if (j + 1U < ny)
                    builder.add(row, index(nx, ny, i, j + 1U, k), -1.0);
                if (k > 0U)
                    builder.add(row, index(nx, ny, i, j, k - 1U), -1.0);
                if (k + 1U < nz)
                    builder.add(row, index(nx, ny, i, j, k + 1U), -1.0);
            }
        }
    }

    const auto csr = builder.build();
    const cfd::core::StructuredSevenPointOperator stencil(nx, ny, nz);
    std::vector<double> x(cells), rhs(cells), csr_y(cells), stencil_y(cells), residual(cells);
    for (std::size_t i = 0; i < cells; ++i)
        x[i] = 0.25 + 0.001 * static_cast<double>((i * 17U) % 101U);

    csr.multiply(x, csr_y);
    stencil.apply(x, stencil_y);
    require(max_difference(csr_y, stencil_y) < 1.0e-12, "structured stencil must match CSR SpMV");

    for (std::size_t i = 0; i < cells; ++i)
        rhs[i] = csr_y[i] + 0.1 * std::sin(static_cast<double>(i));
    const double fused_norm = stencil.apply_residual_l2(rhs, x, residual);
    double explicit_norm = 0.0;
    for (std::size_t i = 0; i < cells; ++i) {
        const double expected = rhs[i] - csr_y[i];
        explicit_norm += expected * expected;
        require(std::abs(residual[i] - expected) < 1.0e-12, "fused residual value mismatch");
    }
    require(std::abs(fused_norm - std::sqrt(explicit_norm)) < 1.0e-12, "fused residual norm mismatch");
    require(stencil.estimated_apply_bytes() == cells * 8U * sizeof(double), "structured traffic estimate");
    std::cout << "structured stencil regression passed\n";
    return 0;
}
