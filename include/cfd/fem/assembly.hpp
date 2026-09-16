#pragma once

#include "cfd/core/csr_matrix.hpp"
#include "cfd/fem/mesh2d.hpp"

#include <functional>
#include <span>

namespace cfd::fem {

[[nodiscard]] cfd::core::CsrMatrix assemble_tri3_laplace_matrix(
    const Mesh2D& mesh,
    const std::function<double(Node2)>& diffusivity);

// Matrix-free full-node Laplace action. No essential boundary conditions are
// applied here; callers can constrain/project the result as appropriate.
void apply_tri3_laplace_matrix_free(const Mesh2D& mesh,
                                    const std::function<double(Node2)>& diffusivity,
                                    std::span<const double> x,
                                    std::span<double> y);

} // namespace cfd::fem
