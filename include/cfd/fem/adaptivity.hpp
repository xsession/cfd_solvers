#pragma once

#include "cfd/fem/mesh2d.hpp"

#include <functional>
#include <span>
#include <vector>

namespace cfd::fem {

struct PoissonErrorEstimate2D {
    std::vector<double> element_indicator;
    double global_indicator{};
};

// Residual/jump estimator for -div(k grad u)=f with P1 triangles. The element
// residual is represented by h_K^2 ||f||^2 because the affine P1 Laplacian is
// zero inside each element; interior-edge normal-flux jumps provide the second
// contribution. This is an estimator, not an exact error norm.
[[nodiscard]] PoissonErrorEstimate2D estimate_poisson_error_tri3(
    const Mesh2D& mesh,
    std::span<const double> nodal_solution,
    const std::function<double(Node2)>& source,
    double diffusivity = 1.0);

// Mark the largest indicators until at least target_fraction of the squared
// global estimator is represented (Dorfler marking). target_fraction in (0,1].
[[nodiscard]] std::vector<unsigned char> mark_dorfler(
    std::span<const double> indicator,
    double target_fraction);

// Conforming longest-edge adaptive refinement. Every marked triangle selects
// its longest edge. All triangles incident to a selected edge are subdivided;
// triangles touched by 1/2/3 split edges are triangulated with conforming
// templates. Boundary patch ids are preserved on split boundary edges.
[[nodiscard]] Mesh2D refine_tri3_longest_edges(
    const Mesh2D& mesh,
    std::span<const unsigned char> marked_element);

} // namespace cfd::fem
