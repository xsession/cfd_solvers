#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <span>
#include <vector>

namespace cfd::fvm {

void gauss_gradient_scalar_into(const PolyMesh& mesh,std::span<const double> cell_values,
                                std::span<Vec3> result,std::span<const double> boundary_face_values = {});
void least_squares_gradient_scalar_into(const PolyMesh& mesh,std::span<const double> cell_values,
                                        std::span<Vec3> result,std::span<const double> boundary_face_values = {});

// Boundary values are indexed by global face index. An empty span requests zero-gradient boundaries.
[[nodiscard]] std::vector<Vec3> gauss_gradient_scalar(const PolyMesh& mesh,
                                                       std::span<const double> cell_values,
                                                       std::span<const double> boundary_face_values = {});

[[nodiscard]] std::vector<double> gauss_divergence_vector(const PolyMesh& mesh,
                                                           std::span<const Vec3> cell_values,
                                                           std::span<const Vec3> boundary_face_values = {});

[[nodiscard]] std::vector<double> orthogonal_laplacian_scalar(const PolyMesh& mesh,
                                                               std::span<const double> cell_values,
                                                               double diffusivity = 1.0,
                                                               std::span<const double> boundary_face_values = {});

[[nodiscard]] std::vector<Vec3> least_squares_gradient_scalar(
    const PolyMesh& mesh,
    std::span<const double> cell_values,
    std::span<const double> boundary_face_values = {});

// Face-normal derivative using owner/neighbour orthogonal difference plus a
// deferred non-orthogonal correction from least-squares cell gradients.
[[nodiscard]] std::vector<double> corrected_face_normal_gradient_scalar(
    const PolyMesh& mesh,
    std::span<const double> cell_values,
    std::span<const double> boundary_face_values = {});

[[nodiscard]] std::vector<double> corrected_laplacian_scalar(
    const PolyMesh& mesh,
    std::span<const double> cell_values,
    double diffusivity = 1.0,
    std::span<const double> boundary_face_values = {});

} // namespace cfd::fvm
