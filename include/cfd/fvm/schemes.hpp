#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <span>
#include <vector>

namespace cfd::fvm {

enum class FaceInterpolationScheme {
    linear,
    upwind,
    // MUSCL-style upwind reconstruction with a Barth-Jespersen limiter.
    // Second-order on smooth linear fields and bounded by the local stencil.
    bounded_linear,
    minmod,
    van_leer
};

struct ScalarSchemeWorkspace {
    std::vector<Vec3> gradient;
    std::vector<double> minimum,maximum,limiter,face_values;
    void resize(const PolyMesh& mesh);
    [[nodiscard]] std::size_t allocated_bytes() const noexcept;
};

// Output and input spans must not overlap. Resizing the workspace once makes
// subsequent scalar reconstruction/divergence calls allocation-free.
void interpolate_scalar_to_faces_into(const PolyMesh& mesh,std::span<const double> cell_values,
    std::span<double> output,ScalarSchemeWorkspace& workspace,
    FaceInterpolationScheme scheme = FaceInterpolationScheme::linear,
    std::span<const double> face_flux = {},std::span<const double> boundary_face_values = {});
void convective_divergence_scalar_into(const PolyMesh& mesh,std::span<const double> cell_values,
    std::span<const double> face_flux,std::span<double> output,ScalarSchemeWorkspace& workspace,
    FaceInterpolationScheme scheme = FaceInterpolationScheme::upwind,
    std::span<const double> boundary_face_values = {});

// face_flux is oriented with Face::area: positive means owner -> neighbour/outside.
// Boundary values are indexed by global face index. Empty boundary spans request
// zero-gradient/outflow behavior.
[[nodiscard]] std::vector<double> interpolate_scalar_to_faces(
    const PolyMesh& mesh,
    std::span<const double> cell_values,
    FaceInterpolationScheme scheme = FaceInterpolationScheme::linear,
    std::span<const double> face_flux = {},
    std::span<const double> boundary_face_values = {});

[[nodiscard]] std::vector<Vec3> interpolate_vector_to_faces(
    const PolyMesh& mesh,
    std::span<const Vec3> cell_values,
    FaceInterpolationScheme scheme = FaceInterpolationScheme::linear,
    std::span<const double> face_flux = {},
    std::span<const Vec3> boundary_face_values = {});

[[nodiscard]] std::vector<double> face_flux_from_velocity(
    const PolyMesh& mesh,
    std::span<const Vec3> cell_velocity,
    std::span<const Vec3> boundary_face_velocity = {});

[[nodiscard]] std::vector<double> convective_divergence_scalar(
    const PolyMesh& mesh,
    std::span<const double> cell_values,
    std::span<const double> face_flux,
    FaceInterpolationScheme scheme = FaceInterpolationScheme::upwind,
    std::span<const double> boundary_face_values = {});

[[nodiscard]] std::vector<Vec3> convective_divergence_vector(
    const PolyMesh& mesh,
    std::span<const Vec3> cell_values,
    std::span<const double> face_flux,
    FaceInterpolationScheme scheme = FaceInterpolationScheme::upwind,
    std::span<const Vec3> boundary_face_values = {});

} // namespace cfd::fvm
