#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <span>
#include <vector>

namespace cfd::fvm {

enum class VelocityBoundaryType {
    fixedValue,
    zeroGradient,
    slip
};

struct VelocityBoundaryCondition {
    VelocityBoundaryType type{VelocityBoundaryType::zeroGradient};
    Vec3 value{};
};

enum class PressureBoundaryType {
    fixedValue,
    zeroGradient
};

struct PressureBoundaryCondition {
    PressureBoundaryType type{PressureBoundaryType::zeroGradient};
    double value{};
};

struct FaceOrthogonalDecomposition {
    Vec3 delta{};                 // owner -> neighbour/boundary face
    Vec3 orthogonal_area{};       // parallel to delta
    Vec3 nonorthogonal_area{};    // remainder
    double orthogonal_metric{};   // dot(Sf,delta)/|delta|^2
};

[[nodiscard]] FaceOrthogonalDecomposition decompose_face_area(
    const PolyMesh& mesh,
    std::size_t face_index);

[[nodiscard]] std::vector<Vec3> boundary_velocity_values(
    const PolyMesh& mesh,
    std::span<const Vec3> cell_velocity,
    std::span<const VelocityBoundaryCondition> patch_conditions);

[[nodiscard]] std::vector<double> boundary_pressure_values(
    const PolyMesh& mesh,
    std::span<const double> cell_pressure,
    std::span<const PressureBoundaryCondition> patch_conditions);

// Rhie-Chow-style collocated face flux. h_by_a is the pressure-free momentum
// predictor and pressure_mobility is V/aP for kinematic pressure.
[[nodiscard]] std::vector<double> rhie_chow_face_flux(
    const PolyMesh& mesh,
    std::span<const Vec3> h_by_a,
    std::span<const double> pressure_mobility,
    std::span<const double> pressure,
    std::span<const VelocityBoundaryCondition> velocity_boundary,
    std::span<const PressureBoundaryCondition> pressure_boundary,
    bool include_nonorthogonal_correction = true);

[[nodiscard]] double flux_divergence_l2(
    const PolyMesh& mesh,
    std::span<const double> face_flux);

[[nodiscard]] double flux_divergence_max(
    const PolyMesh& mesh,
    std::span<const double> face_flux);

} // namespace cfd::fvm
