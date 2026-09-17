#pragma once

#include "cfd/fem/mesh2d.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

namespace cfd::em {

using Complex = std::complex<double>;
using ComplexVec2 = std::array<Complex,2>;

struct NedelecTri3Geometry {
    double area{};
    // Physical gradients of barycentric coordinates lambda_0..lambda_2.
    std::array<std::array<double,2>,3> gradient_lambda{};
};

// Local lowest-order first-kind Nedelec edges are oriented (0->1),(1->2),(2->0).
[[nodiscard]] NedelecTri3Geometry nedelec_tri3_geometry(const cfd::fem::Mesh2D& mesh,
                                                        const cfd::fem::Tri3& triangle);
[[nodiscard]] std::array<double,2> nedelec_tri3_basis(const NedelecTri3Geometry& geometry,
                                                      std::size_t local_edge,
                                                      std::array<double,3> barycentric);
[[nodiscard]] double nedelec_tri3_curl(const NedelecTri3Geometry& geometry,std::size_t local_edge);

struct EdgeMaxwell2DConfig {
    double frequency_hz{1.0e8};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    double conductivity_s_per_m{};
};

struct EdgeMaxwell2DResult {
    // Global edge orientation is always lower node index -> higher node index.
    std::vector<std::array<std::size_t,2>> edges;
    // Line integral of E along each oriented global edge [V].
    std::vector<Complex> edge_voltage_v;
    // Reconstructed in-plane E at each triangle centroid [V/m].
    std::vector<ComplexVec2> electric_centroid_v_per_m;
};

using CurrentDensity2D = std::function<ComplexVec2(cfd::fem::Node2)>;

// Driven 2-D in-plane full-wave Maxwell reference using lowest-order Nedelec
// edge elements on Tri3 meshes. All mesh boundary edges are PEC (zero tangential
// electric field). The phasor convention is exp(+j*omega*t).
[[nodiscard]] EdgeMaxwell2DResult solve_driven_edge_maxwell_2d(
    const cfd::fem::Mesh2D& mesh,const EdgeMaxwell2DConfig& config,
    const CurrentDensity2D& impressed_current_density_a_per_m2);

} // namespace cfd::em
