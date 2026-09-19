#pragma once

#include "cfd/core/complex_sparse.hpp"
#include "cfd/fem/mesh3d.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::em {

using Complex = std::complex<double>;
using ComplexVec3 = std::array<Complex, 3>;
using RealVec3 = std::array<double, 3>;

struct NedelecTet4Geometry {
    double volume_m3{};
    std::array<RealVec3, 4> gradient_lambda{};
};

// Local first-kind lowest-order Nedelec edges:
// (0->1),(0->2),(0->3),(1->2),(1->3),(2->3).
[[nodiscard]] NedelecTet4Geometry nedelec_tet4_geometry(const cfd::fem::Mesh3D& mesh,
                                                        const cfd::fem::Tet4& tetrahedron);
[[nodiscard]] RealVec3 nedelec_tet4_basis(const NedelecTet4Geometry& geometry, std::size_t local_edge,
                                          std::array<double, 4> barycentric);
[[nodiscard]] RealVec3 nedelec_tet4_curl(const NedelecTet4Geometry& geometry, std::size_t local_edge);

struct EdgeMaxwell3DConfig {
    double frequency_hz{1.0e8};
    double relative_permittivity{1.0};
    double relative_permittivity_x{};
    double relative_permittivity_y{};
    double relative_permittivity_z{};
    double relative_permeability{1.0};
    double conductivity_s_per_m{};
    cfd::core::ComplexSparseSolveConfig linear{};
};

struct EdgeCurrentPort3D {
    // Directed mesh edge. The impressed phasor current is positive from node_a
    // to node_b and contributes through the edge basis line integral.
    std::size_t node_a{};
    std::size_t node_b{};
    Complex current_a{};
};

struct EdgeMaxwell3DResult {
    // Global edge orientation is lower node index -> higher node index.
    std::vector<std::array<std::size_t, 2>> edges;
    std::vector<unsigned char> boundary_edge;
    std::vector<Complex> edge_voltage_v;
    std::vector<ComplexVec3> electric_centroid_v_per_m;
    std::vector<ComplexVec3> magnetic_centroid_a_per_m;
    cfd::core::IterativeSolverResult linear_result{};
};

struct EdgePortThevenin {
    Complex open_circuit_voltage_v{};
    Complex port_impedance_ohm{};
};

using CurrentDensity3D = std::function<ComplexVec3(cfd::fem::Point3)>;

// Driven 3-D curl-curl Maxwell solve on Tet4 meshes with lowest-order Nedelec
// edge elements. Exterior mesh faces are PEC. The phasor convention is
// exp(+j*omega*t). A null volume source is allowed when edge ports are used.
[[nodiscard]] EdgeMaxwell3DResult
solve_driven_edge_maxwell_3d(const cfd::fem::Mesh3D& mesh, const EdgeMaxwell3DConfig& config,
                             const CurrentDensity3D& impressed_current_density_a_per_m2 = {},
                             std::span<const EdgeCurrentPort3D> edge_ports = {});

// Converts a solved edge-port field response to a small-signal Thevenin
// point. The caller can connect the returned impedance to Circuit/N-port
// termination code without coupling the FEM matrix to the circuit library.
[[nodiscard]] EdgePortThevenin edge_port_thevenin(const EdgeMaxwell3DResult& field, std::size_t edge_index,
                                                  Complex impressed_current_a);

struct ResonatorQuality3D {
    double electric_energy_j{};
    double magnetic_energy_j{};
    double dielectric_loss_w{};
    double conductor_loss_w{};
    double quality_factor{};
};

// Energy/loss post-processing for a harmonic 3-D edge-FEM field. Volume loss
// uses config.conductivity_s_per_m. Optional PEC-wall conductor loss uses the
// good-conductor surface resistance sqrt(pi*f*mu/sigma_wall).
[[nodiscard]] ResonatorQuality3D resonator_quality_3d(const cfd::fem::Mesh3D& mesh, const EdgeMaxwell3DConfig& config,
                                                      const EdgeMaxwell3DResult& field,
                                                      double wall_conductivity_s_per_m = 0.0,
                                                      double wall_relative_permeability = 1.0);

} // namespace cfd::em
