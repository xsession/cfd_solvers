#pragma once

#include "cfd/fem/reference_element.hpp"
#include "cfd/io/mesh_io.hpp"
#include "cfd/rf/nport.hpp"

#include <complex>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace cfd::rf {

// A triangular PEC surface element. The vertex order defines the local normal
// used by the RWG orientation builder.
struct SurfaceTriangle {
    fem::Point3 a{};
    fem::Point3 b{};
    fem::Point3 c{};
};

struct RwgBasisFunction {
    std::size_t edge{};
    fem::Point3 edge_start{};
    fem::Point3 edge_end{};
    fem::Point3 plus_free_vertex{};
    fem::Point3 minus_free_vertex{};
    std::size_t plus_triangle{};
    std::size_t minus_triangle{std::numeric_limits<std::size_t>::max()};
    double edge_length_m{};
    int plus_sign{};
    int minus_sign{};
};

enum class SurfaceMomExcitation { plane_wave, delta_gap };

struct SurfaceMomConfig {
    double frequency_hz{};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    SurfaceMomExcitation excitation{SurfaceMomExcitation::plane_wave};
    fem::Point3 plane_wave_direction{0.0, 0.0, 1.0};
    fem::Point3 plane_wave_electric{1.0, 0.0, 0.0};
    std::complex<double> feed_voltage_v{1.0, 0.0};
    std::size_t feed_edge{std::numeric_limits<std::size_t>::max()};
    // The compact reference assembler uses a three-point triangle rule and a
    // finite self-interaction distance. This is deliberately explicit: it is
    // useful for small studies, but is not a singular-quadrature or MLFMM claim.
    double singular_regularization_fraction{1.0e-3};
    double mesh_tolerance_m{1.0e-12};
};

struct SurfaceMomResult {
    std::vector<SurfaceTriangle> triangles;
    std::vector<RwgBasisFunction> basis;
    std::vector<std::complex<double>> current_a;
    ComplexMatrix system_matrix;
    std::size_t feed_edge{std::numeric_limits<std::size_t>::max()};
    std::complex<double> feed_current_a{};
    std::complex<double> feed_impedance_ohm{};
    double accepted_power_w{};
    double relative_residual{};
    double wave_number_per_m{};
    double wave_impedance_ohm{};
};

struct SurfaceMomFarField {
    fem::Point3 direction{};
    std::complex<double> e_x{};
    std::complex<double> e_y{};
    std::complex<double> e_z{};
    double power_density_w_m2{};
};

// Dense RWG electric-field integral-equation reference for a PEC triangle
// surface. It supports open surfaces, so boundary half-RWG functions are
// retained. The implementation is intended for small meshes and regression
// studies; dielectric SIE, singular quadrature, fast multipole acceleration,
// and electrically-large production meshes remain separate work.
[[nodiscard]] SurfaceMomResult solve_pec_surface_mom(std::span<const SurfaceTriangle> triangles,
                                                     const SurfaceMomConfig& config);

// Convenience bridge for the repository's OBJ/ASCII-STL triangle reader.
[[nodiscard]] SurfaceMomResult solve_pec_surface_mom(const io::TriangleSurface& surface,
                                                     const SurfaceMomConfig& config);

// Evaluate the radiated far field of a solved surface current using the
// electric-current radiation integral. The direction is normalized internally.
[[nodiscard]] SurfaceMomFarField surface_mom_far_field(const SurfaceMomResult& result,
                                                       fem::Point3 observation_direction, double distance_m);

} // namespace cfd::rf
