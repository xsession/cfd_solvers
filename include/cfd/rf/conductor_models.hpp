#pragma once

#include "cfd/rf/antenna.hpp"
#include "cfd/rf/peec.hpp"

#include <complex>
#include <cstddef>
#include <vector>

namespace cfd::rf {

struct GroundReflection {
    std::complex<double> TE{};
    std::complex<double> TM{};
};

// Fresnel reflection coefficients for a lossy planar ground half-space. This
// is a reusable finite/conductive-ground boundary contract; the thin-wire MoM
// image kernel remains explicitly PEC/free-space today.
[[nodiscard]] GroundReflection finite_ground_reflection(double frequency_hz, double conductivity_s_per_m,
                                                        double relative_permittivity, double incidence_angle_rad);

// Engineering current-crowding factor for two parallel round conductors.
// It is deliberately a bounded PEEC correction, not a replacement for a
// volumetric skin/proximity field solve.
[[nodiscard]] double round_wire_proximity_factor(double center_spacing_m, double radius_m, double frequency_hz,
                                                 double conductivity_s_per_m, double relative_permeability = 1.0);

// Cross-section subdivision used by surface/volume PEEC callers. The source
// segment is assumed to run along x; returned filaments tile its y-z section.
[[nodiscard]] std::vector<FilamentSegment> subdivide_rectangular_peec(const FilamentSegment& centerline,
                                                                      double width_y_m, double height_z_m,
                                                                      std::size_t divisions_y, std::size_t divisions_z);

} // namespace cfd::rf
