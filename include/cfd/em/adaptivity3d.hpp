#pragma once

#include "cfd/em/edge_fem3d.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::em {

// Face-jump estimator from reconstructed element fields. Tangential electric
// jumps and impedance-scaled tangential magnetic jumps are accumulated on the
// two tetrahedra adjacent to each interior face.
[[nodiscard]] std::vector<double> maxwell_face_jump_indicators_3d(
    const cfd::fem::Mesh3D& mesh,const EdgeMaxwell3DResult& field,
    double impedance_scale_ohm=376.730313668);

// Dorfler/bulk marking: smallest descending indicator prefix whose sum reaches
// fraction * total. Returns a byte mask over tetrahedra.
[[nodiscard]] std::vector<unsigned char> mark_maxwell_dorfler(
    std::span<const double> indicators,double fraction=0.5);

// Conforming local refinement by selecting the longest edge of every marked
// Tet4, then bisecting every tetrahedron incident on each selected global edge.
// Splitting the complete edge star prevents hanging nodes on shared faces.
[[nodiscard]] cfd::fem::Mesh3D refine_tet4_marked_longest_edges(
    const cfd::fem::Mesh3D& mesh,std::span<const unsigned char> marked);

} // namespace cfd::em
