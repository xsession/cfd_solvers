#pragma once
#include "cfd/fvm/schemes.hpp"
#include <string_view>
namespace cfd::fvm {
[[nodiscard]] FaceInterpolationScheme parse_face_interpolation_scheme(std::string_view name);
[[nodiscard]] const char* face_interpolation_scheme_name(FaceInterpolationScheme scheme) noexcept;
[[nodiscard]] double courant_limited_timestep(const PolyMesh& mesh,std::span<const double> face_flux,double target_courant,double maximum_dt);
} // namespace cfd::fvm
