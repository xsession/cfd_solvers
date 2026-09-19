#pragma once

#include "cfd/optics/serialization.hpp"

#include <string>
#include <string_view>

namespace cfd::optics {

enum class LensInterchangeFormat { zemax_zmx, codev_seq, oslo_len };

// Clean-room text interchange for the common sequential surface subset shared
// by ZMX, CODE V SEQ and OSLO LEN files: surface type, radius/curvature,
// vertex spacing, aperture, refractive index, conic constant and four even
// asphere coefficients. Unsupported vendor-specific records are ignored.
[[nodiscard]] std::string export_lens_system(const SequentialOpticalSystem& system, LensInterchangeFormat format);
[[nodiscard]] SequentialOpticalSystem import_lens_system(std::string_view text, LensInterchangeFormat format);

} // namespace cfd::optics
