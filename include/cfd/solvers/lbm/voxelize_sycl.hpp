#pragma once

#include "cfd/solvers/lbm/voxelize.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

namespace cfd::lbm {
[[nodiscard]] std::vector<unsigned char> voxelize_surface_sycl(
    const cfd::io::TriangleSurface& surface,
    const VoxelGrid3D& grid,
    sycl::queue& queue);
} // namespace cfd::lbm
#endif
