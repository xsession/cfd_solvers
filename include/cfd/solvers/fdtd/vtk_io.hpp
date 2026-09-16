#pragma once

#include <filesystem>

namespace cfd::fdtd {
class Maxwell3D;

// Portable legacy ASCII VTK output for debugging and ParaView interoperability.
// Fields are written on the solver's common padded logical lattice.
void write_maxwell3d_vtk_ascii(const Maxwell3D& solver,const std::filesystem::path& path);

} // namespace cfd::fdtd
