#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <cstddef>
#include <filesystem>
#include <span>
#include <string>
#include <vector>

namespace cfd::io {

struct FvmScalarFieldView {
    std::string name;
    std::span<const double> values;
};

struct FvmScalarField {
    std::string name;
    std::vector<double> values;
};

struct FvmCheckpoint {
    cfd::fvm::PolyMesh mesh;
    double time{};
    std::size_t step{};
    std::vector<FvmScalarField> scalar_fields;
};

[[nodiscard]] bool fvm_hdf5_checkpoint_available() noexcept;

void write_fvm_hdf5_checkpoint(const std::filesystem::path& path,
                               const cfd::fvm::PolyMesh& mesh,
                               std::span<const FvmScalarFieldView> scalar_fields,
                               double time,
                               std::size_t step);

[[nodiscard]] FvmCheckpoint read_fvm_hdf5_checkpoint(const std::filesystem::path& path);

} // namespace cfd::io
