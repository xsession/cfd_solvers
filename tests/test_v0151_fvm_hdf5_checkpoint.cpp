#include "cfd/io/fvm_checkpoint.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void validate_round_trip() {
    auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(3, 2, 1, 1.0, 0.5, 0.25, 0.17);
    std::vector<double> temperature(mesh.cell_count());
    std::vector<double> pressure(mesh.cell_count());
    for (std::size_t i = 0; i < mesh.cell_count(); ++i) {
        temperature[i] = 300.0 + static_cast<double>(i);
        pressure[i] = 101325.0 - 10.0 * static_cast<double>(i);
    }
    const std::vector<cfd::io::FvmScalarFieldView> fields{
        {"T", temperature}, {"p", pressure}
    };
    const auto path = std::filesystem::temp_directory_path() / "cfd_solvers_v0151_fvm_checkpoint.h5";
    std::filesystem::remove(path);
    cfd::io::write_fvm_hdf5_checkpoint(path, mesh, fields, 0.125, 17U);
    require(std::filesystem::exists(path), "checkpoint file should be created");

    const auto restored = cfd::io::read_fvm_hdf5_checkpoint(path);
    require(restored.mesh.cell_count() == mesh.cell_count(), "checkpoint cell count should round-trip");
    require(restored.mesh.face_count() == mesh.face_count(), "checkpoint face count should round-trip");
    require(restored.mesh.patches().size() == mesh.patches().size(), "checkpoint patches should round-trip");
    require(std::abs(restored.time - 0.125) < 1.0e-15 && restored.step == 17U,
            "checkpoint time and step should round-trip");
    require(restored.scalar_fields.size() == 2U, "checkpoint fields should round-trip");
    require(restored.scalar_fields[0].name == "T" && restored.scalar_fields[0].values == temperature,
            "temperature field should round-trip exactly");
    require(restored.scalar_fields[1].name == "p" && restored.scalar_fields[1].values == pressure,
            "pressure field should round-trip exactly");
    for (std::size_t i = 0; i < mesh.cell_count(); ++i) {
        const auto a = mesh.cells()[i];
        const auto b = restored.mesh.cells()[i];
        require(std::abs(a.center.x-b.center.x) < 1e-15 && std::abs(a.center.y-b.center.y) < 1e-15 &&
                std::abs(a.center.z-b.center.z) < 1e-15 && std::abs(a.volume-b.volume) < 1e-15,
                "cell geometry should round-trip");
    }
    std::filesystem::remove(path);
}

void validate_fail_closed() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(1, 1, 1);
    const std::vector<double> field{42.0};
    const std::vector<cfd::io::FvmScalarFieldView> fields{{"field", field}};
    bool threw = false;
    try {
        cfd::io::write_fvm_hdf5_checkpoint("unused.h5", mesh, fields, 0.0, 0U);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    require(threw, "non-HDF5 build must fail closed for checkpoint output");
}
} // namespace

int main() {
    try {
        if (cfd::io::fvm_hdf5_checkpoint_available()) validate_round_trip();
        else validate_fail_closed();
        std::cout << "v0.15.1 FVM HDF5 checkpoint tests passed (available="
                  << cfd::io::fvm_hdf5_checkpoint_available() << ")\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v0.15.1 FVM HDF5 checkpoint test failed: " << error.what() << '\n';
        return 1;
    }
}
