#include "cfd/io/mesh_io.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/lbm/visualization.hpp"
#include "cfd/solvers/lbm/voxelize.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        const auto model =
            cfd::io::read_ascii_stl(std::filesystem::path(CFD_SOURCE_DIR) / "examples" / "models" / "ahmed_body.stl");
        require(model.triangles.size() == 18U, "Ahmed body STL triangle count");

        constexpr std::size_t nx = 32U;
        constexpr std::size_t ny = 20U;
        constexpr std::size_t nz = 20U;
        const cfd::lbm::VoxelGrid3D grid{nx, ny, nz, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
        const auto solid = cfd::lbm::voxelize_surface(model, grid);
        const auto solid_cells = static_cast<std::size_t>(std::count(solid.begin(), solid.end(), 1U));
        require(solid_cells > 0U && solid_cells < solid.size(),
                "Ahmed body voxelization should produce a bounded solid");

        cfd::lbm::D3Q19Solver solver({nx, ny, nz, 0.58F});
        solver.initialize_taylor_green(0.02F);
        solver.step(4U);
        const auto fields = solver.compute_macroscopic();
        std::vector<cfd::lbm::Vector3> velocity(solver.cells());
        std::vector<double> speed(solver.cells());
        for (std::size_t i = 0; i < solver.cells(); ++i) {
            velocity[i] = {fields.ux[i], fields.uy[i], fields.uz[i]};
            speed[i] = std::sqrt(static_cast<double>(fields.ux[i]) * fields.ux[i] +
                                 static_cast<double>(fields.uy[i]) * fields.uy[i] +
                                 static_cast<double>(fields.uz[i]) * fields.uz[i]);
            require(std::isfinite(speed[i]), "3-D LBM model test produced non-finite speed");
        }
        const auto q = cfd::lbm::q_criterion_3d(nx, ny, nz, velocity);
        require(q.size() == solver.cells(), "3-D Q-criterion field size");

        constexpr std::size_t channel_nx = 20U;
        constexpr std::size_t channel_ny = 10U;
        constexpr std::size_t channel_nz = 8U;
        const auto channel_index = [](std::size_t x, std::size_t y, std::size_t z) {
            return (z * channel_ny + y) * channel_nx + x;
        };
        std::vector<std::uint8_t> channel_solid(channel_nx * channel_ny * channel_nz, 0U);
        for (std::size_t z = 0; z < channel_nz; ++z) {
            for (std::size_t x = 0; x < channel_nx; ++x) {
                channel_solid[channel_index(x, 0U, z)] = 1U;
                channel_solid[channel_index(x, channel_ny - 1U, z)] = 1U;
            }
        }
        cfd::lbm::D3Q19Solver channel({channel_nx, channel_ny, channel_nz, 0.8F, 2.0e-5F, 0.0F, 0.0F});
        channel.set_solid_mask(channel_solid);
        channel.initialize_uniform();
        channel.step(200U);
        const auto channel_fields = channel.compute_macroscopic();
        float fluid_max_ux = 0.0F;
        float solid_max_speed = 0.0F;
        for (std::size_t z = 0; z < channel_nz; ++z) {
            for (std::size_t y = 0; y < channel_ny; ++y) {
                for (std::size_t x = 0; x < channel_nx; ++x) {
                    const auto i = channel_index(x, y, z);
                    if (channel_solid[i] != 0U) {
                        solid_max_speed =
                            std::max(solid_max_speed, std::abs(channel_fields.ux[i]) + std::abs(channel_fields.uy[i]) +
                                                          std::abs(channel_fields.uz[i]));
                    } else {
                        fluid_max_ux = std::max(fluid_max_ux, channel_fields.ux[i]);
                    }
                }
            }
        }
        require(fluid_max_ux > 1.0e-4F, "3-D solid-wall channel responds to body acceleration");
        require(solid_max_speed < 1.0e-7F, "3-D stationary bounce-back keeps solid cells at rest");

        const auto out = std::filesystem::temp_directory_path() / "cfd_solvers_real_3d_model_test";
        std::filesystem::remove_all(out);
        cfd::lbm::write_vtk_structured_3d(out / "model.vtk", nx, ny, nz, speed, velocity, solid, "speed");
        cfd::lbm::write_json_state_3d(out / "state.json", nx, ny, nz, solver.time_step(), speed, solid, "speed");
        {
            std::ifstream vtk(out / "model.vtk");
            const std::string text((std::istreambuf_iterator<char>(vtk)), std::istreambuf_iterator<char>());
            require(text.find("DIMENSIONS 32 20 20") != std::string::npos &&
                        text.find("SCALARS solid unsigned_char 1") != std::string::npos,
                    "3-D VTK model export");
        }
        {
            std::ifstream state(out / "state.json");
            const std::string text((std::istreambuf_iterator<char>(state)), std::istreambuf_iterator<char>());
            require(text.find("\"dimensions\": [32, 20, 20]") != std::string::npos &&
                        text.find("\"solid\"") != std::string::npos,
                    "3-D live JSON export");
        }
        std::filesystem::remove_all(out);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "real 3-D model test failed: " << error.what() << '\n';
        return 1;
    }
}
