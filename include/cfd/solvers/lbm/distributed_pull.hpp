#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/distributed/cartesian.hpp"
#include "cfd/distributed/halo_grid.hpp"
#include "cfd/distributed/selective_halo.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace cfd::lbm {

struct DistributedLbmConfig {
    cfd::distributed::Extent3 global{128, 128, 128};
    float tau{0.6F};
    float acceleration_x{0.0F};
    float acceleration_y{0.0F};
    float acceleration_z{0.0F};
};

template<class Descriptor>
class DistributedPullBlock {
public:
    static constexpr int q = Descriptor::q;
    using Grid = cfd::distributed::HaloSoA3D<float, static_cast<std::size_t>(q)>;

    DistributedPullBlock(DistributedLbmConfig config, cfd::distributed::Brick brick)
        : config_(config), brick_(brick), current_(brick.extent), next_(brick.extent) {
        if (config_.tau <= 0.5F) throw std::invalid_argument("LBM tau must be > 0.5");
        if (brick_.global.x != config_.global.x || brick_.global.y != config_.global.y || brick_.global.z != config_.global.z) {
            throw std::invalid_argument("brick/global LBM extent mismatch");
        }
        if constexpr (Descriptor::dimensions == 2) {
            if (config_.global.z != 1 || brick_.extent.z != 1) throw std::invalid_argument("2-D distributed LBM requires z=1");
        }
        initialize_uniform();
    }

    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F, float uz = 0.0F) {
        initialize([=](std::size_t, std::size_t, std::size_t) {
            return std::array<float, 4>{rho, ux, uy, uz};
        });
    }

    void initialize_taylor_green(float amplitude = 0.03F) {
        const float two_pi = 2.0F * std::numbers::pi_v<float>;
        initialize([=, this](std::size_t gx, std::size_t gy, std::size_t gz) {
            const float xf = (static_cast<float>(gx) + 0.5F) / static_cast<float>(config_.global.x);
            const float yf = (static_cast<float>(gy) + 0.5F) / static_cast<float>(config_.global.y);
            const float sx = std::sin(two_pi * xf);
            const float cx = std::cos(two_pi * xf);
            const float sy = std::sin(two_pi * yf);
            const float cy = std::cos(two_pi * yf);
            float z_factor = 1.0F;
            if constexpr (Descriptor::dimensions == 3) {
                const float zf = (static_cast<float>(gz) + 0.5F) / static_cast<float>(config_.global.z);
                z_factor = std::cos(two_pi * zf);
            }
            return std::array<float, 4>{1.0F, amplitude * sx * cy * z_factor,
                                        -amplitude * cx * sy * z_factor, 0.0F};
        });
    }

    template<class Initializer>
    void initialize(Initializer&& initializer) {
        time_step_ = 0;
        const auto e = brick_.extent;
        for (std::size_t z = 1; z <= e.z; ++z) {
            for (std::size_t y = 1; y <= e.y; ++y) {
                for (std::size_t x = 1; x <= e.x; ++x) {
                    const auto state = initializer(brick_.begin.x + x - 1U,
                                                   brick_.begin.y + y - 1U,
                                                   brick_.begin.z + z - 1U);
                    for (int d = 0; d < q; ++d) {
                        const float value = detail::equilibrium<Descriptor>(d, state[0], state[1], state[2], state[3]);
                        current_(static_cast<std::size_t>(d), x, y, z) = value;
                        next_(static_cast<std::size_t>(d), x, y, z) = value;
                    }
                }
            }
        }
    }

    void compute_strict_interior() {
        const auto e = brick_.extent;
        if (e.x <= 2 || e.y <= 2 || (Descriptor::dimensions == 3 && e.z <= 2)) return;
        cfd::core::parallel_for(e.cells(), [&](std::size_t n) {
            const auto [x, y, z] = local_coordinates(n);
            const bool interior_xy = x > 1 && x < e.x && y > 1 && y < e.y;
            const bool interior_z = Descriptor::dimensions == 2 || (z > 1 && z < e.z);
            if (interior_xy && interior_z) collide_pull_cell(x, y, z);
        });
    }

    void compute_boundary() {
        const auto e = brick_.extent;
        cfd::core::parallel_for(e.cells(), [&](std::size_t n) {
            const auto [x, y, z] = local_coordinates(n);
            const bool interior_xy = x > 1 && x < e.x && y > 1 && y < e.y;
            const bool interior_z = Descriptor::dimensions == 2 || (z > 1 && z < e.z);
            if (!(interior_xy && interior_z)) collide_pull_cell(x, y, z);
        });
    }

    void compute_all() {
        const auto e = brick_.extent;
        cfd::core::parallel_for(e.cells(), [&](std::size_t n) {
            const auto [x, y, z] = local_coordinates(n);
            collide_pull_cell(x, y, z);
        });
    }

    void finish_step() noexcept {
        std::swap(current_, next_);
        ++time_step_;
    }

    [[nodiscard]] MacroscopicFields compute_macroscopic() const {
        const std::size_t cells = brick_.extent.cells();
        MacroscopicFields fields{
            cfd::core::AlignedVector<float>(cells, 0.0F),
            cfd::core::AlignedVector<float>(cells, 0.0F),
            cfd::core::AlignedVector<float>(cells, 0.0F),
            cfd::core::AlignedVector<float>(cells, 0.0F)};
        cfd::core::parallel_for(cells, [&](std::size_t n) {
            const auto [x, y, z] = local_coordinates(n);
            float rho = 0.0F;
            float mx = 0.0F;
            float my = 0.0F;
            float mz = 0.0F;
            for (int d = 0; d < q; ++d) {
                const float value = current_(static_cast<std::size_t>(d), x, y, z);
                rho += value;
                mx += value * static_cast<float>(Descriptor::cx(d));
                my += value * static_cast<float>(Descriptor::cy(d));
                mz += value * static_cast<float>(Descriptor::cz(d));
            }
            fields.rho[n] = rho;
            if (rho > 0.0F) {
                fields.ux[n] = mx / rho + 0.5F * config_.acceleration_x;
                fields.uy[n] = my / rho + 0.5F * config_.acceleration_y;
                fields.uz[n] = mz / rho + 0.5F * config_.acceleration_z;
            }
        });
        return fields;
    }

    [[nodiscard]] double local_mass() const {
        return cfd::core::parallel_sum(brick_.extent.cells(), [&](std::size_t n) {
            const auto [x, y, z] = local_coordinates(n);
            double rho = 0.0;
            for (int d = 0; d < q; ++d) rho += current_(static_cast<std::size_t>(d), x, y, z);
            return rho;
        });
    }

    [[nodiscard]] const DistributedLbmConfig& config() const noexcept { return config_; }
    [[nodiscard]] const cfd::distributed::Brick& brick() const noexcept { return brick_; }
    [[nodiscard]] Grid& current_grid() noexcept { return current_; }
    [[nodiscard]] const Grid& current_grid() const noexcept { return current_; }
    [[nodiscard]] std::uint64_t time_step() const noexcept { return time_step_; }

    void save_checkpoint(const std::string& path) const {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("cannot open checkpoint for writing: " + path);
        const std::array<char, 8> magic{'C','F','D','C','H','K','3','\0'};
        out.write(magic.data(), static_cast<std::streamsize>(magic.size()));
        write_value(out, std::uint32_t{1});
        write_value(out, std::uint32_t{0x01020304U});
        write_value(out, std::uint32_t{static_cast<std::uint32_t>(q)});
        write_extent(out, config_.global);
        write_extent(out, brick_.begin);
        write_extent(out, brick_.extent);
        write_value(out, time_step_);
        write_value(out, config_.tau);
        write_value(out, config_.acceleration_x);
        write_value(out, config_.acceleration_y);
        write_value(out, config_.acceleration_z);
        const auto e = brick_.extent;
        for (int d = 0; d < q; ++d) {
            for (std::size_t z = 1; z <= e.z; ++z) {
                for (std::size_t y = 1; y <= e.y; ++y) {
                    for (std::size_t x = 1; x <= e.x; ++x) {
                        write_value(out, current_(static_cast<std::size_t>(d), x, y, z));
                    }
                }
            }
        }
        if (!out) throw std::runtime_error("checkpoint write failed: " + path);
    }

    void load_checkpoint(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open checkpoint for reading: " + path);
        std::array<char, 8> magic{};
        in.read(magic.data(), static_cast<std::streamsize>(magic.size()));
        const std::array<char, 8> expected{'C','F','D','C','H','K','3','\0'};
        if (magic != expected) throw std::runtime_error("invalid checkpoint magic");
        const auto version = read_value<std::uint32_t>(in);
        const auto endian = read_value<std::uint32_t>(in);
        const auto stored_q = read_value<std::uint32_t>(in);
        if (version != 1 || endian != 0x01020304U || stored_q != static_cast<std::uint32_t>(q)) {
            throw std::runtime_error("unsupported checkpoint format/descriptor");
        }
        const auto global = read_extent(in);
        const auto begin = read_extent(in);
        const auto extent = read_extent(in);
        if (!same_extent(global, config_.global) || !same_extent(begin, brick_.begin) || !same_extent(extent, brick_.extent)) {
            throw std::runtime_error("checkpoint decomposition does not match block");
        }
        time_step_ = read_value<std::uint64_t>(in);
        const float tau = read_value<float>(in);
        const float ax = read_value<float>(in);
        const float ay = read_value<float>(in);
        const float az = read_value<float>(in);
        if (tau != config_.tau || ax != config_.acceleration_x || ay != config_.acceleration_y || az != config_.acceleration_z) {
            throw std::runtime_error("checkpoint LBM configuration mismatch");
        }
        const auto e = brick_.extent;
        for (int d = 0; d < q; ++d) {
            for (std::size_t z = 1; z <= e.z; ++z) {
                for (std::size_t y = 1; y <= e.y; ++y) {
                    for (std::size_t x = 1; x <= e.x; ++x) {
                        const float value = read_value<float>(in);
                        current_(static_cast<std::size_t>(d), x, y, z) = value;
                        next_(static_cast<std::size_t>(d), x, y, z) = value;
                    }
                }
            }
        }
        if (!in) throw std::runtime_error("checkpoint read failed");
    }

private:
    DistributedLbmConfig config_{};
    cfd::distributed::Brick brick_{};
    Grid current_;
    Grid next_;
    std::uint64_t time_step_{};

    [[nodiscard]] std::array<std::size_t, 3> local_coordinates(std::size_t n) const noexcept {
        const auto e = brick_.extent;
        const std::size_t x = n % e.x;
        const std::size_t yz = n / e.x;
        const std::size_t y = yz % e.y;
        const std::size_t z = yz / e.y;
        return {x + 1U, y + 1U, z + 1U};
    }

    void collide_pull_cell(std::size_t x, std::size_t y, std::size_t z) {
        std::array<float, q> fin{};
        float rho = 0.0F;
        float mx = 0.0F;
        float my = 0.0F;
        float mz = 0.0F;
        for (int d = 0; d < q; ++d) {
            const std::size_t xs = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(x) - Descriptor::cx(d));
            const std::size_t ys = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(y) - Descriptor::cy(d));
            const std::size_t zs = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(z) - Descriptor::cz(d));
            const float value = current_(static_cast<std::size_t>(d), xs, ys, zs);
            fin[static_cast<std::size_t>(d)] = value;
            rho += value;
            mx += value * static_cast<float>(Descriptor::cx(d));
            my += value * static_cast<float>(Descriptor::cy(d));
            mz += value * static_cast<float>(Descriptor::cz(d));
        }
        float ux = rho > 0.0F ? mx / rho : 0.0F;
        float uy = rho > 0.0F ? my / rho : 0.0F;
        float uz = rho > 0.0F ? mz / rho : 0.0F;
        ux += 0.5F * config_.acceleration_x;
        uy += 0.5F * config_.acceleration_y;
        uz += 0.5F * config_.acceleration_z;
        const float omega = 1.0F / config_.tau;
        for (int d = 0; d < q; ++d) {
            const std::size_t di = static_cast<std::size_t>(d);
            const float feq = detail::equilibrium<Descriptor>(d, rho, ux, uy, uz);
            next_(di, x, y, z) = fin[di] - omega * (fin[di] - feq) +
                                  detail::guo_force<Descriptor>(d, rho, ux, uy, uz,
                                                               config_.acceleration_x,
                                                               config_.acceleration_y,
                                                               config_.acceleration_z,
                                                               omega);
        }
    }

    template<class T>
    static void write_value(std::ostream& out, const T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        out.write(reinterpret_cast<const char*>(&value), static_cast<std::streamsize>(sizeof(T)));
    }
    static void write_extent(std::ostream& out, cfd::distributed::Extent3 e) {
        write_value(out, static_cast<std::uint64_t>(e.x));
        write_value(out, static_cast<std::uint64_t>(e.y));
        write_value(out, static_cast<std::uint64_t>(e.z));
    }
    template<class T>
    static T read_value(std::istream& in) {
        T value{};
        in.read(reinterpret_cast<char*>(&value), static_cast<std::streamsize>(sizeof(T)));
        return value;
    }
    static cfd::distributed::Extent3 read_extent(std::istream& in) {
        return {static_cast<std::size_t>(read_value<std::uint64_t>(in)),
                static_cast<std::size_t>(read_value<std::uint64_t>(in)),
                static_cast<std::size_t>(read_value<std::uint64_t>(in))};
    }
    static bool same_extent(cfd::distributed::Extent3 a, cfd::distributed::Extent3 b) noexcept {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

template<class Descriptor>
void virtual_exchange(std::vector<DistributedPullBlock<Descriptor>>& blocks,
                      std::array<bool, 3> periodic = {true, true, true}) {
    if (blocks.empty()) return;
    struct Message {
        std::size_t destination{};
        cfd::distributed::NeighborOffset incoming{};
        std::vector<float> payload;
    };
    std::vector<Message> messages;
    messages.reserve(blocks.size() * 26U);
    for (std::size_t rank = 0; rank < blocks.size(); ++rank) {
        auto& block = blocks[rank];
        if (block.brick().rank != rank) throw std::invalid_argument("virtual blocks must be ordered by rank");
        for (const auto offset : cfd::distributed::neighbor_offsets()) {
            const int neighbor = cfd::distributed::neighbor_rank(block.brick(), offset, periodic);
            if (neighbor < 0) continue;
            messages.push_back({static_cast<std::size_t>(neighbor), offset.opposite(),
                                block.current_grid().pack(offset)});
        }
    }
    for (auto& message : messages) {
        blocks[message.destination].current_grid().unpack(message.incoming, message.payload);
    }
}

template<class Descriptor>
void virtual_selective_exchange(std::vector<DistributedPullBlock<Descriptor>>& blocks,
                                std::array<bool, 3> periodic = {true, true, true}) {
    if (blocks.empty()) return;
    struct Message {
        std::size_t destination{};
        cfd::distributed::NeighborOffset incoming{};
        std::vector<float> payload;
    };
    std::vector<Message> messages;
    messages.reserve(blocks.size() * 18U);
    for (std::size_t rank = 0; rank < blocks.size(); ++rank) {
        auto& block = blocks[rank];
        if (block.brick().rank != rank) throw std::invalid_argument("virtual blocks must be ordered by rank");
        for (const auto offset : cfd::distributed::neighbor_offsets()) {
            if (cfd::distributed::selective_halo_value_count<Descriptor>(block.current_grid(), offset) == 0) continue;
            const int neighbor = cfd::distributed::neighbor_rank(block.brick(), offset, periodic);
            if (neighbor < 0) continue;
            messages.push_back({static_cast<std::size_t>(neighbor), offset.opposite(),
                                cfd::distributed::pack_selective<Descriptor>(block.current_grid(), offset)});
        }
    }
    for (auto& message : messages) {
        cfd::distributed::unpack_selective<Descriptor>(blocks[message.destination].current_grid(),
                                                       message.incoming, message.payload);
    }
}

template<class Descriptor>
void virtual_selective_step(std::vector<DistributedPullBlock<Descriptor>>& blocks,
                            std::array<bool, 3> periodic = {true, true, true}) {
    for (auto& block : blocks) block.compute_strict_interior();
    virtual_selective_exchange(blocks, periodic);
    for (auto& block : blocks) block.compute_boundary();
    for (auto& block : blocks) block.finish_step();
}

template<class Descriptor>
void virtual_step(std::vector<DistributedPullBlock<Descriptor>>& blocks,
                  std::array<bool, 3> periodic = {true, true, true}) {
    for (auto& block : blocks) block.compute_strict_interior();
    virtual_exchange(blocks, periodic);
    for (auto& block : blocks) block.compute_boundary();
    for (auto& block : blocks) block.finish_step();
}

template<class Descriptor>
[[nodiscard]] MacroscopicFields gather_virtual_macroscopic(const std::vector<DistributedPullBlock<Descriptor>>& blocks) {
    if (blocks.empty()) return {};
    const auto global = blocks.front().config().global;
    const std::size_t cells = global.cells();
    MacroscopicFields result{
        cfd::core::AlignedVector<float>(cells, 0.0F),
        cfd::core::AlignedVector<float>(cells, 0.0F),
        cfd::core::AlignedVector<float>(cells, 0.0F),
        cfd::core::AlignedVector<float>(cells, 0.0F)};
    for (const auto& block : blocks) {
        const auto macro = block.compute_macroscopic();
        const auto& b = block.brick();
        std::size_t n = 0;
        for (std::size_t z = 0; z < b.extent.z; ++z) {
            for (std::size_t y = 0; y < b.extent.y; ++y) {
                for (std::size_t x = 0; x < b.extent.x; ++x, ++n) {
                    const std::size_t gx = b.begin.x + x;
                    const std::size_t gy = b.begin.y + y;
                    const std::size_t gz = b.begin.z + z;
                    const std::size_t gi = (gz * global.y + gy) * global.x + gx;
                    result.rho[gi] = macro.rho[n];
                    result.ux[gi] = macro.ux[n];
                    result.uy[gi] = macro.uy[n];
                    result.uz[gi] = macro.uz[n];
                }
            }
        }
    }
    return result;
}

} // namespace cfd::lbm
