#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/distributed/cartesian.hpp"
#include "cfd/solvers/lbm/distributed_pull.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace cfd::lbm {

#if defined(CFD_HAS_SYCL)

template<class Descriptor>
class DistributedSyclPullBlock {
public:
    static constexpr int q = Descriptor::q;

    DistributedSyclPullBlock(DistributedLbmConfig config,
                             cfd::distributed::Brick brick,
                             const sycl::device& device)
        : config_(config), brick_(brick),
          total_{brick.extent.x + 2U, brick.extent.y + 2U, brick.extent.z + 2U},
          total_cells_(total_.cells()),
          queue_(device, sycl::property::queue::in_order{}) {
        validate();
        current_ = sycl::malloc_device<float>(static_cast<std::size_t>(q) * total_cells_, queue_);
        next_ = sycl::malloc_device<float>(static_cast<std::size_t>(q) * total_cells_, queue_);
        if (!current_ || !next_) {
            if (current_) sycl::free(current_, queue_);
            if (next_) sycl::free(next_, queue_);
            current_ = nullptr;
            next_ = nullptr;
            throw std::bad_alloc{};
        }
        initialize_uniform();
    }

    ~DistributedSyclPullBlock() noexcept {
        try { queue_.wait_and_throw(); } catch (...) {}
        if (current_) sycl::free(current_, queue_);
        if (next_) sycl::free(next_, queue_);
    }

    DistributedSyclPullBlock(const DistributedSyclPullBlock&) = delete;
    DistributedSyclPullBlock& operator=(const DistributedSyclPullBlock&) = delete;

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
        std::vector<float> host(static_cast<std::size_t>(q) * total_cells_, 0.0F);
        const auto e = brick_.extent;
        for (std::size_t z = 1; z <= e.z; ++z) {
            for (std::size_t y = 1; y <= e.y; ++y) {
                for (std::size_t x = 1; x <= e.x; ++x) {
                    const auto state = initializer(brick_.begin.x + x - 1U,
                                                   brick_.begin.y + y - 1U,
                                                   brick_.begin.z + z - 1U);
                    const std::size_t cell = index(x, y, z);
                    for (int d = 0; d < q; ++d) {
                        host[static_cast<std::size_t>(d) * total_cells_ + cell] =
                            detail::equilibrium<Descriptor>(d, state[0], state[1], state[2], state[3]);
                    }
                }
            }
        }
        const std::size_t bytes = host.size() * sizeof(float);
        queue_.memcpy(current_, host.data(), bytes);
        queue_.memcpy(next_, host.data(), bytes).wait();
    }

    void compute_strict_interior() {
        launch_collision(false);
    }

    void compute_boundary() {
        launch_collision(true);
    }

    void compute_all() {
        launch_collision_all();
    }

    void finish_step() noexcept {
        std::swap(current_, next_);
        ++time_step_;
    }

    void wait() { queue_.wait_and_throw(); }

    [[nodiscard]] MacroscopicFields download_macroscopic() const {
        const auto host = download_storage();
        const std::size_t cells = brick_.extent.cells();
        MacroscopicFields fields{
            cfd::core::AlignedVector<float>(cells, 0.0F),
            cfd::core::AlignedVector<float>(cells, 0.0F),
            cfd::core::AlignedVector<float>(cells, 0.0F),
            cfd::core::AlignedVector<float>(cells, 0.0F)};
        std::size_t n = 0;
        const auto e = brick_.extent;
        for (std::size_t z = 1; z <= e.z; ++z) {
            for (std::size_t y = 1; y <= e.y; ++y) {
                for (std::size_t x = 1; x <= e.x; ++x, ++n) {
                    const std::size_t cell = index(x, y, z);
                    float rho = 0.0F;
                    float mx = 0.0F;
                    float my = 0.0F;
                    float mz = 0.0F;
                    for (int d = 0; d < q; ++d) {
                        const float value = host[static_cast<std::size_t>(d) * total_cells_ + cell];
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
                }
            }
        }
        return fields;
    }

    [[nodiscard]] double local_mass() const {
        const auto fields = download_macroscopic();
        double mass = 0.0;
        for (const float rho : fields.rho) mass += static_cast<double>(rho);
        return mass;
    }

    void save_checkpoint(const std::string& path) const {
        const auto host = download_storage();
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
                        write_value(out, host[static_cast<std::size_t>(d) * total_cells_ + index(x, y, z)]);
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
        std::vector<float> host(static_cast<std::size_t>(q) * total_cells_, 0.0F);
        const auto e = brick_.extent;
        for (int d = 0; d < q; ++d) {
            for (std::size_t z = 1; z <= e.z; ++z) {
                for (std::size_t y = 1; y <= e.y; ++y) {
                    for (std::size_t x = 1; x <= e.x; ++x) {
                        host[static_cast<std::size_t>(d) * total_cells_ + index(x, y, z)] = read_value<float>(in);
                    }
                }
            }
        }
        if (!in) throw std::runtime_error("checkpoint read failed");
        const std::size_t bytes = host.size() * sizeof(float);
        queue_.memcpy(current_, host.data(), bytes);
        queue_.memcpy(next_, host.data(), bytes).wait();
    }

    [[nodiscard]] const DistributedLbmConfig& config() const noexcept { return config_; }
    [[nodiscard]] const cfd::distributed::Brick& brick() const noexcept { return brick_; }
    [[nodiscard]] cfd::distributed::Extent3 total_extent() const noexcept { return total_; }
    [[nodiscard]] std::size_t total_cells() const noexcept { return total_cells_; }
    [[nodiscard]] float* current_data() noexcept { return current_; }
    [[nodiscard]] const float* current_data() const noexcept { return current_; }
    [[nodiscard]] sycl::queue& queue() noexcept { return queue_; }
    [[nodiscard]] const sycl::queue& queue() const noexcept { return queue_; }
    [[nodiscard]] std::uint64_t time_step() const noexcept { return time_step_; }

    [[nodiscard]] std::size_t index(std::size_t x, std::size_t y, std::size_t z) const noexcept {
        return (z * total_.y + y) * total_.x + x;
    }

private:
    DistributedLbmConfig config_{};
    cfd::distributed::Brick brick_{};
    cfd::distributed::Extent3 total_{};
    std::size_t total_cells_{};
    mutable sycl::queue queue_;
    float* current_{nullptr};
    float* next_{nullptr};
    std::uint64_t time_step_{};

    void validate() const {
        if (config_.tau <= 0.5F) throw std::invalid_argument("LBM tau must be > 0.5");
        if (brick_.global.x != config_.global.x || brick_.global.y != config_.global.y || brick_.global.z != config_.global.z) {
            throw std::invalid_argument("brick/global LBM extent mismatch");
        }
        if constexpr (Descriptor::dimensions == 2) {
            if (config_.global.z != 1 || brick_.extent.z != 1) throw std::invalid_argument("2-D distributed LBM requires z=1");
        }
    }

    void launch_collision(bool boundary_only) {
        const auto e = brick_.extent;
        const auto total = total_;
        const std::size_t total_cells = total_cells_;
        const std::size_t cells = e.cells();
        const float tau = config_.tau;
        const float ax = config_.acceleration_x;
        const float ay = config_.acceleration_y;
        const float az = config_.acceleration_z;
        const float* current = current_;
        float* next = next_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
            const std::size_t n = gid[0];
            const std::size_t lx = n % e.x;
            const std::size_t yz = n / e.x;
            const std::size_t ly = yz % e.y;
            const std::size_t lz = yz / e.y;
            const std::size_t x = lx + 1U;
            const std::size_t y = ly + 1U;
            const std::size_t z = lz + 1U;
            const bool strict_xy = x > 1U && x < e.x && y > 1U && y < e.y;
            const bool strict_z = Descriptor::dimensions == 2 || (z > 1U && z < e.z);
            const bool strict = strict_xy && strict_z;
            if (boundary_only == strict) return;
            collide_device_cell(current, next, total, total_cells, x, y, z, tau, ax, ay, az);
        });
    }

    void launch_collision_all() {
        const auto e = brick_.extent;
        const auto total = total_;
        const std::size_t total_cells = total_cells_;
        const std::size_t cells = e.cells();
        const float tau = config_.tau;
        const float ax = config_.acceleration_x;
        const float ay = config_.acceleration_y;
        const float az = config_.acceleration_z;
        const float* current = current_;
        float* next = next_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
            const std::size_t n = gid[0];
            const std::size_t x = n % e.x + 1U;
            const std::size_t yz = n / e.x;
            const std::size_t y = yz % e.y + 1U;
            const std::size_t z = yz / e.y + 1U;
            collide_device_cell(current, next, total, total_cells, x, y, z, tau, ax, ay, az);
        });
    }

    static void collide_device_cell(const float* current, float* next,
                                    cfd::distributed::Extent3 total, std::size_t total_cells,
                                    std::size_t x, std::size_t y, std::size_t z,
                                    float tau, float ax, float ay, float az) {
        float fin[Descriptor::q];
        float rho = 0.0F;
        float mx = 0.0F;
        float my = 0.0F;
        float mz = 0.0F;
        for (int d = 0; d < Descriptor::q; ++d) {
            const std::size_t xs = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(x) - Descriptor::cx(d));
            const std::size_t ys = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(y) - Descriptor::cy(d));
            const std::size_t zs = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(z) - Descriptor::cz(d));
            const std::size_t source = (zs * total.y + ys) * total.x + xs;
            const float value = current[static_cast<std::size_t>(d) * total_cells + source];
            fin[d] = value;
            rho += value;
            mx += value * static_cast<float>(Descriptor::cx(d));
            my += value * static_cast<float>(Descriptor::cy(d));
            mz += value * static_cast<float>(Descriptor::cz(d));
        }
        float ux = rho > 0.0F ? mx / rho : 0.0F;
        float uy = rho > 0.0F ? my / rho : 0.0F;
        float uz = rho > 0.0F ? mz / rho : 0.0F;
        ux += 0.5F * ax;
        uy += 0.5F * ay;
        uz += 0.5F * az;
        const float omega = 1.0F / tau;
        const float uu = 1.5F * (ux * ux + uy * uy + uz * uz);
        const float fx = rho * ax;
        const float fy = rho * ay;
        const float fz = rho * az;
        const float uf = ux * fx + uy * fy + uz * fz;
        const std::size_t cell = (z * total.y + y) * total.x + x;
        for (int d = 0; d < Descriptor::q; ++d) {
            const float cx = static_cast<float>(Descriptor::cx(d));
            const float cy = static_cast<float>(Descriptor::cy(d));
            const float cz = static_cast<float>(Descriptor::cz(d));
            const float cu_raw = cx * ux + cy * uy + cz * uz;
            const float cu = 3.0F * cu_raw;
            const float feq = Descriptor::weight(d) * rho * (1.0F + cu + 0.5F * cu * cu - uu);
            const float cf = cx * fx + cy * fy + cz * fz;
            const float force = Descriptor::weight(d) * (1.0F - 0.5F * omega) *
                                (3.0F * (cf - uf) + 9.0F * cu_raw * cf);
            next[static_cast<std::size_t>(d) * total_cells + cell] =
                fin[d] - omega * (fin[d] - feq) + force;
        }
    }

    [[nodiscard]] std::vector<float> download_storage() const {
        std::vector<float> host(static_cast<std::size_t>(q) * total_cells_);
        queue_.memcpy(host.data(), current_, host.size() * sizeof(float)).wait();
        return host;
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

#endif

} // namespace cfd::lbm
