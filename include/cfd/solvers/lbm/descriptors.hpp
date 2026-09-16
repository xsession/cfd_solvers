#pragma once

#include <array>
#include <cstddef>

namespace cfd::lbm {

namespace detail {

template<class Descriptor>
consteval bool descriptor_pairs_are_opposites() {
    if (Descriptor::q < 3 || (Descriptor::q % 2) == 0) return false;
    if (Descriptor::cx(0) != 0 || Descriptor::cy(0) != 0 || Descriptor::cz(0) != 0) return false;
    for (int i = 1; i < Descriptor::q; i += 2) {
        if (Descriptor::cx(i) != -Descriptor::cx(i + 1)) return false;
        if (Descriptor::cy(i) != -Descriptor::cy(i + 1)) return false;
        if (Descriptor::cz(i) != -Descriptor::cz(i + 1)) return false;
    }
    return true;
}

template<class Descriptor>
consteval bool descriptor_weights_sum_to_one() {
    double sum = 0.0;
    for (int i = 0; i < Descriptor::q; ++i) sum += static_cast<double>(Descriptor::weight(i));
    return sum > 0.999999 && sum < 1.000001;
}

} // namespace detail

struct D2Q9InPlaceDescriptor {
    static constexpr int dimensions = 2;
    static constexpr int q = 9;

    inline static constexpr std::array<int, q> cx_values{
        0,  1, -1,  0,  0,  1, -1, -1,  1
    };
    inline static constexpr std::array<int, q> cy_values{
        0,  0,  0,  1, -1,  1, -1,  1, -1
    };
    inline static constexpr std::array<int, q> cz_values{};
    inline static constexpr std::array<float, q> weights{
        4.0F/9.0F,
        1.0F/9.0F, 1.0F/9.0F,
        1.0F/9.0F, 1.0F/9.0F,
        1.0F/36.0F, 1.0F/36.0F,
        1.0F/36.0F, 1.0F/36.0F
    };

    static constexpr int cx(int i) noexcept { return cx_values[static_cast<std::size_t>(i)]; }
    static constexpr int cy(int i) noexcept { return cy_values[static_cast<std::size_t>(i)]; }
    static constexpr int cz(int i) noexcept { return cz_values[static_cast<std::size_t>(i)]; }
    static constexpr float weight(int i) noexcept { return weights[static_cast<std::size_t>(i)]; }
    static constexpr int opposite(int i) noexcept {
        return i == 0 ? 0 : ((i & 1) != 0 ? i + 1 : i - 1);
    }
};

struct D3Q19Descriptor {
    static constexpr int dimensions = 3;
    static constexpr int q = 19;

    inline static constexpr std::array<int, q> cx_values{
        0,
         1, -1,  0,  0,  0,  0,
         1, -1,  1, -1,
         1, -1,  1, -1,
         0,  0,  0,  0
    };
    inline static constexpr std::array<int, q> cy_values{
        0,
         0,  0,  1, -1,  0,  0,
         1, -1, -1,  1,
         0,  0,  0,  0,
         1, -1,  1, -1
    };
    inline static constexpr std::array<int, q> cz_values{
        0,
         0,  0,  0,  0,  1, -1,
         0,  0,  0,  0,
         1, -1, -1,  1,
         1, -1, -1,  1
    };
    inline static constexpr std::array<float, q> weights{
        1.0F/3.0F,
        1.0F/18.0F, 1.0F/18.0F, 1.0F/18.0F, 1.0F/18.0F, 1.0F/18.0F, 1.0F/18.0F,
        1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F,
        1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F,
        1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F
    };

    static constexpr int cx(int i) noexcept { return cx_values[static_cast<std::size_t>(i)]; }
    static constexpr int cy(int i) noexcept { return cy_values[static_cast<std::size_t>(i)]; }
    static constexpr int cz(int i) noexcept { return cz_values[static_cast<std::size_t>(i)]; }
    static constexpr float weight(int i) noexcept { return weights[static_cast<std::size_t>(i)]; }
    static constexpr int opposite(int i) noexcept {
        return i == 0 ? 0 : ((i & 1) != 0 ? i + 1 : i - 1);
    }
};

struct D3Q27Descriptor {
    static constexpr int dimensions = 3;
    static constexpr int q = 27;

    inline static constexpr std::array<int, q> cx_values{
        0,
         1, -1,  0,  0,  0,  0,
         1, -1,  1, -1,
         1, -1,  1, -1,
         0,  0,  0,  0,
         1, -1,  1, -1,  1, -1,  1, -1
    };
    inline static constexpr std::array<int, q> cy_values{
        0,
         0,  0,  1, -1,  0,  0,
         1, -1, -1,  1,
         0,  0,  0,  0,
         1, -1,  1, -1,
         1, -1,  1, -1, -1,  1, -1,  1
    };
    inline static constexpr std::array<int, q> cz_values{
        0,
         0,  0,  0,  0,  1, -1,
         0,  0,  0,  0,
         1, -1, -1,  1,
         1, -1, -1,  1,
         1, -1, -1,  1,  1, -1, -1,  1
    };
    inline static constexpr std::array<float, q> weights{
        8.0F/27.0F,
        2.0F/27.0F, 2.0F/27.0F, 2.0F/27.0F, 2.0F/27.0F, 2.0F/27.0F, 2.0F/27.0F,
        1.0F/54.0F, 1.0F/54.0F, 1.0F/54.0F, 1.0F/54.0F,
        1.0F/54.0F, 1.0F/54.0F, 1.0F/54.0F, 1.0F/54.0F,
        1.0F/54.0F, 1.0F/54.0F, 1.0F/54.0F, 1.0F/54.0F,
        1.0F/216.0F, 1.0F/216.0F, 1.0F/216.0F, 1.0F/216.0F,
        1.0F/216.0F, 1.0F/216.0F, 1.0F/216.0F, 1.0F/216.0F
    };

    static constexpr int cx(int i) noexcept { return cx_values[static_cast<std::size_t>(i)]; }
    static constexpr int cy(int i) noexcept { return cy_values[static_cast<std::size_t>(i)]; }
    static constexpr int cz(int i) noexcept { return cz_values[static_cast<std::size_t>(i)]; }
    static constexpr float weight(int i) noexcept { return weights[static_cast<std::size_t>(i)]; }
    static constexpr int opposite(int i) noexcept {
        return i == 0 ? 0 : ((i & 1) != 0 ? i + 1 : i - 1);
    }
};

static_assert(detail::descriptor_pairs_are_opposites<D2Q9InPlaceDescriptor>());
static_assert(detail::descriptor_pairs_are_opposites<D3Q19Descriptor>());
static_assert(detail::descriptor_pairs_are_opposites<D3Q27Descriptor>());
static_assert(detail::descriptor_weights_sum_to_one<D2Q9InPlaceDescriptor>());
static_assert(detail::descriptor_weights_sum_to_one<D3Q19Descriptor>());
static_assert(detail::descriptor_weights_sum_to_one<D3Q27Descriptor>());

} // namespace cfd::lbm
