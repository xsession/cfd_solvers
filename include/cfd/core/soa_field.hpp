#pragma once

#include "cfd/core/aligned_allocator.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>

namespace cfd::core {

template<class T, std::size_t Components>
class StaticSoA {
public:
    explicit StaticSoA(std::size_t size = 0) : size_(size), data_(size * Components) {}

    void resize(std::size_t size) {
        size_ = size;
        data_.resize(size * Components);
    }

    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] static constexpr std::size_t components() noexcept { return Components; }

    [[nodiscard]] T* component(std::size_t c) {
        if (c >= Components) throw std::out_of_range("StaticSoA component");
        return data_.data() + c * size_;
    }

    [[nodiscard]] const T* component(std::size_t c) const {
        if (c >= Components) throw std::out_of_range("StaticSoA component");
        return data_.data() + c * size_;
    }

    [[nodiscard]] T& operator()(std::size_t c, std::size_t i) noexcept {
        return data_[c * size_ + i];
    }

    [[nodiscard]] const T& operator()(std::size_t c, std::size_t i) const noexcept {
        return data_[c * size_ + i];
    }

    [[nodiscard]] AlignedVector<T>& raw() noexcept { return data_; }
    [[nodiscard]] const AlignedVector<T>& raw() const noexcept { return data_; }

private:
    std::size_t size_{};
    AlignedVector<T> data_;
};

} // namespace cfd::core
