#pragma once

#include <cstddef>
#include <new>
#include <vector>

namespace cfd::core {

template<class T, std::size_t Alignment = 64>
class AlignedAllocator {
public:
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template<class U>
    struct rebind { using other = AlignedAllocator<U, Alignment>; };

    AlignedAllocator() noexcept = default;
    template<class U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    [[nodiscard]] T* allocate(std::size_t n) {
        if (n > static_cast<std::size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length{};
        }
        return static_cast<T*>(::operator new(n * sizeof(T), std::align_val_t{Alignment}));
    }

    void deallocate(T* ptr, std::size_t) noexcept {
        ::operator delete(ptr, std::align_val_t{Alignment});
    }
};

template<class T, class U, std::size_t Alignment>
constexpr bool operator==(const AlignedAllocator<T, Alignment>&,
                          const AlignedAllocator<U, Alignment>&) noexcept {
    return true;
}

template<class T, class U, std::size_t Alignment>
constexpr bool operator!=(const AlignedAllocator<T, Alignment>&,
                          const AlignedAllocator<U, Alignment>&) noexcept {
    return false;
}

template<class T, std::size_t Alignment = 64>
using AlignedVector = std::vector<T, AlignedAllocator<T, Alignment>>;

} // namespace cfd::core
