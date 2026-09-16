#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>

#if defined(CFD_HAS_OPENMP)
#include <omp.h>
#endif

namespace cfd::core {

inline unsigned cpu_thread_capacity() noexcept {
#if defined(CFD_HAS_OPENMP)
    return static_cast<unsigned>(omp_get_max_threads());
#else
    const unsigned n = std::thread::hardware_concurrency();
    return n == 0U ? 1U : n;
#endif
}

inline void set_cpu_threads(unsigned count) noexcept {
#if defined(CFD_HAS_OPENMP)
    if (count > 0U) {
        omp_set_num_threads(static_cast<int>(count));
    }
#else
    (void)count;
#endif
}

template<class Function>
inline void parallel_for(std::size_t count, Function&& fn) {
#if defined(CFD_HAS_OPENMP)
    #pragma omp parallel for schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(count); ++i) {
        fn(static_cast<std::size_t>(i));
    }
#else
    for (std::size_t i = 0; i < count; ++i) {
        fn(i);
    }
#endif
}

template<class Function>
inline double parallel_sum(std::size_t count, Function&& fn) {
    double value = 0.0;
#if defined(CFD_HAS_OPENMP)
    #pragma omp parallel for reduction(+:value) schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(count); ++i) {
        value += static_cast<double>(fn(static_cast<std::size_t>(i)));
    }
#else
    for (std::size_t i = 0; i < count; ++i) {
        value += static_cast<double>(fn(i));
    }
#endif
    return value;
}


template<class Function>
inline double parallel_max(std::size_t count, Function&& fn) {
    double value = 0.0;
#if defined(CFD_HAS_OPENMP)
    #pragma omp parallel for reduction(max:value) schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(count); ++i) {
        value = std::max(value, static_cast<double>(fn(static_cast<std::size_t>(i))));
    }
#else
    for (std::size_t i = 0; i < count; ++i) {
        value = std::max(value, static_cast<double>(fn(i)));
    }
#endif
    return value;
}

} // namespace cfd::core
