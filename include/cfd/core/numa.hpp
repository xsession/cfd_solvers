#pragma once
#include "cfd/core/parallel.hpp"
#include <algorithm>
#include <cstddef>
#include <span>
namespace cfd::core {
// Portable NUMA first-touch baseline: pages are initialized by the worker that
// will typically consume that range. Binding threads to NUMA nodes remains OS/runtime specific.
template<class T> void numa_first_touch(std::span<T> values,const T& value=T{}){
    parallel_for(values.size(),[&](std::size_t i){values[i]=value;});
}
} // namespace cfd::core
