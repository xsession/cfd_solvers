#pragma once
#include <array>
#include <cstddef>
#include <type_traits>
namespace cfd::core {
template<class T,std::size_t N> struct SimdPack {static_assert(N>0);std::array<T,N> lane{};SimdPack()=default;explicit SimdPack(T x){lane.fill(x);}T& operator[](std::size_t i){return lane[i];}const T& operator[](std::size_t i)const{return lane[i];}
    friend SimdPack operator+(const SimdPack&a,const SimdPack&b){SimdPack r;for(std::size_t i=0;i<N;++i)r[i]=a[i]+b[i];return r;}friend SimdPack operator-(const SimdPack&a,const SimdPack&b){SimdPack r;for(std::size_t i=0;i<N;++i)r[i]=a[i]-b[i];return r;}friend SimdPack operator*(const SimdPack&a,const SimdPack&b){SimdPack r;for(std::size_t i=0;i<N;++i)r[i]=a[i]*b[i];return r;}friend SimdPack operator*(const SimdPack&a,T b){SimdPack r;for(std::size_t i=0;i<N;++i)r[i]=a[i]*b;return r;}[[nodiscard]] T sum()const{T r{};for(auto x:lane)r+=x;return r;}};
} // namespace cfd::core
