#pragma once

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::spectral {

using Complex = std::complex<double>;

[[nodiscard]] bool is_power_of_two(std::size_t n) noexcept;
void fft_inplace(std::span<Complex> values,bool inverse=false);

void fft2_inplace(std::span<Complex> values,std::size_t nx,std::size_t ny,bool inverse=false);

[[nodiscard]] std::vector<double> wavenumbers(std::size_t n,double spacing);
[[nodiscard]] std::vector<double> derivative_periodic_1d(std::span<const double> field,double spacing);
[[nodiscard]] std::vector<double> fractional_laplacian_periodic_1d(std::span<const double> field,double spacing,double order);

} // namespace cfd::spectral
