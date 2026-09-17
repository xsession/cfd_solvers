#pragma once
#include <complex>
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::optics {
struct DiffractionSample {double coordinate{};double intensity{};};
[[nodiscard]] std::vector<DiffractionSample> fraunhofer_1d(std::span<const std::complex<double>> pupil,double sample_pitch,double wavelength,double propagation_distance);
[[nodiscard]] std::vector<double> normalized_psf_1d(std::span<const std::complex<double>> pupil);
[[nodiscard]] std::vector<double> mtf_from_psf_1d(std::span<const double> psf);
} // namespace cfd::optics
