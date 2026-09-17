#pragma once
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::optics {
struct WavefrontSample {double rho{},theta{},opd{};};
struct ZernikeCoefficient {int n{},m{};double value{};};
[[nodiscard]] double zernike(int n,int m,double rho,double theta);
[[nodiscard]] std::vector<ZernikeCoefficient> fit_zernike(std::span<const WavefrontSample> samples,unsigned max_radial_order);
[[nodiscard]] double wavefront_rms(std::span<const WavefrontSample> samples,bool remove_piston=true);
[[nodiscard]] double longitudinal_chromatic_shift(std::span<const double> focal_positions);
}
