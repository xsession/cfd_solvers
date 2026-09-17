#pragma once
#include "cfd/fem/reference_element.hpp"
#include <complex>
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::fdtd {
struct ComplexVec3 {std::complex<double>x{},y{},z{};};
struct NearFieldSurfaceSample {cfd::fem::Point3 position{};cfd::fem::Point3 normal{};double area{};ComplexVec3 electric{},magnetic{};};
struct FarFieldSample {ComplexVec3 electric{};double radiation_intensity{};};
[[nodiscard]] FarFieldSample nf2ff(std::span<const NearFieldSurfaceSample> surface,cfd::fem::Point3 direction,double frequency_hz,double wave_impedance=376.730313668);
struct SarGrid3D {std::size_t nx{},ny{},nz{};double dx{},dy{},dz{};std::vector<double> ex,ey,ez,conductivity,density;void validate()const;};
[[nodiscard]] std::vector<double> local_sar(const SarGrid3D& grid);
[[nodiscard]] std::vector<double> mass_averaged_sar(const SarGrid3D& grid,double target_mass_kg);
[[nodiscard]] std::vector<double> sar_1g(const SarGrid3D& grid);
[[nodiscard]] std::vector<double> sar_10g(const SarGrid3D& grid);
} // namespace cfd::fdtd
