#pragma once

#include <complex>
#include <vector>

namespace cfd::optics {

struct JonesVector {
    std::complex<double> x{1.0,0.0};
    std::complex<double> y{0.0,0.0};
};

struct StokesVector {
    double i{};
    double q{};
    double u{};
    double v{};
};

[[nodiscard]] StokesVector stokes_from_jones(JonesVector j);

struct FresnelCoefficients {
    double cos_transmitted{};
    double rs{};
    double rp{};
    double ts{};
    double tp{};
    double reflectance_s{};
    double reflectance_p{};
    bool total_internal_reflection{};
};

[[nodiscard]] FresnelCoefficients fresnel_dielectric(double n_from,double n_to,double cos_incident);

struct ThinFilmLayer {
    double refractive_index{};
    double thickness_nm{};
};

[[nodiscard]] double thin_film_reflectance_normal(double n_incident,
                                                   double n_substrate,
                                                   double wavelength_nm,
                                                   const std::vector<ThinFilmLayer>& layers);

} // namespace cfd::optics
