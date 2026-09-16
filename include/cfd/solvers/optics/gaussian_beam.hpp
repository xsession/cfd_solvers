#pragma once
#include <complex>

namespace cfd::optics {

// All lengths use metres. q=z+i*z_R; radius is the 1/e field (1/e^2 intensity) radius.
class GaussianBeam {
public:
    GaussianBeam(double vacuum_wavelength,double waist_radius,double distance_from_waist=0.0,
                 double refractive_index=1.0);
    void propagate(double distance);
    void thin_lens(double focal_length);
    // ABCD matrix acts on (height, angle); output index accounts for refraction.
    void transform(double a,double b,double c,double d,double output_refractive_index);
    [[nodiscard]] double radius() const;
    [[nodiscard]] double curvature_radius() const;
    [[nodiscard]] std::complex<double> q() const noexcept { return q_; }
    [[nodiscard]] double refractive_index() const noexcept { return index_; }
private:
    double wavelength_{};
    double index_{};
    std::complex<double> q_{};
};
} // namespace cfd::optics
