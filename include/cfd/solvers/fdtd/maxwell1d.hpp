#pragma once

#include "cfd/core/aligned_allocator.hpp"

#include <cstddef>
#include <vector>

namespace cfd::fdtd {

enum class Boundary1D { pec, mur1, pml };
enum class DispersionModel1D : unsigned char { none, debye, drude, lorentz };

struct Maxwell1DConfig {
    std::size_t cells{1024};
    double dx{1.0e-3};
    double courant{0.99};
    double epsilon_r{1.0};
    double mu_r{1.0};
    Boundary1D boundary{Boundary1D::pec};
    std::size_t pml_cells{16U};
    double pml_order{3.0};
    double pml_target_reflection{1.0e-8};
};

class Maxwell1D {
public:
    explicit Maxwell1D(Maxwell1DConfig config);

    void initialize_gaussian(double center_fraction = 0.3, double width_fraction = 0.04);
    void set_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_r,double conductivity_s_per_m=0.0);
    void set_debye_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                            double delta_epsilon_r,double relaxation_time_s,double conductivity_s_per_m=0.0);
    void set_drude_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                            double plasma_frequency_rad_s,double collision_frequency_rad_s,double conductivity_s_per_m=0.0);
    void set_lorentz_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                              double delta_epsilon_r,double resonance_frequency_rad_s,
                              double damping_rad_s,double conductivity_s_per_m=0.0);
    void add_soft_source(std::size_t cell,double value);
    void set_hard_source(std::size_t cell,double value);
    void step(std::size_t count = 1);

    // Electromagnetic field energy; excludes energy stored in ADE oscillators.
    [[nodiscard]] double energy() const;
    [[nodiscard]] const cfd::core::AlignedVector<double>& electric() const noexcept { return ez_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& magnetic() const noexcept { return hy_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& polarization() const noexcept { return polarization_; }
    [[nodiscard]] DispersionModel1D dispersion_model(std::size_t cell) const;
    [[nodiscard]] double dt() const noexcept { return dt_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] double local_wave_speed(std::size_t cell) const;

private:
    Maxwell1DConfig config_;
    double dt_{};
    double mu_{};
    double time_{};
    cfd::core::AlignedVector<double> epsilon_;
    cfd::core::AlignedVector<double> conductivity_;
    cfd::core::AlignedVector<double> e_decay_;
    cfd::core::AlignedVector<double> e_curl_;
    cfd::core::AlignedVector<double> pml_sigma_e_;
    cfd::core::AlignedVector<double> h_decay_;
    cfd::core::AlignedVector<double> h_curl_;
    cfd::core::AlignedVector<double> ez_;
    cfd::core::AlignedVector<double> hy_;
    std::vector<DispersionModel1D> dispersion_;
    cfd::core::AlignedVector<double> polarization_;
    cfd::core::AlignedVector<double> polarization_previous_;
    cfd::core::AlignedVector<double> dispersion_a_;
    cfd::core::AlignedVector<double> dispersion_b_;
    cfd::core::AlignedVector<double> dispersion_c_;

    void preprocess_material_coefficients();
    void initialize_pml_profile();
    void clear_dispersion(std::size_t begin_cell,std::size_t end_cell);
    void validate_material(std::size_t begin_cell,std::size_t end_cell,
                           double epsilon_r,double conductivity,bool dispersive) const;
    void apply_boundary(double old_left,double old_left_adj,double old_right,double old_right_adj);
};

} // namespace cfd::fdtd
