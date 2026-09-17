#pragma once

#include "cfd/core/aligned_allocator.hpp"

#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

namespace cfd::fdtd {

enum class Boundary1D { pec, mur1, pml, cpml };
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
    double cpml_kappa_max{5.0};
    double cpml_alpha_fraction{0.0};
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

    // +x total-field/scattered-field interface. The analytical incident E field
    // is sampled at the Yee E/H half steps and injected through equivalent
    // corrections at the left edge of the total-field region.
    void set_tfsf_source(std::size_t first_total_cell,std::function<double(double)> incident_e,
                         double impedance_ohm=0.0);
    void clear_tfsf_source() noexcept;

    // Parallel RLC element represented by a field-aligned voltage V=E*length
    // and terminal current distributed over area. Infinity disables R or L.
    void set_parallel_lumped_rlc(std::size_t cell,double length_m,double area_m2,
                                 double resistance_ohm=std::numeric_limits<double>::infinity(),
                                 double inductance_h=std::numeric_limits<double>::infinity(),
                                 double capacitance_f=0.0);
    void clear_lumped_rlc(std::size_t cell);

    void add_soft_source(std::size_t cell,double value);
    void set_hard_source(std::size_t cell,double value);
    void step(std::size_t count = 1);

    [[nodiscard]] double energy() const;
    [[nodiscard]] const cfd::core::AlignedVector<double>& electric() const noexcept { return ez_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& magnetic() const noexcept { return hy_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& polarization() const noexcept { return polarization_; }
    [[nodiscard]] DispersionModel1D dispersion_model(std::size_t cell) const;
    [[nodiscard]] double lumped_inductor_current_density(std::size_t cell) const;
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

    cfd::core::AlignedVector<double> cpml_kappa_e_,cpml_b_e_,cpml_c_e_,cpml_psi_e_;
    cfd::core::AlignedVector<double> cpml_kappa_h_,cpml_b_h_,cpml_c_h_,cpml_psi_h_;

    cfd::core::AlignedVector<double> lumped_epsilon_add_,lumped_sigma_add_;
    cfd::core::AlignedVector<double> lumped_inductor_drive_,lumped_inductor_current_density_;

    bool tfsf_enabled_{false};
    std::size_t tfsf_cell_{};
    std::function<double(double)> tfsf_incident_e_;
    double tfsf_impedance_{};

    void preprocess_material_coefficients();
    void initialize_pml_profile();
    void initialize_cpml_profile();
    void clear_dispersion(std::size_t begin_cell,std::size_t end_cell);
    void apply_boundary(double old_left,double old_left_adj,double old_right,double old_right_adj);
};

} // namespace cfd::fdtd
