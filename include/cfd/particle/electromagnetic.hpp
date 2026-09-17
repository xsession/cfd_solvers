#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::particle {

struct Vec3 { double x{},y{},z{}; };
struct ElectromagneticField { Vec3 electric_v_per_m{},magnetic_t{}; };

struct ChargedParticle {
    Vec3 position_m{};
    Vec3 velocity_m_per_s{};
    double charge_c{1.0};
    double mass_kg{1.0};
    double weight{1.0};
};

using FieldSampler = std::function<ElectromagneticField(Vec3,double)>;

// Non-relativistic Boris pusher. It is volume preserving and exactly preserves
// speed for a uniform magnetic field in exact arithmetic.
void boris_push(ChargedParticle& particle,const ElectromagneticField& field,double dt_s);
void track_particle(ChargedParticle& particle,const FieldSampler& field,double start_time_s,
                    double dt_s,std::size_t steps);

struct PicParticle1D {
    double position_m{};
    double velocity_m_per_s{};
    double charge_c{1.0};
    double mass_kg{1.0};
    double weight{1.0};
};

struct ElectrostaticPic1DConfig {
    std::size_t grid_points{64U};
    double length_m{1.0};
    double relative_permittivity{1.0};
    double dt_s{1.0e-3};
    bool neutralize_mean_charge{true};
};

// Periodic spectral Poisson solve. rho is nodal charge density for a unit
// cross-sectional area. A non-zero mean is removed because periodic Poisson
// requires a neutral net charge.
[[nodiscard]] std::vector<double> periodic_electric_field_from_charge_density(
    std::span<const double> charge_density_c_per_m3,double length_m,double relative_permittivity,
    bool remove_mean = true);

class ElectrostaticPic1D {
public:
    explicit ElectrostaticPic1D(ElectrostaticPic1DConfig config = {});
    void set_particles(std::vector<PicParticle1D> particles);
    [[nodiscard]] const std::vector<PicParticle1D>& particles() const noexcept { return particles_; }
    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& electric_field() const noexcept { return electric_field_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void deposit_and_solve();
    void step(std::size_t steps = 1U);

private:
    ElectrostaticPic1DConfig config_;
    std::vector<PicParticle1D> particles_;
    std::vector<double> charge_density_;
    std::vector<double> electric_field_;
    double time_s_{};

    [[nodiscard]] double interpolate_field(double position_m) const;
    [[nodiscard]] double wrap_position(double position_m) const;
};

} // namespace cfd::particle
