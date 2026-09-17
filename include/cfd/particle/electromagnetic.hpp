#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <vector>

#if defined(CFD_HAS_MPI)
#include <mpi.h>
#endif

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

// Relativistic Boris update using proper velocity u=gamma*v. It preserves
// gamma (kinetic energy) for a magnetic-only field up to roundoff.
void relativistic_boris_push(ChargedParticle& particle,const ElectromagneticField& field,double dt_s);

// Vay relativistic pusher for ultra-relativistic crossed E/B cases. It uses
// proper velocity internally and reduces the spurious E-cross-B drift error that
// ordinary relativistic Boris can accumulate at high gamma.
void vay_push(ChargedParticle& particle,const ElectromagneticField& field,double dt_s);

struct NeutralCollisionModel {
    double neutral_density_per_m3{};
    double elastic_cross_section_m2{};
    double ionization_cross_section_m2{};
    double ionization_energy_ev{15.0};
    unsigned long long random_seed{0xC0FFEEULL};
};
struct CollisionStatistics {
    std::size_t elastic_events{};
    std::size_t ionization_events{};
    double energy_loss_j{};
};
// Null-collision style Monte-Carlo collision baseline for neutral gas/plasma
// work. Elastic events isotropically scatter particle velocity at fixed kinetic
// energy. Ionization events remove the configured threshold energy and append a
// low-energy secondary macro-particle with the same charge/mass/weight.
[[nodiscard]] CollisionStatistics apply_monte_carlo_collisions(
    std::vector<ChargedParticle>& particles,double dt_s,const NeutralCollisionModel& model);
[[nodiscard]] double lorentz_gamma(const ChargedParticle& particle);
[[nodiscard]] double kinetic_energy_j(const ChargedParticle& particle);

struct AxisAlignedParticleBox { Vec3 minimum_m{},maximum_m{1.0,1.0,1.0}; };
enum class ParticleWallMode { absorb,specular_reflect };
struct SecondaryEmissionModel {
    double maximum_yield{1.5};
    double energy_at_maximum_ev{300.0};
    double threshold_energy_ev{20.0};
};
struct ParticleWallInteraction {
    bool alive{true};
    bool impacted{};
    Vec3 impact_normal{};
    double secondary_macro_weight{};
};
[[nodiscard]] double secondary_electron_yield(double incident_energy_ev,const SecondaryEmissionModel& model={});
// Applies an axis-aligned wall after a particle push. Absorbing walls report
// an expected emitted secondary macro-weight; reflection mirrors both position
// and the normal velocity component.
[[nodiscard]] ParticleWallInteraction apply_particle_box_boundary(
    ChargedParticle& particle,const AxisAlignedParticleBox& box,ParticleWallMode mode,
    const SecondaryEmissionModel* secondary=nullptr);

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


struct ChargeConservingCurrent1D {
    std::vector<double> charge_density_old_c_per_m3;
    std::vector<double> charge_density_new_c_per_m3;
    std::vector<double> current_x_a_per_m2;
    double continuity_linf_residual{};
};

// Deposits periodic CIC charge before and after a particle move and reconstructs
// the longitudinal current whose spectral divergence exactly satisfies the
// discrete continuity equation. The helper is intentionally independent from the
// field update so it can validate future higher-dimensional deposition schemes.
[[nodiscard]] ChargeConservingCurrent1D deposit_charge_conserving_current_1d(
    std::span<const ChargedParticle> particles_old,std::span<const ChargedParticle> particles_new,
    std::size_t grid_points,double length_m,double dt_s);

struct ElectromagneticPic1DConfig {
    std::size_t grid_points{64U};
    double length_m{1.0};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    double dt_s{1.0e-12};
    bool solve_longitudinal_poisson{true};
    bool neutralize_mean_charge{true};
};

struct ElectromagneticPicDiagnostics {
    double field_energy_j{};
    double particle_kinetic_energy_j{};
    double total_energy_j{};
    double charge_continuity_linf_residual{};
    CollisionStatistics collision_statistics{};
    std::size_t particle_count{};
};

// 1-D/3-velocity electromagnetic PIC baseline. The grid is periodic along x.
// Transverse fields use a Yee-style staggered update (Ey/Ez at nodes,
// By/Bz at half-cells), Ex is recovered from periodic electrostatic Poisson,
// particles are advanced with the Vay pusher, and longitudinal current is
// reconstructed from charge continuity.
class ElectromagneticPic1D {
public:
    explicit ElectromagneticPic1D(ElectromagneticPic1DConfig config = {});

    void set_particles(std::vector<ChargedParticle> particles);
    [[nodiscard]] const std::vector<ChargedParticle>& particles() const noexcept { return particles_; }

    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& current_x() const noexcept { return current_x_; }
    [[nodiscard]] const std::vector<double>& current_y() const noexcept { return current_y_; }
    [[nodiscard]] const std::vector<double>& current_z() const noexcept { return current_z_; }
    [[nodiscard]] const std::vector<double>& electric_x() const noexcept { return electric_x_; }
    [[nodiscard]] const std::vector<double>& electric_y() const noexcept { return electric_y_; }
    [[nodiscard]] const std::vector<double>& electric_z() const noexcept { return electric_z_; }
    [[nodiscard]] const std::vector<double>& magnetic_y() const noexcept { return magnetic_y_; }
    [[nodiscard]] const std::vector<double>& magnetic_z() const noexcept { return magnetic_z_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void set_transverse_fields(std::span<const double> electric_y_v_per_m,
                               std::span<const double> electric_z_v_per_m,
                               std::span<const double> magnetic_y_t,
                               std::span<const double> magnetic_z_t);
    void initialize_right_traveling_mode(double electric_amplitude_v_per_m,std::size_t mode_index = 1U);
    void deposit_sources();
    [[nodiscard]] ElectromagneticField gather_field(double position_m) const;
    [[nodiscard]] ElectromagneticPicDiagnostics diagnostics() const;
    void step(std::size_t steps = 1U,const NeutralCollisionModel* collision_model = nullptr);

private:
    ElectromagneticPic1DConfig config_;
    std::vector<ChargedParticle> particles_;
    std::vector<double> charge_density_;
    std::vector<double> current_x_;
    std::vector<double> current_y_;
    std::vector<double> current_z_;
    std::vector<double> electric_x_;
    std::vector<double> electric_y_;
    std::vector<double> electric_z_;
    std::vector<double> magnetic_y_;
    std::vector<double> magnetic_z_;
    double time_s_{};
    double last_continuity_residual_{};
    CollisionStatistics accumulated_collisions_{};

    [[nodiscard]] double dx() const noexcept { return config_.length_m/static_cast<double>(config_.grid_points); }
    [[nodiscard]] double epsilon() const noexcept;
    [[nodiscard]] double mu() const noexcept;
    [[nodiscard]] double medium_light_speed() const noexcept;
    [[nodiscard]] double wrap_position(double position_m) const;
    [[nodiscard]] double interpolate_node_field(std::span<const double> nodal,double position_m) const;
    [[nodiscard]] double interpolate_half_field(std::span<const double> half_cell,double position_m) const;
    void advance_magnetic(double dt_s);
    void deposit_transverse_current(std::span<const ChargedParticle> particles_old);
};


struct Vec2 { double x{},y{}; };
struct PicParticle2D {
    Vec2 position_m{};
    Vec3 velocity_m_per_s{};
    double charge_c{1.0};
    double mass_kg{1.0};
    double weight{1.0};
};

struct ElectrostaticField2D {
    std::vector<double> electric_x_v_per_m;
    std::vector<double> electric_y_v_per_m;
};

struct ChargeConservingCurrent2D {
    std::vector<double> charge_density_old_c_per_m3;
    std::vector<double> charge_density_new_c_per_m3;
    std::vector<double> current_x_a_per_m2;
    std::vector<double> current_y_a_per_m2;
    double continuity_linf_residual{};
};

struct ElectrostaticPic2DConfig {
    std::size_t nx{32U};
    std::size_t ny{32U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double relative_permittivity{1.0};
    double dt_s{1.0e-3};
    bool neutralize_mean_charge{true};
};

[[nodiscard]] std::vector<double> deposit_cic_charge_density_2d(
    std::span<const PicParticle2D> particles,std::size_t nx,std::size_t ny,
    double length_x_m,double length_y_m);

[[nodiscard]] ElectrostaticField2D periodic_electric_field_from_charge_density_2d(
    std::span<const double> charge_density_c_per_m3,std::size_t nx,std::size_t ny,
    double length_x_m,double length_y_m,double relative_permittivity,bool remove_mean = true);

[[nodiscard]] ChargeConservingCurrent2D deposit_charge_conserving_current_2d(
    std::span<const PicParticle2D> particles_old,std::span<const PicParticle2D> particles_new,
    std::size_t nx,std::size_t ny,double length_x_m,double length_y_m,double dt_s);

class ElectrostaticPic2D {
public:
    explicit ElectrostaticPic2D(ElectrostaticPic2DConfig config = {});
    void set_particles(std::vector<PicParticle2D> particles);
    [[nodiscard]] const std::vector<PicParticle2D>& particles() const noexcept { return particles_; }
    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& electric_x() const noexcept { return electric_x_; }
    [[nodiscard]] const std::vector<double>& electric_y() const noexcept { return electric_y_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void deposit_and_solve();
    [[nodiscard]] Vec2 gather_electric_field(Vec2 position_m) const;
    void step(std::size_t steps = 1U);

private:
    ElectrostaticPic2DConfig config_;
    std::vector<PicParticle2D> particles_;
    std::vector<double> charge_density_;
    std::vector<double> electric_x_;
    std::vector<double> electric_y_;
    double time_s_{};

    [[nodiscard]] double dx() const noexcept { return config_.length_x_m/static_cast<double>(config_.nx); }
    [[nodiscard]] double dy() const noexcept { return config_.length_y_m/static_cast<double>(config_.ny); }
    [[nodiscard]] Vec2 wrap_position(Vec2 position_m) const;
    [[nodiscard]] double bilinear(std::span<const double> nodal,Vec2 position_m) const;
};



struct ElectromagneticPic2DConfig {
    std::size_t nx{32U};
    std::size_t ny{32U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    double dt_s{1.0e-12};
    bool solve_longitudinal_poisson{false};
    bool neutralize_mean_charge{true};
};

struct ElectromagneticPic2DDiagnostics {
    double field_energy_j{};
    double particle_kinetic_energy_j{};
    double total_energy_j{};
    double charge_continuity_linf_residual{};
    CollisionStatistics collision_statistics{};
    std::size_t particle_count{};
};

// Periodic 2-D/3-velocity electromagnetic PIC baseline. All field components
// are stored on a periodic nodal grid and updated with centered curl operators;
// particles gather bilinear E/B, advance with the Vay pusher, and deposit in-
// plane charge-conserving Jx/Jy plus CIC Jz for transverse current coupling.
class ElectromagneticPic2D {
public:
    explicit ElectromagneticPic2D(ElectromagneticPic2DConfig config = {});

    void set_particles(std::vector<PicParticle2D> particles);
    [[nodiscard]] const std::vector<PicParticle2D>& particles() const noexcept { return particles_; }

    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& current_x() const noexcept { return current_x_; }
    [[nodiscard]] const std::vector<double>& current_y() const noexcept { return current_y_; }
    [[nodiscard]] const std::vector<double>& current_z() const noexcept { return current_z_; }
    [[nodiscard]] const std::vector<double>& electric_x() const noexcept { return electric_x_; }
    [[nodiscard]] const std::vector<double>& electric_y() const noexcept { return electric_y_; }
    [[nodiscard]] const std::vector<double>& electric_z() const noexcept { return electric_z_; }
    [[nodiscard]] const std::vector<double>& magnetic_x() const noexcept { return magnetic_x_; }
    [[nodiscard]] const std::vector<double>& magnetic_y() const noexcept { return magnetic_y_; }
    [[nodiscard]] const std::vector<double>& magnetic_z() const noexcept { return magnetic_z_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void set_fields(std::span<const double> electric_x_v_per_m,
                    std::span<const double> electric_y_v_per_m,
                    std::span<const double> electric_z_v_per_m,
                    std::span<const double> magnetic_x_t,
                    std::span<const double> magnetic_y_t,
                    std::span<const double> magnetic_z_t);
    void initialize_tm_z_mode(double electric_z_amplitude_v_per_m,std::size_t mode_x = 1U,std::size_t mode_y = 0U);
    void deposit_sources();
    [[nodiscard]] ElectromagneticField gather_field(Vec2 position_m) const;
    [[nodiscard]] ElectromagneticPic2DDiagnostics diagnostics() const;
    void step(std::size_t steps = 1U,const NeutralCollisionModel* collision_model = nullptr);

private:
    ElectromagneticPic2DConfig config_;
    std::vector<PicParticle2D> particles_;
    std::vector<double> charge_density_;
    std::vector<double> current_x_;
    std::vector<double> current_y_;
    std::vector<double> current_z_;
    std::vector<double> electric_x_;
    std::vector<double> electric_y_;
    std::vector<double> electric_z_;
    std::vector<double> magnetic_x_;
    std::vector<double> magnetic_y_;
    std::vector<double> magnetic_z_;
    double time_s_{};
    double last_continuity_residual_{};
    CollisionStatistics accumulated_collisions_{};

    [[nodiscard]] double dx() const noexcept { return config_.length_x_m/static_cast<double>(config_.nx); }
    [[nodiscard]] double dy() const noexcept { return config_.length_y_m/static_cast<double>(config_.ny); }
    [[nodiscard]] double epsilon() const noexcept;
    [[nodiscard]] double mu() const noexcept;
    [[nodiscard]] double medium_light_speed() const noexcept;
    [[nodiscard]] Vec2 wrap_position(Vec2 position_m) const;
    [[nodiscard]] double bilinear(std::span<const double> nodal,Vec2 position_m) const;
    [[nodiscard]] double ddx(std::span<const double> values,std::size_t ix,std::size_t iy) const;
    [[nodiscard]] double ddy(std::span<const double> values,std::size_t ix,std::size_t iy) const;
    void advance_magnetic(double dt_s);
    void deposit_transverse_current_z(std::span<const PicParticle2D> particles_old);
};

struct PicBoundaryBox2D { Vec2 minimum_m{},maximum_m{1.0,1.0}; };
struct ParticleBoundary2DResult {
    bool alive{true};
    bool impacted{};
    Vec2 impact_normal{};
    double secondary_macro_weight{};
};
[[nodiscard]] ParticleBoundary2DResult apply_particle_box_boundary_2d(
    PicParticle2D& particle,const PicBoundaryBox2D& box,ParticleWallMode mode,
    const SecondaryEmissionModel* secondary=nullptr);

enum class GridBoundaryMode2D { periodic,electric_wall,absorbing_sponge };
struct GridBoundary2DConfig {
    GridBoundaryMode2D mode{GridBoundaryMode2D::periodic};
    std::size_t sponge_cells{2U};
    double sponge_strength{3.0};
    double sponge_polynomial_order{2.0};
};
void apply_electrostatic_field_boundary_2d(std::vector<double>& electric_x_v_per_m,
                                           std::vector<double>& electric_y_v_per_m,
                                           std::size_t nx,std::size_t ny,
                                           const GridBoundary2DConfig& boundary);


struct StaggeredElectromagneticPic2DConfig {
    std::size_t nx{32U};
    std::size_t ny{32U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    double dt_s{1.0e-12};
    GridBoundary2DConfig field_boundary{};
    ParticleWallMode particle_boundary{ParticleWallMode::specular_reflect};
    bool periodic_particles{true};
    bool remove_absorbed_particles{true};
    bool enable_secondary_yield_report{false};
    SecondaryEmissionModel secondary_model{};
};

struct StaggeredElectromagneticPic2DDiagnostics {
    double field_energy_j{};
    double particle_kinetic_energy_j{};
    double total_energy_j{};
    double charge_continuity_linf_residual{};
    std::size_t particle_count{};
    std::size_t absorbed_particles{};
    double secondary_macro_weight{};
    CollisionStatistics collision_statistics{};
};

// True 2-D/3-V Yee-style electromagnetic PIC baseline. Electric and magnetic
// components live on staggered Yee locations but are exposed as same-sized
// arrays for simple diagnostics. Periodic boundaries use wrapped differences;
// electric-wall and sponge modes apply reusable non-periodic boundary handling.
class StaggeredElectromagneticPic2D {
public:
    explicit StaggeredElectromagneticPic2D(StaggeredElectromagneticPic2DConfig config = {});

    void set_particles(std::vector<PicParticle2D> particles);
    [[nodiscard]] const std::vector<PicParticle2D>& particles() const noexcept { return particles_; }

    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& current_x() const noexcept { return current_x_; }
    [[nodiscard]] const std::vector<double>& current_y() const noexcept { return current_y_; }
    [[nodiscard]] const std::vector<double>& current_z() const noexcept { return current_z_; }
    [[nodiscard]] const std::vector<double>& electric_x() const noexcept { return electric_x_; }
    [[nodiscard]] const std::vector<double>& electric_y() const noexcept { return electric_y_; }
    [[nodiscard]] const std::vector<double>& electric_z() const noexcept { return electric_z_; }
    [[nodiscard]] const std::vector<double>& magnetic_x() const noexcept { return magnetic_x_; }
    [[nodiscard]] const std::vector<double>& magnetic_y() const noexcept { return magnetic_y_; }
    [[nodiscard]] const std::vector<double>& magnetic_z() const noexcept { return magnetic_z_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void set_fields(std::span<const double> electric_x_v_per_m,
                    std::span<const double> electric_y_v_per_m,
                    std::span<const double> electric_z_v_per_m,
                    std::span<const double> magnetic_x_t,
                    std::span<const double> magnetic_y_t,
                    std::span<const double> magnetic_z_t);
    void initialize_tm_z_mode(double electric_z_amplitude_v_per_m,std::size_t mode_x = 1U,std::size_t mode_y = 0U);
    [[nodiscard]] ElectromagneticField gather_field(Vec2 position_m) const;
    [[nodiscard]] StaggeredElectromagneticPic2DDiagnostics diagnostics() const;
    void step(std::size_t steps = 1U,const NeutralCollisionModel* collision_model = nullptr);

private:
    StaggeredElectromagneticPic2DConfig config_;
    std::vector<PicParticle2D> particles_;
    std::vector<double> charge_density_;
    std::vector<double> current_x_;
    std::vector<double> current_y_;
    std::vector<double> current_z_;
    std::vector<double> electric_x_;
    std::vector<double> electric_y_;
    std::vector<double> electric_z_;
    std::vector<double> magnetic_x_;
    std::vector<double> magnetic_y_;
    std::vector<double> magnetic_z_;
    double time_s_{};
    double last_continuity_residual_{};
    std::size_t absorbed_particles_{};
    double secondary_macro_weight_{};
    CollisionStatistics accumulated_collisions_{};

    [[nodiscard]] double dx() const noexcept { return config_.length_x_m/static_cast<double>(config_.nx); }
    [[nodiscard]] double dy() const noexcept { return config_.length_y_m/static_cast<double>(config_.ny); }
    [[nodiscard]] double epsilon() const noexcept;
    [[nodiscard]] double mu() const noexcept;
    [[nodiscard]] double medium_light_speed() const noexcept;
    [[nodiscard]] bool periodic_fields() const noexcept;
    [[nodiscard]] Vec2 place_particle(Vec2 position_m) const;
    [[nodiscard]] double sample(std::span<const double> values,Vec2 position_m,double offset_x,double offset_y) const;
    [[nodiscard]] double diff_x_forward(std::span<const double> values,std::size_t ix,std::size_t iy) const;
    [[nodiscard]] double diff_y_forward(std::span<const double> values,std::size_t ix,std::size_t iy) const;
    [[nodiscard]] double diff_x_backward(std::span<const double> values,std::size_t ix,std::size_t iy) const;
    [[nodiscard]] double diff_y_backward(std::span<const double> values,std::size_t ix,std::size_t iy) const;
    void advance_magnetic(double dt_s);
    void deposit_transverse_current_z(std::span<const PicParticle2D> particles_old);
    void update_electric(double dt_s);
    void apply_field_boundary();
    void apply_particle_boundary();
};


struct PicParticle3D {
    Vec3 position_m{};
    Vec3 velocity_m_per_s{};
    double charge_c{1.0};
    double mass_kg{1.0};
    double weight{1.0};
};


struct ParticleCellSort3DResult {
    std::vector<PicParticle3D> sorted_particles;
    std::vector<std::size_t> original_indices;
    std::vector<std::size_t> cell_offsets;
    std::vector<std::size_t> cell_counts;
    std::size_t occupied_cells{};
};

[[nodiscard]] std::size_t particle_cell_index_3d(
    Vec3 position_m,std::size_t nx,std::size_t ny,std::size_t nz,
    double length_x_m,double length_y_m,double length_z_m);

// Stable cell ordering for CPU cache locality and future particle-domain
// decomposition. The helper keeps the original particle indices so diagnostics,
// restart and current-deposition code can correlate reordered particles with
// their source ordering.
[[nodiscard]] ParticleCellSort3DResult sort_particles_by_cell_3d(
    std::span<const PicParticle3D> particles,std::size_t nx,std::size_t ny,std::size_t nz,
    double length_x_m,double length_y_m,double length_z_m);

struct ParticleGuardHalo3DConfig {
    std::size_t nx{1U};
    std::size_t ny{1U};
    std::size_t nz{1U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double length_z_m{1.0};
    std::size_t guard_cells{1U};
    bool periodic{true};
};

struct ParticleGuardHalo3DResult {
    std::vector<std::size_t> x_minus;
    std::vector<std::size_t> x_plus;
    std::vector<std::size_t> y_minus;
    std::vector<std::size_t> y_plus;
    std::vector<std::size_t> z_minus;
    std::vector<std::size_t> z_plus;
    std::vector<std::size_t> interior;
    std::vector<std::size_t> outside_domain;
};

// Classifies particles near sub-domain faces. Edge/corner particles can appear
// in multiple face lists, matching the later guard-cell exchange requirement.
[[nodiscard]] ParticleGuardHalo3DResult classify_particle_guard_halos_3d(
    std::span<const PicParticle3D> particles,const ParticleGuardHalo3DConfig& config);

inline constexpr std::size_t invalid_pic_domain_id_3d = static_cast<std::size_t>(-1);

struct PicDomainGrid3DConfig {
    std::size_t global_nx{1U};
    std::size_t global_ny{1U};
    std::size_t global_nz{1U};
    std::size_t domains_x{1U};
    std::size_t domains_y{1U};
    std::size_t domains_z{1U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double length_z_m{1.0};
    bool periodic{true};
};

struct PicDomain3D {
    std::size_t domain_id{};
    std::size_t ix{};
    std::size_t iy{};
    std::size_t iz{};
    std::size_t first_x{};
    std::size_t first_y{};
    std::size_t first_z{};
    std::size_t cells_x{};
    std::size_t cells_y{};
    std::size_t cells_z{};
    Vec3 minimum_m{};
    Vec3 maximum_m{};
};

[[nodiscard]] std::vector<PicDomain3D> make_pic_domain_grid_3d(const PicDomainGrid3DConfig& config);
[[nodiscard]] std::size_t locate_pic_domain_3d(Vec3 position_m,const PicDomainGrid3DConfig& config);

struct ParticleDomainMigration3D {
    std::vector<std::vector<PicParticle3D>> particles_by_domain;
    std::vector<std::vector<std::size_t>> source_indices_by_domain;
    std::vector<std::size_t> destination_domain_by_particle;
    std::vector<std::size_t> outside_indices;
};

// Serial particle migration bucketing for future MPI domain decomposition.
// Periodic positions are wrapped before locating the destination domain;
// nonperiodic particles outside the global box are reported separately.
[[nodiscard]] ParticleDomainMigration3D plan_particle_domain_migration_3d(
    std::span<const PicParticle3D> particles,const PicDomainGrid3DConfig& config);

struct ScalarGuardedBlock3D {
    std::size_t domain_id{};
    std::size_t interior_nx{};
    std::size_t interior_ny{};
    std::size_t interior_nz{};
    std::size_t guard_cells{};
    std::vector<double> values;

    [[nodiscard]] std::size_t padded_nx() const noexcept { return interior_nx + 2U*guard_cells; }
    [[nodiscard]] std::size_t padded_ny() const noexcept { return interior_ny + 2U*guard_cells; }
    [[nodiscard]] std::size_t padded_nz() const noexcept { return interior_nz + 2U*guard_cells; }
};

// Builds ghost-padded scalar blocks by copying guard cells from neighbouring
// serial subdomains. This is a deterministic pack/scatter contract for future
// MPI guard-cell exchange. Nonperiodic exterior ghost cells receive
// exterior_value.
[[nodiscard]] std::vector<ScalarGuardedBlock3D> exchange_scalar_guard_cells_3d(
    const std::vector<std::vector<double>>& cell_values_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    std::size_t guard_cells,double exterior_value = 0.0);


struct ElectromagneticFieldBlocks3D {
    std::vector<std::vector<double>> electric_x_by_domain;
    std::vector<std::vector<double>> electric_y_by_domain;
    std::vector<std::vector<double>> electric_z_by_domain;
    std::vector<std::vector<double>> magnetic_x_by_domain;
    std::vector<std::vector<double>> magnetic_y_by_domain;
    std::vector<std::vector<double>> magnetic_z_by_domain;
};

struct ElectromagneticGuardedFieldBlocks3D {
    std::vector<ScalarGuardedBlock3D> electric_x;
    std::vector<ScalarGuardedBlock3D> electric_y;
    std::vector<ScalarGuardedBlock3D> electric_z;
    std::vector<ScalarGuardedBlock3D> magnetic_x;
    std::vector<ScalarGuardedBlock3D> magnetic_y;
    std::vector<ScalarGuardedBlock3D> magnetic_z;
};

// Convenience wrapper around the scalar guard-cell exchange contract for the six
// electromagnetic field components used by the 3-D PIC kernels. It keeps the
// physics layer independent from MPI while giving future communicator/GPU code a
// single field-vector guard contract to validate.
[[nodiscard]] ElectromagneticGuardedFieldBlocks3D exchange_electromagnetic_field_guard_cells_3d(
    const ElectromagneticFieldBlocks3D& fields_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    std::size_t guard_cells,double exterior_value = 0.0);


struct ParticleMigrationRecord3D {
    std::size_t source_domain_id{};
    std::size_t source_index{};
};

struct ParticleMigrationMessage3D {
    std::size_t source_domain_id{};
    std::size_t destination_domain_id{};
    std::vector<PicParticle3D> particles;
    std::vector<std::size_t> source_indices;
};

struct PackedParticleMigration3D {
    std::vector<std::vector<PicParticle3D>> retained_by_domain;
    std::vector<std::vector<std::size_t>> retained_source_indices_by_domain;
    std::vector<ParticleMigrationMessage3D> messages;
    std::vector<ParticleMigrationRecord3D> outside_particles;
};

// Communication-ready particle migration contract. This does not require MPI:
// each message is an ordered payload that can later be sent from source_domain_id
// to destination_domain_id. Applying the messages to retained_by_domain must
// produce the same particle ownership as a serial migration bucket.
[[nodiscard]] PackedParticleMigration3D pack_particle_migration_messages_3d(
    const std::vector<std::vector<PicParticle3D>>& particles_by_source_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config);
[[nodiscard]] std::vector<std::vector<PicParticle3D>> apply_particle_migration_messages_3d(
    std::vector<std::vector<PicParticle3D>> retained_by_domain,
    std::span<const ParticleMigrationMessage3D> messages);

struct ScalarGuardCellValue3D {
    std::size_t padded_x{};
    std::size_t padded_y{};
    std::size_t padded_z{};
    double value{};
};

struct ScalarGuardCellMessage3D {
    std::size_t source_domain_id{};
    std::size_t destination_domain_id{};
    int offset_x{};
    int offset_y{};
    int offset_z{};
    std::vector<ScalarGuardCellValue3D> values;
};

// Communication-ready scalar guard exchange. Messages contain the destination
// padded coordinates to update, so an MPI layer only needs to transport the
// payload; it does not need to understand the global decomposition geometry.
[[nodiscard]] std::vector<ScalarGuardCellMessage3D> pack_scalar_guard_cell_messages_3d(
    const std::vector<std::vector<double>>& cell_values_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    std::size_t guard_cells);
[[nodiscard]] std::vector<ScalarGuardedBlock3D> apply_scalar_guard_cell_messages_3d(
    const std::vector<std::vector<double>>& cell_values_by_domain,
    std::span<const PicDomain3D> domains,std::size_t guard_cells,double exterior_value,
    std::span<const ScalarGuardCellMessage3D> messages);

enum class PicTransportPayloadKind3D { particle_migration, scalar_guard_cells };

struct PicRankTopology3D {
    std::size_t rank_count{};
    std::vector<std::size_t> domain_to_rank;
};

// Maps one or more PIC subdomains to deterministic logical ranks. The mapping is
// block-structured in domain-index space and intentionally independent from MPI
// so it can be validated in serial before a real communicator is attached.
[[nodiscard]] PicRankTopology3D make_pic_rank_topology_3d(
    std::span<const PicDomain3D> domains,std::size_t ranks_x,std::size_t ranks_y,std::size_t ranks_z);
[[nodiscard]] std::size_t rank_for_pic_domain_3d(const PicRankTopology3D& topology,std::size_t domain_id);

struct PicTransportEnvelope3D {
    PicTransportPayloadKind3D kind{PicTransportPayloadKind3D::particle_migration};
    std::size_t source_domain_id{};
    std::size_t destination_domain_id{};
    std::size_t source_rank{};
    std::size_t destination_rank{};
    ParticleMigrationMessage3D particle_migration;
    ScalarGuardCellMessage3D scalar_guard;
};

struct PicTransportDiagnostics3D {
    std::size_t envelopes{};
    std::size_t particle_messages{};
    std::size_t scalar_guard_messages{};
    std::size_t particle_payload_count{};
    std::size_t scalar_guard_value_count{};
    std::size_t same_rank_messages{};
    std::size_t remote_rank_messages{};
};

// In-memory deterministic transport used to validate the same rank-addressed
// message contract that an MPI layer will later send/receive.
class InMemoryPicTransport3D {
public:
    explicit InMemoryPicTransport3D(std::size_t rank_count);
    void post(std::span<const PicTransportEnvelope3D> envelopes);
    [[nodiscard]] std::vector<PicTransportEnvelope3D> receive_rank(std::size_t rank);
    [[nodiscard]] std::size_t pending_count() const noexcept;
    void clear() noexcept;
private:
    std::size_t rank_count_{};
    std::vector<std::vector<PicTransportEnvelope3D>> inbox_by_rank_;
};

[[nodiscard]] std::vector<PicTransportEnvelope3D> make_pic_transport_envelopes_3d(
    std::span<const ParticleMigrationMessage3D> particle_messages,
    std::span<const ScalarGuardCellMessage3D> scalar_guard_messages,
    const PicRankTopology3D& topology);
[[nodiscard]] PicTransportDiagnostics3D summarize_pic_transport_3d(
    std::span<const PicTransportEnvelope3D> envelopes);

struct PicTransportedExchange3D {
    std::vector<std::vector<PicParticle3D>> particles_by_domain;
    std::vector<ScalarGuardedBlock3D> guarded_blocks;
    PicTransportDiagnostics3D diagnostics;
};

// Executes a complete communicator-shaped exchange through the deterministic
// in-memory transport. It is intentionally serial but validates ordering, rank
// addressing and reconstruction semantics before MPI is introduced.
[[nodiscard]] PicTransportedExchange3D exchange_pic_messages_in_memory_3d(
    std::vector<std::vector<PicParticle3D>> retained_by_domain,
    std::span<const ParticleMigrationMessage3D> particle_messages,
    const std::vector<std::vector<double>>& cell_values_by_domain,
    std::span<const ScalarGuardCellMessage3D> scalar_guard_messages,
    std::span<const PicDomain3D> domains,std::size_t guard_cells,double exterior_value,
    const PicRankTopology3D& topology);

struct PicSerializedEnvelope3D {
    std::vector<std::uint8_t> bytes;
};

struct PicSerializedExchangePlan3D {
    std::size_t rank_count{};
    std::vector<std::vector<PicSerializedEnvelope3D>> envelopes_by_destination_rank;
    std::vector<std::size_t> message_count_by_destination_rank;
    std::vector<std::size_t> byte_count_by_destination_rank;
};

// ABI-local binary serialization for communicator backends. Each envelope is
// self-contained and versioned; the v0.10.5 MPI wrapper can transport these
// opaque records without knowing particle/guard internals. The format is
// intentionally deterministic for regression tests but not a portable file
// format across endian/ABI boundaries.
[[nodiscard]] PicSerializedEnvelope3D serialize_pic_transport_envelope_3d(
    const PicTransportEnvelope3D& envelope);
[[nodiscard]] PicTransportEnvelope3D deserialize_pic_transport_envelope_3d(
    const PicSerializedEnvelope3D& serialized);
[[nodiscard]] std::vector<PicSerializedEnvelope3D> serialize_pic_transport_envelopes_3d(
    std::span<const PicTransportEnvelope3D> envelopes);
[[nodiscard]] std::vector<PicTransportEnvelope3D> deserialize_pic_transport_envelopes_3d(
    std::span<const PicSerializedEnvelope3D> serialized);
[[nodiscard]] PicSerializedExchangePlan3D plan_serialized_pic_exchange_3d(
    std::span<const PicTransportEnvelope3D> envelopes,std::size_t rank_count);
[[nodiscard]] PicTransportDiagnostics3D summarize_pic_serialized_transport_3d(
    std::span<const PicSerializedEnvelope3D> serialized);


struct DistributedPicExchangeRound3D {
    std::vector<std::vector<PicParticle3D>> particles_by_domain;
    std::vector<ScalarGuardedBlock3D> scalar_guarded_blocks;
    PicTransportDiagnostics3D transport_diagnostics;
    PicSerializedExchangePlan3D serialized_plan;
    std::size_t serialized_message_count{};
    std::size_t serialized_byte_count{};
    std::vector<ParticleMigrationRecord3D> outside_particles;
};

// End-to-end rank-local distributed PIC exchange round. This packs particles and
// scalar guard cells, addresses envelopes to logical ranks, serializes per-rank
// payloads, deserializes them again and applies the results. It is intentionally
// deterministic and MPI-free so the complete distributed exchange contract can
// be regression-tested before a real communicator or GPU backend is attached.
[[nodiscard]] DistributedPicExchangeRound3D run_serialized_distributed_pic_exchange_round_3d(
    const std::vector<std::vector<PicParticle3D>>& particles_by_domain,
    const std::vector<std::vector<double>>& scalar_cell_values_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    std::size_t guard_cells,double exterior_value,const PicRankTopology3D& topology);



struct DistributedStaggeredPicStep3DConfig {
    std::size_t guard_cells{1U};
    double exterior_value{};
    double dt_s{1.0e-12};
    bool exchange_em_field_guards{true};
};

struct DistributedStaggeredPicStep3D {
    std::vector<std::vector<PicParticle3D>> particles_by_domain;
    ElectromagneticGuardedFieldBlocks3D guarded_fields;
    DistributedPicExchangeRound3D exchange_round;
    std::size_t particle_count_before{};
    std::size_t particle_count_after_push{};
    std::size_t particle_count_after_migration{};
    double max_particle_displacement_m{};
};

// Rank-local distributed 3-D EM-PIC step bridge. Particles are pushed with
// local domain electromagnetic fields, then passed through the serialized
// particle/guard exchange contract. Field guards are also exchanged through the
// six-component EM guard wrapper. This validates the timestep/exchange ordering
// before the same messages are wired to a real MPI runtime.
[[nodiscard]] DistributedStaggeredPicStep3D run_serialized_distributed_staggered_pic_step_3d(
    const std::vector<std::vector<PicParticle3D>>& particles_by_domain,
    const ElectromagneticFieldBlocks3D& fields_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    const DistributedStaggeredPicStep3DConfig& step_config,
    const PicRankTopology3D& topology);

#if defined(CFD_HAS_MPI)
// Real MPI transport for the same serialized envelope contract. It uses an
// allgather baseline so the interface is validated before adding optimized
// sparse neighbour schedules. Returned envelopes are exactly those addressed to
// this rank.
[[nodiscard]] std::vector<PicSerializedEnvelope3D> exchange_pic_serialized_envelopes_mpi_3d(
    std::span<const PicSerializedEnvelope3D> outbound,MPI_Comm communicator);
#endif

struct ElectrostaticField3D {
    std::vector<double> electric_x_v_per_m;
    std::vector<double> electric_y_v_per_m;
    std::vector<double> electric_z_v_per_m;
};

struct ChargeConservingCurrent3D {
    std::vector<double> charge_density_old_c_per_m3;
    std::vector<double> charge_density_new_c_per_m3;
    std::vector<double> current_x_a_per_m2;
    std::vector<double> current_y_a_per_m2;
    std::vector<double> current_z_a_per_m2;
    double continuity_linf_residual{};
};

struct ElectrostaticPic3DConfig {
    std::size_t nx{16U};
    std::size_t ny{16U};
    std::size_t nz{16U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double length_z_m{1.0};
    double relative_permittivity{1.0};
    double dt_s{1.0e-3};
    bool neutralize_mean_charge{true};
};

[[nodiscard]] std::vector<double> deposit_cic_charge_density_3d(
    std::span<const PicParticle3D> particles,std::size_t nx,std::size_t ny,std::size_t nz,
    double length_x_m,double length_y_m,double length_z_m);

[[nodiscard]] ElectrostaticField3D periodic_electric_field_from_charge_density_3d(
    std::span<const double> charge_density_c_per_m3,std::size_t nx,std::size_t ny,std::size_t nz,
    double length_x_m,double length_y_m,double length_z_m,double relative_permittivity,
    bool remove_mean = true);

[[nodiscard]] ChargeConservingCurrent3D deposit_charge_conserving_current_3d(
    std::span<const PicParticle3D> particles_old,std::span<const PicParticle3D> particles_new,
    std::size_t nx,std::size_t ny,std::size_t nz,
    double length_x_m,double length_y_m,double length_z_m,double dt_s);

// Finite-volume local current reconstruction. Unlike the spectral helper above,
// this uses compact cumulative face fluxes along x, y and z to satisfy the
// periodic backward-difference continuity equation. It is intended as the first
// production-path bridge toward local Villasenor-Buneman/Esirkepov deposition.
[[nodiscard]] ChargeConservingCurrent3D deposit_charge_conserving_current_3d_local(
    std::span<const PicParticle3D> particles_old,std::span<const PicParticle3D> particles_new,
    std::size_t nx,std::size_t ny,std::size_t nz,
    double length_x_m,double length_y_m,double length_z_m,double dt_s);

enum class CurrentDeposition3DMode { spectral_continuity, local_finite_volume };

// Periodic 3-D electrostatic PIC foundation. This deliberately starts with
// spectral Poisson, trilinear CIC deposition/gather and a spectral continuity
// current reconstruction so future 3-D EM-PIC work has validated 3-D particle
// geometry before adding Yee staggered fields and PML/wall boundaries.
class ElectrostaticPic3D {
public:
    explicit ElectrostaticPic3D(ElectrostaticPic3DConfig config = {});
    void set_particles(std::vector<PicParticle3D> particles);
    [[nodiscard]] const std::vector<PicParticle3D>& particles() const noexcept { return particles_; }
    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& electric_x() const noexcept { return electric_x_; }
    [[nodiscard]] const std::vector<double>& electric_y() const noexcept { return electric_y_; }
    [[nodiscard]] const std::vector<double>& electric_z() const noexcept { return electric_z_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void deposit_and_solve();
    [[nodiscard]] Vec3 gather_electric_field(Vec3 position_m) const;
    void step(std::size_t steps = 1U);

private:
    ElectrostaticPic3DConfig config_;
    std::vector<PicParticle3D> particles_;
    std::vector<double> charge_density_;
    std::vector<double> electric_x_;
    std::vector<double> electric_y_;
    std::vector<double> electric_z_;
    double time_s_{};

    [[nodiscard]] double dx() const noexcept { return config_.length_x_m/static_cast<double>(config_.nx); }
    [[nodiscard]] double dy() const noexcept { return config_.length_y_m/static_cast<double>(config_.ny); }
    [[nodiscard]] double dz() const noexcept { return config_.length_z_m/static_cast<double>(config_.nz); }
    [[nodiscard]] Vec3 wrap_position(Vec3 position_m) const;
    [[nodiscard]] double trilinear(std::span<const double> nodal,Vec3 position_m) const;
};


struct ElectromagneticPic3DConfig {
    std::size_t nx{12U};
    std::size_t ny{12U};
    std::size_t nz{12U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double length_z_m{1.0};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    double dt_s{1.0e-12};
    bool solve_longitudinal_poisson{false};
    bool neutralize_mean_charge{true};
};

struct ElectromagneticPic3DDiagnostics {
    double field_energy_j{};
    double particle_kinetic_energy_j{};
    double total_energy_j{};
    double charge_continuity_linf_residual{};
    CollisionStatistics collision_statistics{};
    std::size_t particle_count{};
};


enum class GridBoundaryMode3D { periodic,electric_wall,absorbing_sponge };
struct GridBoundary3DConfig {
    GridBoundaryMode3D mode{GridBoundaryMode3D::periodic};
    std::size_t sponge_cells{2U};
    double sponge_strength{3.0};
    double sponge_polynomial_order{2.0};
};
void apply_electromagnetic_field_boundary_3d(std::vector<double>& electric_x_v_per_m,
                                             std::vector<double>& electric_y_v_per_m,
                                             std::vector<double>& electric_z_v_per_m,
                                             std::vector<double>& magnetic_x_t,
                                             std::vector<double>& magnetic_y_t,
                                             std::vector<double>& magnetic_z_t,
                                             std::size_t nx,std::size_t ny,std::size_t nz,
                                             const GridBoundary3DConfig& boundary);

struct StaggeredElectromagneticPic3DConfig {
    std::size_t nx{10U};
    std::size_t ny{10U};
    std::size_t nz{10U};
    double length_x_m{1.0};
    double length_y_m{1.0};
    double length_z_m{1.0};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    double dt_s{1.0e-12};
    GridBoundary3DConfig field_boundary{};
    ParticleWallMode particle_boundary{ParticleWallMode::specular_reflect};
    bool periodic_particles{true};
    bool remove_absorbed_particles{true};
    bool enable_secondary_yield_report{false};
    SecondaryEmissionModel secondary_model{};
    CurrentDeposition3DMode current_deposition{CurrentDeposition3DMode::spectral_continuity};
    bool sort_particles_by_cell{false};
    std::size_t particle_sort_interval{1U};
};

struct StaggeredElectromagneticPic3DDiagnostics {
    double field_energy_j{};
    double particle_kinetic_energy_j{};
    double total_energy_j{};
    double charge_continuity_linf_residual{};
    std::size_t particle_count{};
    std::size_t absorbed_particles{};
    double secondary_macro_weight{};
    CollisionStatistics collision_statistics{};
    std::size_t particle_sort_passes{};
    std::size_t occupied_particle_cells{};
};

// True 3-D/3-V Yee-style electromagnetic PIC reference. Component arrays are
// stored with explicit Yee semantics: Ex/Ey/Ez live on oriented electric edges
// and Bx/By/Bz on staggered magnetic faces. For this compact CPU baseline the
// arrays keep equal extents, while gather/update methods apply the component
// offsets. Current coupling reuses the validated 3-D spectral continuity
// reconstruction; local Esirkepov/Villasenor-Buneman deposition is a future
// production upgrade.
class StaggeredElectromagneticPic3D {
public:
    explicit StaggeredElectromagneticPic3D(StaggeredElectromagneticPic3DConfig config = {});

    void set_particles(std::vector<PicParticle3D> particles);
    [[nodiscard]] const std::vector<PicParticle3D>& particles() const noexcept { return particles_; }

    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& current_x() const noexcept { return current_x_; }
    [[nodiscard]] const std::vector<double>& current_y() const noexcept { return current_y_; }
    [[nodiscard]] const std::vector<double>& current_z() const noexcept { return current_z_; }
    [[nodiscard]] const std::vector<double>& electric_x() const noexcept { return electric_x_; }
    [[nodiscard]] const std::vector<double>& electric_y() const noexcept { return electric_y_; }
    [[nodiscard]] const std::vector<double>& electric_z() const noexcept { return electric_z_; }
    [[nodiscard]] const std::vector<double>& magnetic_x() const noexcept { return magnetic_x_; }
    [[nodiscard]] const std::vector<double>& magnetic_y() const noexcept { return magnetic_y_; }
    [[nodiscard]] const std::vector<double>& magnetic_z() const noexcept { return magnetic_z_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void set_fields(std::span<const double> electric_x_v_per_m,
                    std::span<const double> electric_y_v_per_m,
                    std::span<const double> electric_z_v_per_m,
                    std::span<const double> magnetic_x_t,
                    std::span<const double> magnetic_y_t,
                    std::span<const double> magnetic_z_t);
    void initialize_z_polarized_mode(double electric_z_amplitude_v_per_m,
                                     std::size_t mode_x = 1U,
                                     std::size_t mode_y = 0U,
                                     std::size_t mode_z = 0U);
    [[nodiscard]] ElectromagneticField gather_field(Vec3 position_m) const;
    [[nodiscard]] StaggeredElectromagneticPic3DDiagnostics diagnostics() const;
    void step(std::size_t steps = 1U,const NeutralCollisionModel* collision_model = nullptr);

private:
    StaggeredElectromagneticPic3DConfig config_;
    std::vector<PicParticle3D> particles_;
    std::vector<double> charge_density_;
    std::vector<double> current_x_;
    std::vector<double> current_y_;
    std::vector<double> current_z_;
    std::vector<double> electric_x_;
    std::vector<double> electric_y_;
    std::vector<double> electric_z_;
    std::vector<double> magnetic_x_;
    std::vector<double> magnetic_y_;
    std::vector<double> magnetic_z_;
    double time_s_{};
    double last_continuity_residual_{};
    std::size_t absorbed_particles_{};
    double secondary_macro_weight_{};
    CollisionStatistics accumulated_collisions_{};
    std::size_t particle_sort_passes_{};
    std::size_t occupied_particle_cells_{};
    std::size_t step_counter_{};

    [[nodiscard]] double dx() const noexcept { return config_.length_x_m/static_cast<double>(config_.nx); }
    [[nodiscard]] double dy() const noexcept { return config_.length_y_m/static_cast<double>(config_.ny); }
    [[nodiscard]] double dz() const noexcept { return config_.length_z_m/static_cast<double>(config_.nz); }
    [[nodiscard]] double epsilon() const noexcept;
    [[nodiscard]] double mu() const noexcept;
    [[nodiscard]] double medium_light_speed() const noexcept;
    [[nodiscard]] bool periodic_fields() const noexcept;
    [[nodiscard]] Vec3 place_particle(Vec3 position_m) const;
    [[nodiscard]] double sample(std::span<const double> values,Vec3 position_m,double offset_x,double offset_y,double offset_z) const;
    [[nodiscard]] double diff_x_forward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    [[nodiscard]] double diff_y_forward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    [[nodiscard]] double diff_z_forward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    [[nodiscard]] double diff_x_backward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    [[nodiscard]] double diff_y_backward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    [[nodiscard]] double diff_z_backward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    void advance_magnetic(double dt_s);
    void update_electric(double dt_s);
    void apply_field_boundary();
    void apply_particle_boundary();
    void maybe_sort_particles();
};

// Periodic 3-D/3-V electromagnetic PIC baseline. This is a compact nodal-grid
// Maxwell/PIC reference: E/B are stored on the same periodic grid, curls use
// centered periodic differences, particles gather trilinear E/B and advance with
// the Vay pusher, and Jx/Jy/Jz are reconstructed from old/new CIC charge
// densities so the spectral continuity equation is satisfied by construction.
class ElectromagneticPic3D {
public:
    explicit ElectromagneticPic3D(ElectromagneticPic3DConfig config = {});

    void set_particles(std::vector<PicParticle3D> particles);
    [[nodiscard]] const std::vector<PicParticle3D>& particles() const noexcept { return particles_; }

    [[nodiscard]] const std::vector<double>& charge_density() const noexcept { return charge_density_; }
    [[nodiscard]] const std::vector<double>& current_x() const noexcept { return current_x_; }
    [[nodiscard]] const std::vector<double>& current_y() const noexcept { return current_y_; }
    [[nodiscard]] const std::vector<double>& current_z() const noexcept { return current_z_; }
    [[nodiscard]] const std::vector<double>& electric_x() const noexcept { return electric_x_; }
    [[nodiscard]] const std::vector<double>& electric_y() const noexcept { return electric_y_; }
    [[nodiscard]] const std::vector<double>& electric_z() const noexcept { return electric_z_; }
    [[nodiscard]] const std::vector<double>& magnetic_x() const noexcept { return magnetic_x_; }
    [[nodiscard]] const std::vector<double>& magnetic_y() const noexcept { return magnetic_y_; }
    [[nodiscard]] const std::vector<double>& magnetic_z() const noexcept { return magnetic_z_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    void set_fields(std::span<const double> electric_x_v_per_m,
                    std::span<const double> electric_y_v_per_m,
                    std::span<const double> electric_z_v_per_m,
                    std::span<const double> magnetic_x_t,
                    std::span<const double> magnetic_y_t,
                    std::span<const double> magnetic_z_t);
    void initialize_z_polarized_mode(double electric_z_amplitude_v_per_m,
                                     std::size_t mode_x = 1U,
                                     std::size_t mode_y = 0U,
                                     std::size_t mode_z = 0U);
    void deposit_sources();
    [[nodiscard]] ElectromagneticField gather_field(Vec3 position_m) const;
    [[nodiscard]] ElectromagneticPic3DDiagnostics diagnostics() const;
    void step(std::size_t steps = 1U,const NeutralCollisionModel* collision_model = nullptr);

private:
    ElectromagneticPic3DConfig config_;
    std::vector<PicParticle3D> particles_;
    std::vector<double> charge_density_;
    std::vector<double> current_x_;
    std::vector<double> current_y_;
    std::vector<double> current_z_;
    std::vector<double> electric_x_;
    std::vector<double> electric_y_;
    std::vector<double> electric_z_;
    std::vector<double> magnetic_x_;
    std::vector<double> magnetic_y_;
    std::vector<double> magnetic_z_;
    double time_s_{};
    double last_continuity_residual_{};
    CollisionStatistics accumulated_collisions_{};

    [[nodiscard]] double dx() const noexcept { return config_.length_x_m/static_cast<double>(config_.nx); }
    [[nodiscard]] double dy() const noexcept { return config_.length_y_m/static_cast<double>(config_.ny); }
    [[nodiscard]] double dz() const noexcept { return config_.length_z_m/static_cast<double>(config_.nz); }
    [[nodiscard]] double epsilon() const noexcept;
    [[nodiscard]] double mu() const noexcept;
    [[nodiscard]] double medium_light_speed() const noexcept;
    [[nodiscard]] Vec3 wrap_position(Vec3 position_m) const;
    [[nodiscard]] double trilinear(std::span<const double> nodal,Vec3 position_m) const;
    [[nodiscard]] double ddx(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    [[nodiscard]] double ddy(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    [[nodiscard]] double ddz(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const;
    void advance_magnetic(double dt_s);
    void update_electric(double dt_s);
};

} // namespace cfd::particle
