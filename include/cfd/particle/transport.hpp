#pragma once

#include "cfd/particle/electromagnetic.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace cfd::particle {

enum class TrackStatus { alive, stopped, escaped, absorbed };
enum class StepLimiter { user, geometry_boundary, continuous_process, discrete_process, energy_cut };
enum class TransportProcessKind { continuous_energy_loss, discrete_interaction, absorber };

struct TransportParticleDefinition {
    std::string name{"particle"};
    double charge_c{};
    double mass_kg{1.0};
};

struct TransportMaterial {
    std::string name{"material"};
    double density_kg_per_m3{1.0};
    double production_cut_energy_ev{0.0};
};

struct TransportRegion {
    std::string name{"region"};
    AxisAlignedParticleBox box{};
    std::size_t material_index{};
    std::size_t sensitive_detector_index{static_cast<std::size_t>(-1)};
};

struct TransportProcess {
    std::string name{"process"};
    TransportProcessKind kind{TransportProcessKind::continuous_energy_loss};
    std::size_t particle_index{static_cast<std::size_t>(-1)}; // all particles when out of range
    std::size_t material_index{static_cast<std::size_t>(-1)}; // all materials when out of range
    double stopping_power_ev_per_m{};                         // continuous dE/dx
    double physical_interaction_length_m{};                   // discrete GPIL analogue
    double energy_loss_fraction{};                            // discrete parent energy loss
    double secondary_energy_fraction{};                       // emitted secondary energy share
    double absorption_threshold_ev{};                         // absorber threshold
};

struct TransportTrack {
    Vec3 position_m{};
    Vec3 direction{}; // normalized on entry to the stepper
    double kinetic_energy_ev{};
    double time_s{};
    double track_length_m{};
    double weight{1.0};
    std::size_t particle_index{};
    TrackStatus status{TrackStatus::alive};
};

struct TransportSecondary {
    TransportTrack track{};
    std::string creator_process;
};

struct TransportStep {
    Vec3 pre_position_m{};
    Vec3 post_position_m{};
    double step_length_m{};
    double pre_energy_ev{};
    double post_energy_ev{};
    double energy_deposit_ev{};
    double global_time_s{};
    std::size_t material_index{static_cast<std::size_t>(-1)};
    std::size_t region_index{static_cast<std::size_t>(-1)};
    StepLimiter limiter{StepLimiter::user};
    std::size_t limiting_process_index{static_cast<std::size_t>(-1)};
    std::vector<TransportSecondary> secondaries;
};

struct TransportScoring {
    double total_energy_deposit_ev{};
    double total_track_length_m{};
    std::size_t steps{};
    std::size_t secondaries{};
    std::size_t boundary_crossings{};
};

struct TransportSensitiveDetector {
    std::string name{"detector"};
    double minimum_energy_deposit_ev{};
    bool record_zero_deposit_steps{false};
};

struct TransportWorld {
    std::vector<TransportParticleDefinition> particles;
    std::vector<TransportMaterial> materials;
    std::vector<TransportRegion> regions;
    std::vector<TransportProcess> processes;
    std::vector<TransportSensitiveDetector> sensitive_detectors;
};

struct TransportPhysicsList {
    std::string name{"physics-list"};
    std::vector<TransportProcess> processes;
    double production_cut_energy_ev{-1.0};
};

struct TransportRandom {
    std::uint64_t state{0x9e3779b97f4a7c15ULL};
    [[nodiscard]] double uniform_open01();
};

struct TransportBvhNode {
    AxisAlignedParticleBox box{};
    std::size_t left{static_cast<std::size_t>(-1)};
    std::size_t right{static_cast<std::size_t>(-1)};
    std::size_t first{};
    std::size_t count{};
};

struct TransportRegionBvh {
    std::vector<TransportBvhNode> nodes;
    std::vector<std::size_t> region_indices;
};

struct TransportHit {
    std::size_t detector_index{};
    std::size_t region_index{};
    std::size_t material_index{};
    std::size_t particle_index{};
    Vec3 position_m{};
    double energy_deposit_ev{};
    double step_length_m{};
    double global_time_s{};
    double weight{1.0};
    StepLimiter limiter{StepLimiter::user};
    std::string process_name;
};

struct TransportHitCollection {
    std::vector<TransportHit> hits;
    double total_energy_deposit_ev{};
    double total_weighted_energy_deposit_ev{};
};

struct TransportDoseGrid3D {
    std::size_t nx{},ny{},nz{};
    AxisAlignedParticleBox box{};
    double density_kg_per_m3{1000.0};
    std::vector<double> energy_deposit_ev;
    std::vector<double> dose_gy;
};

struct TransportConfig {
    double max_step_m{1.0};
    double energy_cut_ev{0.0};
    double boundary_epsilon_m{1.0e-12};
    double light_speed_m_per_s{299792458.0};
};

struct TransportRunResult {
    TransportTrack primary{};
    std::vector<TransportStep> steps;
    std::vector<TransportSecondary> secondaries;
    TransportScoring scoring{};
};

[[nodiscard]] Vec3 normalize_direction(Vec3 direction);
[[nodiscard]] bool contains(const AxisAlignedParticleBox& box,Vec3 position,double tolerance_m=0.0);
[[nodiscard]] std::size_t locate_region(const TransportWorld& world,Vec3 position,double tolerance_m=0.0);
[[nodiscard]] double distance_to_box_boundary(const AxisAlignedParticleBox& box,Vec3 position,Vec3 direction);
[[nodiscard]] bool process_applies(const TransportProcess& process,const TransportTrack& track,std::size_t material_index);
[[nodiscard]] TransportWorld apply_physics_list(TransportWorld world,const TransportPhysicsList& physics_list);
[[nodiscard]] double sample_exponential_interaction_length_m(double mean_free_path_m,TransportRandom& rng);
[[nodiscard]] TransportRegionBvh build_region_bvh(const TransportWorld& world,std::size_t leaf_size=4U);
[[nodiscard]] std::size_t locate_region_bvh(const TransportWorld& world,const TransportRegionBvh& bvh,Vec3 position,double tolerance_m=0.0);
[[nodiscard]] TransportHitCollection collect_transport_hits(const TransportWorld& world,std::span<const TransportStep> steps);
[[nodiscard]] TransportDoseGrid3D make_transport_dose_grid(std::size_t nx,std::size_t ny,std::size_t nz,AxisAlignedParticleBox box,double density_kg_per_m3);
void score_hits_to_dose_grid(const TransportHitCollection& hits,TransportDoseGrid3D& grid);
[[nodiscard]] std::vector<double> dose_rate_to_sar_w_per_kg(const TransportDoseGrid3D& grid,double exposure_time_s);
[[nodiscard]] std::vector<double> project_dose_grid_to_pennes2d_sar(const TransportDoseGrid3D& grid,double exposure_time_s,std::size_t nx,std::size_t ny);
[[nodiscard]] TransportStep transport_one_step(const TransportWorld& world,TransportTrack& track,const TransportConfig& config);
[[nodiscard]] TransportRunResult transport_track(const TransportWorld& world,TransportTrack primary,const TransportConfig& config,std::size_t max_steps=1024U);
[[nodiscard]] TransportStep transport_one_step_stochastic(const TransportWorld& world,TransportTrack& track,const TransportConfig& config,TransportRandom& rng);
[[nodiscard]] TransportRunResult transport_track_stochastic(const TransportWorld& world,TransportTrack primary,const TransportConfig& config,TransportRandom& rng,std::size_t max_steps=1024U);

} // namespace cfd::particle
