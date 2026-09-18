#pragma once
#include "cfd/battery/pack.hpp"
#include <array>
namespace cfd::battery {
struct ThermalVoxel {
    double volumetric_heat_capacity_j_per_m3_k{2.0e6};
    std::array<double,3> conductivity_w_per_m_k{1.0,1.0,1.0};
};
struct ThermalBoundary {
    double convection_w_per_m2_k{0.0};
    double ambient_temperature_k{298.15};
};
struct ThermalGridConfig {
    std::array<std::size_t,3> cells{8,8,8};
    std::array<double,3> spacing_m{0.005,0.005,0.005};
    // x-, x+, y-, y+, z-, z+. Zero convection means insulated.
    std::array<ThermalBoundary,6> boundaries{};
    double initial_temperature_k{298.15};
    std::size_t maximum_iterations{2000};
    double relative_tolerance{1e-12};
};
struct ThermalStepResult {
    double supplied_power_w{}, outward_cooling_power_w{}, stored_energy_change_j{};
    double energy_balance_error_j{}, maximum_temperature_k{};
    std::size_t iterations{};
};
// Cell-centred finite volumes, x fastest. Backward Euler, harmonic face
// conductance, diagonal anisotropy. One material broadcasts to all voxels.
class ThermalGrid3D {
public:
    explicit ThermalGrid3D(ThermalGridConfig config = {}, std::vector<ThermalVoxel> materials = {ThermalVoxel{}});
    void reset(std::span<const double> temperature_k);
    [[nodiscard]] ThermalStepResult step(double dt_s, std::span<const double> voxel_power_w);
    [[nodiscard]] const std::vector<double>& temperature() const noexcept { return temperature_; }
    [[nodiscard]] const std::vector<double>& heat_capacity_j_per_k() const noexcept { return capacity_; }
private:
    struct Face { std::size_t a,b; double conductance; };
    std::vector<Face> faces_;
    ThermalGridConfig config_;
    std::vector<double> temperature_, capacity_, cooling_, ambient_source_;
};
struct ElectrothermalPackResult {
    PackStepResult electrical;
    ThermalStepResult thermal;
};
// Explicit partitioned electrothermal coupling. Each cell maps to disjoint
// voxel indices; unassigned voxels may represent passive pack material.
class ElectrothermalBatteryPack {
public:
    ElectrothermalBatteryPack(std::span<const LithiumIonCellConfig> cells, PackControl control,
                              ThermalGrid3D thermal, std::vector<std::vector<std::size_t>> cell_voxels);
    [[nodiscard]] ElectrothermalPackResult step(double current_a, double dt_s);
    [[nodiscard]] const SeriesBatteryPack& electrical() const noexcept { return pack_; }
    [[nodiscard]] const ThermalGrid3D& thermal() const noexcept { return thermal_; }
private:
    [[nodiscard]] std::vector<double> cell_temperatures(const ThermalGrid3D& field) const;
    SeriesBatteryPack pack_;
    ThermalGrid3D thermal_;
    std::vector<std::vector<std::size_t>> regions_;
    double maximum_temperature_k_;
};
}
