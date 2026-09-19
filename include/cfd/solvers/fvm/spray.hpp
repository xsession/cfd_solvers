#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace cfd::fvm {

struct SprayParcel {
    Vec3 position{};
    Vec3 velocity{};
    double diameter{1.0e-4};
    double mass{1.0e-9};
    double temperature{300.0};
};

struct SprayInjection {
    Vec3 position{};
    Vec3 velocity{};
    double diameter{1.0e-4};
    double temperature{300.0};
    double mass_flow_rate{1.0e-6};
    double start_time{};
    double duration{1.0e-2};
};

// Coupled parcel spray baseline. Parcels receive Stokes drag and D^2-law
// evaporation; the carrier receives vapor-mass, drag-reaction momentum and
// latent-heat source fields on the same PolyMesh.
struct SprayConfig {
    double dt{1.0e-4};
    double liquid_density{1000.0};
    double carrier_dynamic_viscosity{1.8e-5};
    double evaporation_constant{1.0e-8};
    double latent_heat{2.2e6};
    double thermal_relaxation_time{1.0e-3};
    double minimum_diameter{1.0e-8};
    std::size_t maximum_parcels{100000U};
};

class SprayInjectionEvaporation {
public:
    explicit SprayInjectionEvaporation(PolyMesh mesh, SprayConfig config = {});

    void set_injection(SprayInjection injection);
    void clear_injection() noexcept { injection_.reset(); }
    void add_parcel(SprayParcel parcel);
    void clear_parcels() noexcept { parcels_.clear(); }

    void step(std::span<const Vec3> carrier_velocity, std::span<const double> carrier_temperature);
    void step(std::span<const Vec3> carrier_velocity, double carrier_temperature);

    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const SprayConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::vector<SprayParcel>& parcels() const noexcept { return parcels_; }
    [[nodiscard]] const std::vector<double>& vapor_mass_source() const noexcept { return vapor_mass_source_; }
    [[nodiscard]] const std::vector<Vec3>& carrier_momentum_source() const noexcept { return carrier_momentum_source_; }
    [[nodiscard]] const std::vector<double>& carrier_energy_source() const noexcept { return carrier_energy_source_; }
    [[nodiscard]] double total_liquid_mass() const noexcept;
    [[nodiscard]] double injected_mass() const noexcept { return injected_mass_; }
    [[nodiscard]] double evaporated_mass() const noexcept { return evaporated_mass_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }

private:
    PolyMesh mesh_;
    SprayConfig config_;
    std::optional<SprayInjection> injection_;
    std::vector<SprayParcel> parcels_;
    std::vector<double> vapor_mass_source_, carrier_energy_source_;
    std::vector<Vec3> carrier_momentum_source_;
    double injected_mass_{};
    double evaporated_mass_{};
    double time_{};
    std::size_t steps_{};

    void validate_config() const;
    static void validate_vec(Vec3 value, const char* message);
    void validate_parcel(const SprayParcel& parcel) const;
    void validate_injection(const SprayInjection& injection) const;
};

} // namespace cfd::fvm
