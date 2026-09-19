#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::fvm {

// Isothermal/barotropic compressible VOF baseline. The phase fraction and
// mixture mass are transported conservatively with a prescribed volumetric
// face flux; pressure is recovered from the mixture barotropic EOS.
struct CompressibleVofConfig {
    double dt{1.0e-4};
    double liquid_density{1000.0};
    double gas_density{1.2};
    double liquid_bulk_modulus{2.2e9};
    double gas_gamma{1.4};
    double reference_pressure{1.0e5};
    double minimum_pressure{1.0};
    double minimum_density{1.0e-12};
};

class CompressibleVofTransport {
public:
    explicit CompressibleVofTransport(PolyMesh mesh, CompressibleVofConfig config = {});

    void initialize(double liquid_fraction, double pressure);
    void initialize(std::span<const double> liquid_fraction, std::span<const double> pressure);

    template <class F> void initialize(F&& fraction_function, double pressure = 0.0) {
        const double initial_pressure = pressure == 0.0 ? config_.reference_pressure : pressure;
        std::vector<double> fraction(mesh_.cell_count());
        for (std::size_t cell = 0; cell < fraction.size(); ++cell)
            fraction[cell] = static_cast<double>(fraction_function(mesh_.cells()[cell].center));
        std::vector<double> pressures(mesh_.cell_count(), initial_pressure);
        initialize(std::span<const double>(fraction), std::span<const double>(pressures));
    }

    void set_face_flux(std::vector<double> volumetric_flux);
    void step();

    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const CompressibleVofConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::vector<double>& liquid_fraction() const noexcept { return liquid_fraction_; }
    [[nodiscard]] const std::vector<double>& density() const noexcept { return density_; }
    [[nodiscard]] const std::vector<double>& pressure() const noexcept { return pressure_; }
    [[nodiscard]] const std::vector<double>& liquid_density() const noexcept { return liquid_density_; }
    [[nodiscard]] const std::vector<double>& gas_density() const noexcept { return gas_density_; }
    [[nodiscard]] const std::vector<double>& face_flux() const noexcept { return face_flux_; }
    [[nodiscard]] double liquid_volume() const;
    [[nodiscard]] double total_mass() const;
    [[nodiscard]] double minimum_pressure() const;
    [[nodiscard]] double maximum_pressure() const;
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }

private:
    PolyMesh mesh_;
    CompressibleVofConfig config_;
    std::vector<double> liquid_fraction_, mass_, density_, pressure_;
    std::vector<double> liquid_density_, gas_density_, face_flux_;
    double time_{};
    std::size_t steps_{};

    void validate_config() const;
    [[nodiscard]] double mixture_density(double liquid_fraction, double pressure) const;
    [[nodiscard]] double pressure_for_density(double liquid_fraction, double density) const;
    void update_properties();
};

} // namespace cfd::fvm
