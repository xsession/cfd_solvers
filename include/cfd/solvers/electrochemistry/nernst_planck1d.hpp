#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

namespace cfd::electrochemistry {

enum class TransportBoundary1D {
    periodic,
    no_flux
};

struct NernstPlanck1DConfig {
    std::size_t cells{128};
    double length{1.0};
    double temperature{298.15};
    double dt{1.0e-4};
    double advection_velocity{0.0};
    TransportBoundary1D boundary{TransportBoundary1D::periodic};
};

struct NernstPlanckSpecies1D {
    std::string name;
    int charge{};
    double diffusivity{};
    std::vector<double> concentration;
};

class NernstPlanck1D {
public:
    explicit NernstPlanck1D(NernstPlanck1DConfig config);

    [[nodiscard]] const NernstPlanck1DConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::vector<double>& potential() const noexcept { return potential_; }
    [[nodiscard]] const std::vector<NernstPlanckSpecies1D>& species() const noexcept { return species_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] double dx() const noexcept { return dx_; }

    void set_potential(const std::function<double(double)>& fn);
    void set_potential(std::span<const double> values);
    std::size_t add_species(std::string name,
                            int charge,
                            double diffusivity,
                            const std::function<double(double)>& initial_concentration);
    std::size_t add_species(std::string name,
                            int charge,
                            double diffusivity,
                            double uniform_concentration);

    void step(std::size_t count = 1U);
    [[nodiscard]] double total_amount(std::size_t species_index) const;
    [[nodiscard]] double minimum_concentration(std::size_t species_index) const;

private:
    NernstPlanck1DConfig config_;
    double dx_{};
    double time_{};
    std::vector<double> potential_;
    std::vector<NernstPlanckSpecies1D> species_;
    std::vector<double> flux_;
    std::vector<double> next_;

    void step_species(NernstPlanckSpecies1D& species);
};

} // namespace cfd::electrochemistry
