#pragma once

#include <span>
#include <string>
#include <vector>
#include <filesystem>

namespace cfd::chemistry {

struct ThermoSpecies {
    std::string name;
    double standard_chemical_potential{}; // J/mol at the phase reference state
};

class PhaseThermo {
public:
    virtual ~PhaseThermo() = default;
    [[nodiscard]] virtual std::size_t species_count() const noexcept = 0;
    [[nodiscard]] virtual std::vector<double> activities(std::span<const double> mole_fractions,
                                                          double temperature,
                                                          double pressure) const = 0;
    [[nodiscard]] virtual std::vector<double> chemical_potentials(std::span<const double> mole_fractions,
                                                                   double temperature,
                                                                   double pressure) const = 0;
};

class IdealGasPhase final : public PhaseThermo {
public:
    explicit IdealGasPhase(std::vector<ThermoSpecies> species,double reference_pressure=101325.0);
    [[nodiscard]] std::size_t species_count() const noexcept override { return species_.size(); }
    [[nodiscard]] std::vector<double> activities(std::span<const double> mole_fractions,
                                                  double temperature,double pressure) const override;
    [[nodiscard]] std::vector<double> chemical_potentials(std::span<const double> mole_fractions,
                                                           double temperature,double pressure) const override;
private:
    std::vector<ThermoSpecies> species_;
    double reference_pressure_{};
};

class IdealSolutionPhase final : public PhaseThermo {
public:
    explicit IdealSolutionPhase(std::vector<ThermoSpecies> species);
    [[nodiscard]] std::size_t species_count() const noexcept override { return species_.size(); }
    [[nodiscard]] std::vector<double> activities(std::span<const double> mole_fractions,
                                                  double temperature,double pressure) const override;
    [[nodiscard]] std::vector<double> chemical_potentials(std::span<const double> mole_fractions,
                                                           double temperature,double pressure) const override;
private:
    std::vector<ThermoSpecies> species_;
};

[[nodiscard]] std::vector<ThermoSpecies> load_thermo_csv(const std::filesystem::path& path);

} // namespace cfd::chemistry
