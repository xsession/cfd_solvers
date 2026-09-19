#pragma once

#include "cfd/chemistry/kinetics.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <vector>
#include <utility>

namespace cfd::chemistry {

struct Nasa7Polynomial {
    std::array<double, 7> a{};
    double minimum_temperature{200.0};
    double maximum_temperature{6000.0};
    [[nodiscard]] double cp_molar(double temperature) const;
    [[nodiscard]] double h_molar(double temperature) const;
    [[nodiscard]] double s_molar(double temperature) const;
};

struct ShomatePolynomial {
    // NIST form, t=T/1000. cp in J/(mol K), h in kJ/mol, s in J/(mol K).
    double A{}, B{}, C{}, D{}, E{}, F{}, G{}, H{};
    double minimum_temperature{200.0}, maximum_temperature{6000.0};
    [[nodiscard]] double cp_molar(double temperature) const;
    [[nodiscard]] double h_molar_kj(double temperature) const;
    [[nodiscard]] double s_molar(double temperature) const;
};
struct BinaryPitzer11 {
    double beta0{}, beta1{}, c_phi{};
    [[nodiscard]] double mean_activity_coefficient(double molality, double temperature = 298.15) const;
};

class ActivityModel {
public:
    virtual ~ActivityModel() = default;
    [[nodiscard]] virtual std::vector<double>
    activity_coefficients(std::span<const int> charge, std::span<const double> molality, double temperature) const = 0;
};

struct PitzerInteraction {
    std::size_t first{};
    std::size_t second{};
    double beta0{};
    double beta1{};
    double c_phi{};
};

class MulticomponentPitzer final : public ActivityModel {
public:
    explicit MulticomponentPitzer(std::vector<PitzerInteraction> interactions, double a_phi = 0.392, double b = 1.2)
        : interactions_(std::move(interactions)), a_phi_(a_phi), b_(b) {}
    [[nodiscard]] std::vector<double> activity_coefficients(std::span<const int> charge,
                                                            std::span<const double> molality,
                                                            double temperature) const override;

private:
    std::vector<PitzerInteraction> interactions_;
    double a_phi_{};
    double b_{};
};

class DaviesActivity final : public ActivityModel {
public:
    explicit DaviesActivity(double a_parameter = 0.509) : a_(a_parameter) {}
    [[nodiscard]] std::vector<double> activity_coefficients(std::span<const int> charge,
                                                            std::span<const double> molality,
                                                            double temperature) const override;

private:
    double a_{};
};

struct SitPair {
    std::size_t first{};
    std::size_t second{};
    double epsilon{}; // kg/mol
};

class SitActivity final : public ActivityModel {
public:
    explicit SitActivity(std::vector<SitPair> interactions, double a_parameter = 0.509)
        : interactions_(std::move(interactions)), a_(a_parameter) {}
    [[nodiscard]] std::vector<double> activity_coefficients(std::span<const int> charge,
                                                            std::span<const double> molality,
                                                            double temperature) const override;

private:
    std::vector<SitPair> interactions_;
    double a_{};
};

struct BinaryIdealEquilibrium {
    double mole_fraction_a{};
    double mole_fraction_b{};
    double gibbs_molar{};
};

[[nodiscard]] BinaryIdealEquilibrium minimize_binary_ideal_gibbs(double standard_mu_a, double standard_mu_b,
                                                                 double temperature);

[[nodiscard]] double langmuir_surface_coverage(double equilibrium_constant, double activity);
[[nodiscard]] double binary_ion_exchange_fraction(double selectivity, double activity_a, double activity_b);

} // namespace cfd::chemistry
