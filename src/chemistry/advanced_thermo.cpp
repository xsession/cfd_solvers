#include "cfd/chemistry/advanced_thermo.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::chemistry {
namespace {
void validate_temperature(double t, double low, double high) {
    if (!(t > 0.0) || !std::isfinite(t) || t < low || t > high) {
        throw std::out_of_range("temperature outside NASA polynomial range");
    }
}
void validate_solution(std::span<const int> charge, std::span<const double> molality) {
    if (charge.size() != molality.size())
        throw std::invalid_argument("activity-model size mismatch");
    for (double m : molality)
        if (!(m >= 0.0) || !std::isfinite(m))
            throw std::invalid_argument("invalid molality");
}
double ionic_strength(std::span<const int> charge, std::span<const double> molality) {
    double i = 0.0;
    for (std::size_t k = 0; k < molality.size(); ++k) {
        const double z = static_cast<double>(charge[k]);
        i += 0.5 * molality[k] * z * z;
    }
    return i;
}
} // namespace

double Nasa7Polynomial::cp_molar(double t) const {
    validate_temperature(t, minimum_temperature, maximum_temperature);
    const double cp_r = a[0] + a[1] * t + a[2] * t * t + a[3] * t * t * t + a[4] * t * t * t * t;
    return gas_constant * cp_r;
}
double Nasa7Polynomial::h_molar(double t) const {
    validate_temperature(t, minimum_temperature, maximum_temperature);
    const double h_rt =
        a[0] + a[1] * t / 2.0 + a[2] * t * t / 3.0 + a[3] * t * t * t / 4.0 + a[4] * t * t * t * t / 5.0 + a[5] / t;
    return gas_constant * t * h_rt;
}
double Nasa7Polynomial::s_molar(double t) const {
    validate_temperature(t, minimum_temperature, maximum_temperature);
    const double s_r =
        a[0] * std::log(t) + a[1] * t + a[2] * t * t / 2.0 + a[3] * t * t * t / 3.0 + a[4] * t * t * t * t / 4.0 + a[6];
    return gas_constant * s_r;
}

std::vector<double> DaviesActivity::activity_coefficients(std::span<const int> charge, std::span<const double> molality,
                                                          double temperature) const {
    validate_solution(charge, molality);
    if (!(temperature > 0.0) || !std::isfinite(temperature) || !(a_ > 0.0))
        throw std::invalid_argument("invalid Davies controls");
    const double i = ionic_strength(charge, molality);
    const double root = std::sqrt(i);
    std::vector<double> gamma(molality.size(), 1.0);
    for (std::size_t k = 0; k < gamma.size(); ++k) {
        const double z = static_cast<double>(charge[k]);
        const double log10_gamma = -a_ * z * z * (root / (1.0 + root) - 0.3 * i);
        gamma[k] = std::pow(10.0, log10_gamma);
    }
    return gamma;
}

std::vector<double> SitActivity::activity_coefficients(std::span<const int> charge, std::span<const double> molality,
                                                       double temperature) const {
    validate_solution(charge, molality);
    if (!(temperature > 0.0) || !std::isfinite(temperature) || !(a_ > 0.0))
        throw std::invalid_argument("invalid SIT controls");
    const double i = ionic_strength(charge, molality);
    const double root = std::sqrt(i);
    std::vector<double> log10_gamma(molality.size(), 0.0);
    for (std::size_t k = 0; k < molality.size(); ++k) {
        const double z = static_cast<double>(charge[k]);
        log10_gamma[k] = -a_ * z * z * root / (1.0 + 1.5 * root);
    }
    for (const auto& pair : interactions_) {
        if (pair.first >= molality.size() || pair.second >= molality.size() || !std::isfinite(pair.epsilon)) {
            throw std::invalid_argument("invalid SIT interaction");
        }
        log10_gamma[pair.first] += pair.epsilon * molality[pair.second];
        if (pair.first != pair.second)
            log10_gamma[pair.second] += pair.epsilon * molality[pair.first];
    }
    std::vector<double> gamma(log10_gamma.size());
    for (std::size_t k = 0; k < gamma.size(); ++k)
        gamma[k] = std::pow(10.0, log10_gamma[k]);
    return gamma;
}

BinaryIdealEquilibrium minimize_binary_ideal_gibbs(double mu_a, double mu_b, double t) {
    if (!(t > 0.0) || !std::isfinite(t) || !std::isfinite(mu_a) || !std::isfinite(mu_b)) {
        throw std::invalid_argument("invalid binary Gibbs equilibrium input");
    }
    const double delta = (mu_b - mu_a) / (gas_constant * t);
    double xb{};
    if (delta > 50.0)
        xb = std::exp(-delta);
    else if (delta < -50.0)
        xb = 1.0 - std::exp(delta);
    else
        xb = 1.0 / (1.0 + std::exp(delta));
    xb = std::clamp(xb, 1.0e-15, 1.0 - 1.0e-15);
    const double xa = 1.0 - xb;
    const double g = xa * mu_a + xb * mu_b + gas_constant * t * (xa * std::log(xa) + xb * std::log(xb));
    return {xa, xb, g};
}

double langmuir_surface_coverage(double k, double activity) {
    if (!(k >= 0.0) || !(activity >= 0.0) || !std::isfinite(k) || !std::isfinite(activity)) {
        throw std::invalid_argument("invalid Langmuir controls");
    }
    const double ka = k * activity;
    return ka / (1.0 + ka);
}

double binary_ion_exchange_fraction(double selectivity, double activity_a, double activity_b) {
    if (!(selectivity > 0.0) || !(activity_a >= 0.0) || !(activity_b >= 0.0) || !std::isfinite(selectivity) ||
        !std::isfinite(activity_a) || !std::isfinite(activity_b)) {
        throw std::invalid_argument("invalid ion-exchange controls");
    }
    if (activity_b == 0.0)
        return 0.0;
    if (activity_a == 0.0)
        return 1.0;
    const double ratio = selectivity * activity_b / activity_a;
    return ratio / (1.0 + ratio);
}

double ShomatePolynomial::cp_molar(double T) const {
    if (!(T >= minimum_temperature && T <= maximum_temperature))
        throw std::out_of_range("Shomate temperature out of range");
    double t = T / 1000.0;
    return A + B * t + C * t * t + D * t * t * t + E / (t * t);
}
double ShomatePolynomial::h_molar_kj(double T) const {
    if (!(T >= minimum_temperature && T <= maximum_temperature))
        throw std::out_of_range("Shomate temperature out of range");
    double t = T / 1000.0;
    return A * t + B * t * t / 2 + C * t * t * t / 3 + D * t * t * t * t / 4 - E / t + F - H;
}
double ShomatePolynomial::s_molar(double T) const {
    if (!(T >= minimum_temperature && T <= maximum_temperature))
        throw std::out_of_range("Shomate temperature out of range");
    double t = T / 1000.0;
    return A * std::log(t) + B * t + C * t * t / 2 + D * t * t * t / 3 - E / (2 * t * t) + G;
}
double BinaryPitzer11::mean_activity_coefficient(double m, double T) const {
    if (!(m >= 0) || !(T > 0))
        throw std::invalid_argument("invalid Pitzer state");
    if (m == 0)
        return 1.0;
    const double I = m, s = std::sqrt(I), Aphi = 0.392, b = 1.2, alpha = 2.0, x = alpha * s;
    double g = x < 1e-8 ? 1.0 : 2.0 * (1.0 - (1.0 + x) * std::exp(-x)) / (x * x);
    double f = -Aphi * (s / (1.0 + b * s) + (2.0 / b) * std::log(1.0 + b * s));
    double ln_gamma = f + m * (2.0 * beta0 + 2.0 * beta1 * g) + 1.5 * m * m * c_phi;
    return std::exp(ln_gamma);
}

std::vector<double> MulticomponentPitzer::activity_coefficients(std::span<const int> charge,
                                                                std::span<const double> molality,
                                                                double temperature) const {
    if (charge.size() != molality.size() || charge.empty() || !(temperature > 0.0) || !(a_phi_ > 0.0) || !(b_ > 0.0))
        throw std::invalid_argument("invalid multicomponent Pitzer state");
    double ionic_strength = 0.0;
    for (std::size_t i = 0; i < charge.size(); ++i) {
        if (!(molality[i] >= 0.0))
            throw std::invalid_argument("negative Pitzer molality");
        ionic_strength += 0.5 * molality[i] * static_cast<double>(charge[i] * charge[i]);
    }
    const double root = std::sqrt(ionic_strength);
    const double f = -a_phi_ * (root / (1.0 + b_ * root) + (2.0 / b_) * std::log(1.0 + b_ * root));
    std::vector<double> result(charge.size(), 1.0);
    for (std::size_t i = 0; i < charge.size(); ++i) {
        double ln_gamma = static_cast<double>(charge[i] * charge[i]) * f;
        for (const auto& interaction : interactions_) {
            if (interaction.first >= charge.size() || interaction.second >= charge.size() ||
                !(interaction.beta0 >= 0.0) || !(interaction.beta1 >= 0.0) || !(interaction.c_phi >= 0.0))
                throw std::invalid_argument("invalid Pitzer interaction");
            if (i != interaction.first && i != interaction.second)
                continue;
            const std::size_t partner = i == interaction.first ? interaction.second : interaction.first;
            const double x = 2.0 * root;
            const double g = x < 1.0e-12 ? 1.0 : 2.0 * (1.0 - (1.0 + x) * std::exp(-x)) / (x * x);
            ln_gamma += molality[partner] * (2.0 * interaction.beta0 + 2.0 * interaction.beta1 * g) +
                        interaction.c_phi * molality[partner] * molality[partner];
        }
        result[i] = std::exp(ln_gamma);
    }
    return result;
}

} // namespace cfd::chemistry
