#include "cfd/electrochemistry/electrochemistry.hpp"
#include "cfd/chemistry/kinetics.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::electrochemistry {

double nernst_potential(double standard_potential,
                         double temperature,
                         double electrons,
                         double reaction_quotient) {
    if (!(temperature > 0.0) || !(electrons > 0.0) || !(reaction_quotient > 0.0)) {
        throw std::invalid_argument("invalid Nernst parameters");
    }
    return standard_potential
         - (cfd::chemistry::gas_constant * temperature / (electrons * faraday_constant))
         * std::log(reaction_quotient);
}

double butler_volmer_current_density(double exchange_current_density,
                                      double overpotential,
                                      double temperature,
                                      double electrons,
                                      double anodic_transfer,
                                      double cathodic_transfer) {
    if (exchange_current_density < 0.0 || !(temperature > 0.0) || !(electrons > 0.0)
        || anodic_transfer < 0.0 || cathodic_transfer < 0.0
        || !(anodic_transfer + cathodic_transfer > 0.0)) {
        throw std::invalid_argument("invalid Butler-Volmer parameters");
    }
    const double scale = electrons * faraday_constant
                       / (cfd::chemistry::gas_constant * temperature);
    const double anodic = std::exp(anodic_transfer * scale * overpotential);
    const double cathodic = std::exp(-cathodic_transfer * scale * overpotential);
    return exchange_current_density * (anodic - cathodic);
}

double faradaic_molar_flux(double current_density, double electrons) {
    if (!(electrons > 0.0)) throw std::invalid_argument("electron count must be positive");
    return current_density / (electrons * faraday_constant);
}

double corrosion_penetration_rate(double current_density,
                                   double molar_mass,
                                   double metal_density,
                                   double electrons) {
    if (!(molar_mass > 0.0) || !(metal_density > 0.0)) {
        throw std::invalid_argument("invalid corrosion material parameters");
    }
    return faradaic_molar_flux(current_density, electrons) * molar_mass / metal_density;
}

} // namespace cfd::electrochemistry
