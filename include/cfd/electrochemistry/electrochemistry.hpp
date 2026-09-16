#pragma once

namespace cfd::electrochemistry {

inline constexpr double faraday_constant = 96485.33212; // C/mol

[[nodiscard]] double nernst_potential(double standard_potential,
                                      double temperature,
                                      double electrons,
                                      double reaction_quotient);

[[nodiscard]] double butler_volmer_current_density(double exchange_current_density,
                                                    double overpotential,
                                                    double temperature,
                                                    double electrons = 1.0,
                                                    double anodic_transfer = 0.5,
                                                    double cathodic_transfer = 0.5);

[[nodiscard]] double faradaic_molar_flux(double current_density, double electrons = 1.0);

[[nodiscard]] double corrosion_penetration_rate(double current_density,
                                                 double molar_mass,
                                                 double metal_density,
                                                 double electrons = 1.0);

} // namespace cfd::electrochemistry
