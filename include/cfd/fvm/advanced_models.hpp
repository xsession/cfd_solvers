#pragma once
#include "cfd/fvm/poly_mesh.hpp"
namespace cfd::fvm {
[[nodiscard]] double spalart_allmaras_eddy_viscosity(double molecular_viscosity,double nu_tilde,double cv1=7.1);
[[nodiscard]] double k_epsilon_eddy_viscosity(double density,double k,double epsilon,double c_mu=0.09);
[[nodiscard]] double k_omega_sst_eddy_viscosity(double density,double k,double omega,double strain_rate,double f2=1.0,double a1=0.31);
[[nodiscard]] double des_length_scale(double wall_distance,double grid_scale,double c_des=0.65);
struct ReynoldsStress2D{double xx{},yy{},xy{};};
[[nodiscard]] ReynoldsStress2D boussinesq_reynolds_stress(double density,double k,double eddy_viscosity,double du_dx,double dv_dy,double du_dy_plus_dv_dx);
struct PhaseChangeState{double liquid_fraction{},effective_cp{},enthalpy{};};
[[nodiscard]] PhaseChangeState linear_mushy_phase_change(double temperature,double solidus,double liquidus,double cp_solid,double cp_liquid,double latent_heat);
[[nodiscard]] double drift_flux_slip_velocity(double density_continuous,double density_disperse,double particle_diameter,double viscosity,double gravity=9.81);
[[nodiscard]] double euler_euler_drag_source(double continuous_velocity,double dispersed_velocity,double drag_coefficient,double interfacial_area);
[[nodiscard]] double lubrication_thin_film_flux(double thickness,double pressure_gradient,double viscosity,double surface_velocity=0.0);
} // namespace cfd::fvm
