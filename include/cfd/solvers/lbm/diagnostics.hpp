#pragma once
#include "cfd/solvers/lbm/d2q9.hpp"
namespace cfd::lbm {struct ForceTorque2D {double force_x{},force_y{},torque_z{};};[[nodiscard]] ForceTorque2D momentum_exchange_force(const D2Q9Solver& solver,double origin_x,double origin_y);[[nodiscard]] std::array<double,2> darcy_forchheimer_acceleration(double ux,double uy,double kinematic_viscosity,double permeability,double forchheimer_coefficient=0.0);}
