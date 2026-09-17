#include "cfd/solvers/lbm/diagnostics.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::lbm {
ForceTorque2D momentum_exchange_force(const D2Q9Solver&s,double ox,double oy){const auto&g=s.grid();const auto&solid=s.solid_mask();const auto&f=s.populations();ForceTorque2D out;for(std::size_t y=0;y<g.ny;++y)for(std::size_t x=0;x<g.nx;++x){const auto n=g.index(x,y);if(solid[n])continue;for(int q=1;q<9;++q){const long long xx=static_cast<long long>(x)+D2Q9Descriptor::cx(q),yy=static_cast<long long>(y)+D2Q9Descriptor::cy(q);if(xx<0||yy<0||xx>=static_cast<long long>(g.nx)||yy>=static_cast<long long>(g.ny))continue;const auto nb=g.index(static_cast<std::size_t>(xx),static_cast<std::size_t>(yy));if(!solid[nb])continue;const double fx=2.0*f(static_cast<std::size_t>(q),n)*D2Q9Descriptor::cx(q),fy=2.0*f(static_cast<std::size_t>(q),n)*D2Q9Descriptor::cy(q);out.force_x+=fx;out.force_y+=fy;const double rx=static_cast<double>(x)+.5-ox,ry=static_cast<double>(y)+.5-oy;out.torque_z+=rx*fy-ry*fx;}}return out;}
std::array<double,2> darcy_forchheimer_acceleration(double ux,double uy,double nu,double K,double cf){if(!(nu>=0)||!(K>0)||!(cf>=0))throw std::invalid_argument("invalid porous drag controls");const double speed=std::hypot(ux,uy),factor=-(nu/K+cf*speed/std::sqrt(K));return {factor*ux,factor*uy};}
} // namespace cfd::lbm
