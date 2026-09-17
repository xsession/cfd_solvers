#include "cfd/fem/nonlinear_geometry.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::fem {NonlinearTrussState green_lagrange_truss(double L,double l,double A,double E){if(!(L>0)||!(l>0)||!(A>0)||!(E>0))throw std::invalid_argument("invalid nonlinear truss state");double lambda=l/L,eps=.5*(lambda*lambda-1.0),S=E*eps,force=A*S*lambda,tangent=A*E*(1.5*lambda*lambda-.5)/L;return {lambda,eps,S,force,tangent};}}
