#include "cfd/multiphysics/cht.hpp"
#include <stdexcept>
namespace cfd::multiphysics {InterfaceHeatTransfer conjugate_interface_heat_flux(double ta,double ka,double da,double tb,double kb,double db){if(!(ka>0)||!(kb>0)||!(da>0)||!(db>0))throw std::invalid_argument("invalid CHT interface");double ra=da/ka,rb=db/kb,q=(ta-tb)/(ra+rb);return {q,ta-q*ra};}}
