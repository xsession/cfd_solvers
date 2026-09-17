#include "cfd/fvm/radiation.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {double gray_surface_exchange(double ta,double ea,double tb,double eb,double vf){if(ta<0||tb<0||!(ea>0&&ea<=1)||!(eb>0&&eb<=1)||vf<0||vf>1)throw std::invalid_argument("invalid gray radiation");double effective=1.0/(1.0/ea+1.0/eb-1.0);return vf*effective*stefan_boltzmann*(std::pow(ta,4)-std::pow(tb,4));}double optically_thin_radiation_source(double t,double env,double k){if(t<0||env<0||k<0)throw std::invalid_argument("invalid optically-thin radiation");return -4.0*k*stefan_boltzmann*(std::pow(t,4)-std::pow(env,4));}}
