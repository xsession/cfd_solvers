#include "cfd/multiphysics/electrochemical_heat.hpp"
#include <stdexcept>
namespace cfd::multiphysics {std::vector<double> electrochemical_heat_source(std::span<const double>j,std::span<const double>eta,std::span<const double>ohm){if(j.size()!=eta.size()||(!ohm.empty()&&ohm.size()!=j.size()))throw std::invalid_argument("electrochemical heat size mismatch");std::vector<double>q(j.size());for(std::size_t i=0;i<q.size();++i)q[i]=j[i]*eta[i]+(ohm.empty()?0.0:ohm[i]);return q;}}
