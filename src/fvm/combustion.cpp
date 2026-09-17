#include "cfd/fvm/combustion.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
namespace {void check(std::span<const double> y,std::span<const IdealGasSpeciesThermo>s){if(y.empty()||y.size()!=s.size())throw std::invalid_argument("mixture size mismatch");double sum=0;for(std::size_t i=0;i<y.size();++i){if(y[i]<0||!std::isfinite(y[i])||!(s[i].molar_mass>0)||!(s[i].cp_mass>0))throw std::invalid_argument("invalid mixture thermo");sum+=y[i];}if(std::abs(sum-1.0)>1e-10)throw std::invalid_argument("mass fractions must sum to one");}}
CombustionThermoState ideal_gas_mixture_state(double T,double p,std::span<const double> y,std::span<const IdealGasSpeciesThermo>s){check(y,s);if(!(T>0)||!(p>0))throw std::invalid_argument("invalid thermodynamic state");double invM=0,cp=0,hf=0;for(std::size_t i=0;i<y.size();++i){invM+=y[i]/s[i].molar_mass;cp+=y[i]*s[i].cp_mass;hf+=y[i]*s[i].formation_enthalpy_mass;}const double M=1.0/invM,R=cfd::chemistry::gas_constant/M,cv=cp-R;if(!(cv>0))throw std::invalid_argument("mixture cp must exceed gas constant");return {T,p,p/(R*T),M,R,cp,cv,cp/cv,hf+cp*T};}
double temperature_from_mixture_enthalpy(double h,std::span<const double> y,std::span<const IdealGasSpeciesThermo>s){check(y,s);double cp=0,hf=0;for(std::size_t i=0;i<y.size();++i){cp+=y[i]*s[i].cp_mass;hf+=y[i]*s[i].formation_enthalpy_mass;}const double T=(h-hf)/cp;if(!(T>0)||!std::isfinite(T))throw std::invalid_argument("mixture enthalpy implies invalid temperature");return T;}
} // namespace cfd::fvm
