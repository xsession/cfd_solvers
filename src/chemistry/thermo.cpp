#include "cfd/chemistry/thermo.hpp"
#include "cfd/chemistry/kinetics.hpp"

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>

namespace cfd::chemistry {
namespace {
std::vector<double> normalized_fractions(std::span<const double> x,std::size_t expected) {
    if(x.size()!=expected || x.empty()) throw std::invalid_argument("thermodynamic composition size mismatch");
    double sum=0.0;
    for(double v:x){if(v<0.0||!std::isfinite(v))throw std::invalid_argument("invalid mole fraction");sum+=v;}
    if(!(sum>0.0)||!std::isfinite(sum))throw std::invalid_argument("empty thermodynamic composition");
    std::vector<double> out(x.size());
    for(std::size_t i=0;i<x.size();++i)out[i]=x[i]/sum;
    return out;
}
void validate_state(double temperature,double pressure){
    if(!(temperature>0.0)||!(pressure>0.0)||!std::isfinite(temperature)||!std::isfinite(pressure))
        throw std::invalid_argument("invalid thermodynamic state");
}
std::vector<double> chemical_potentials_from_activities(const std::vector<ThermoSpecies>& species,
                                                         std::span<const double> activity,
                                                         double temperature){
    std::vector<double> mu(species.size());
    for(std::size_t i=0;i<species.size();++i){
        if(!(activity[i]>0.0))throw std::domain_error("chemical potential undefined for zero activity");
        mu[i]=species[i].standard_chemical_potential+gas_constant*temperature*std::log(activity[i]);
    }
    return mu;
}
}

IdealGasPhase::IdealGasPhase(std::vector<ThermoSpecies> species,double reference_pressure)
    :species_(std::move(species)),reference_pressure_(reference_pressure){
    if(species_.empty()||!(reference_pressure_>0.0)||!std::isfinite(reference_pressure_))
        throw std::invalid_argument("invalid ideal-gas phase");
}
std::vector<double> IdealGasPhase::activities(std::span<const double>x,double temperature,double pressure)const{
    validate_state(temperature,pressure);auto a=normalized_fractions(x,species_.size());
    for(double& v:a)v*=pressure/reference_pressure_;return a;
}
std::vector<double> IdealGasPhase::chemical_potentials(std::span<const double>x,double temperature,double pressure)const{
    const auto a=activities(x,temperature,pressure);return chemical_potentials_from_activities(species_,a,temperature);
}
IdealSolutionPhase::IdealSolutionPhase(std::vector<ThermoSpecies> species):species_(std::move(species)){
    if(species_.empty())throw std::invalid_argument("invalid ideal-solution phase");
}
std::vector<double> IdealSolutionPhase::activities(std::span<const double>x,double temperature,double pressure)const{
    validate_state(temperature,pressure);return normalized_fractions(x,species_.size());
}
std::vector<double> IdealSolutionPhase::chemical_potentials(std::span<const double>x,double temperature,double pressure)const{
    const auto a=activities(x,temperature,pressure);return chemical_potentials_from_activities(species_,a,temperature);
}
std::vector<ThermoSpecies> load_thermo_csv(const std::filesystem::path& path){
    std::ifstream in(path);if(!in)throw std::runtime_error("cannot open thermochemical CSV");std::vector<ThermoSpecies> out;std::string line;
    while(std::getline(in,line)){if(line.empty()||line[0]=='#')continue;std::stringstream ss(line);std::string name,value;if(!std::getline(ss,name,',')||!std::getline(ss,value))throw std::runtime_error("invalid thermochemical CSV row");double mu=std::stod(value);if(name.empty()||!std::isfinite(mu))throw std::runtime_error("invalid thermochemical CSV value");out.push_back({name,mu});}
    if(out.empty())throw std::runtime_error("empty thermochemical CSV");return out;
}

} // namespace cfd::chemistry
