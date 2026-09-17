#include "cfd/chemistry/reaction_path.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::chemistry {
std::vector<ReactionPathEdge> reaction_path(const ReactionNetwork& network,std::span<const double> c,double t,double minflux){
    if(minflux<0)throw std::invalid_argument("negative reaction-path cutoff");const auto rates=network.reaction_rates(c,t);std::vector<ReactionPathEdge> out;
    for(std::size_t r=0;r<rates.size();++r){const auto&rx=network.reactions()[r];const bool forward=rates[r]>=0;const auto&from=forward?rx.reactants:rx.products;const auto&to=forward?rx.products:rx.reactants;const double rate=std::abs(rates[r]);for(const auto&a:from)for(const auto&b:to){const double flux=rate*a.coefficient*b.coefficient;if(flux>=minflux&&flux>0)out.push_back({a.species,b.species,flux});}}
    return out;
}
}
