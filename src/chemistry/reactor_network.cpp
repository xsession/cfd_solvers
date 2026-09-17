#include "cfd/chemistry/reactor_network.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::chemistry {
WellStirredReactorNetwork::WellStirredReactorNetwork(ReactionNetwork chemistry):chemistry_(std::move(chemistry)){}
std::size_t WellStirredReactorNetwork::add_reactor(WellStirredReactor reactor){
    if(!(reactor.volume>0.0)||!(reactor.temperature>0.0)||!std::isfinite(reactor.volume)||!std::isfinite(reactor.temperature)
       ||reactor.concentrations.size()!=chemistry_.species().size())throw std::invalid_argument("invalid well-stirred reactor");
    for(double c:reactor.concentrations)if(c<0.0||!std::isfinite(c))throw std::invalid_argument("invalid reactor concentration");
    reactors_.push_back(std::move(reactor));return reactors_.size()-1U;
}
void WellStirredReactorNetwork::connect(ReactorConnection c){
    if(c.first>=reactors_.size()||c.second>=reactors_.size()||c.first==c.second||c.exchange_flow<0.0||!std::isfinite(c.exchange_flow))
        throw std::invalid_argument("invalid reactor connection");connections_.push_back(c);
}
void WellStirredReactorNetwork::advance(double duration,double maximum_step,const ReactorConfig& controls){
    if(duration<0.0||!(maximum_step>0.0)||!std::isfinite(duration)||!std::isfinite(maximum_step))throw std::invalid_argument("invalid reactor-network step");
    double time=0.0;const std::size_t ns=chemistry_.species().size();
    while(time<duration){const double dt=std::min(maximum_step,duration-time);std::vector<std::vector<double>> delta(reactors_.size(),std::vector<double>(ns));
        for(const auto& c:connections_){const auto&a=reactors_[c.first];const auto&b=reactors_[c.second];
            for(std::size_t s=0;s<ns;++s){const double amount_rate=c.exchange_flow*(b.concentrations[s]-a.concentrations[s]);delta[c.first][s]+=amount_rate*dt/a.volume;delta[c.second][s]-=amount_rate*dt/b.volume;}}
        for(std::size_t r=0;r<reactors_.size();++r)for(std::size_t s=0;s<ns;++s){reactors_[r].concentrations[s]+=delta[r][s];if(reactors_[r].concentrations[s]<-1e-14)throw std::runtime_error("reactor mixing produced negative concentration");reactors_[r].concentrations[s]=std::max(0.0,reactors_[r].concentrations[s]);}
        for(auto& r:reactors_)if(ns>0U)integrate_isothermal(chemistry_,r.concentrations,r.temperature,dt,controls);
        time+=dt;
    }
}
std::vector<double> WellStirredReactorNetwork::total_species_amounts()const{
    std::vector<double> amount(chemistry_.species().size());for(const auto&r:reactors_)for(std::size_t s=0;s<amount.size();++s)amount[s]+=r.volume*r.concentrations[s];return amount;
}
} // namespace cfd::chemistry
