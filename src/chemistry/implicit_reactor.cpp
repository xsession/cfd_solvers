#include "cfd/chemistry/implicit_reactor.hpp"
#include "cfd/core/newton.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfd::chemistry {
namespace {
bool backward_euler(const ReactionNetwork& network,const std::vector<double>& old,
                    double temperature,double dt,std::vector<double>& next,std::size_t iterations) {
    const auto n=old.size();
    next=old;
    const auto residual=[&](std::span<const double> x,std::span<double> r){
        for(double value:x) if(value<0.0||!std::isfinite(value)) {
            std::fill(r.begin(),r.end(),std::numeric_limits<double>::infinity()); return;
        }
        const auto source=network.source_terms(x,temperature);
        for(std::size_t i=0;i<n;++i) r[i]=x[i]-old[i]-dt*source[i];
    };
    const auto jacobian=[&](std::span<const double> x){
        cfd::core::CsrBuilder builder(n,n);
        const auto base=network.source_terms(x,temperature);
        std::vector<double> perturbed(x.begin(),x.end());
        for(std::size_t j=0;j<n;++j){
            const double h=1.0e-7*std::max(1.0,std::abs(x[j]));
            perturbed[j]=x[j]+h;
            const auto source=network.source_terms(perturbed,temperature);
            for(std::size_t i=0;i<n;++i) builder.add(i,j,(i==j?1.0:0.0)-dt*(source[i]-base[i])/h);
            perturbed[j]=x[j];
        }
        return builder.build();
    };
    cfd::core::NewtonConfig control;
    control.max_iterations=iterations;
    control.absolute_tolerance=1.0e-13;
    control.relative_tolerance=1.0e-11;
    control.linear_relative_tolerance=1.0e-13;
    control.line_search_steps=30U;
    control.minimum_step=1.0e-12;
    return cfd::core::newton_solve(next,residual,jacobian,control).converged;
}
}

ReactorResult integrate_isothermal(const ReactionNetwork& network,std::span<double> concentrations,
    double temperature,double duration,const ReactorConfig& config) {
    if(concentrations.empty()||concentrations.size()!=network.species().size()
       ||!(duration>=0.0)||!std::isfinite(duration)||!(config.initial_step>0.0)||!std::isfinite(config.initial_step)
       ||!(config.minimum_step>0.0)||!std::isfinite(config.minimum_step)
       ||config.minimum_step>config.initial_step||!(config.relative_tolerance>0.0)
       ||!std::isfinite(config.relative_tolerance)||!(config.absolute_tolerance>0.0)
       ||!std::isfinite(config.absolute_tolerance)||config.max_steps==0U||config.newton_iterations==0U)
        throw std::invalid_argument("invalid implicit reactor configuration");
    (void)network.source_terms(concentrations,temperature);
    std::vector<double> state(concentrations.begin(),concentrations.end()),full,half,next;
    ReactorResult result;
    double step=std::min(config.initial_step,duration);
    for(std::size_t attempt=0;result.time_advanced<duration;++attempt){
        if(attempt>=config.max_steps) throw std::runtime_error("implicit reactor step budget exhausted");
        step=std::min(step,duration-result.time_advanced);
        if(result.time_advanced+step==result.time_advanced) throw std::runtime_error("implicit reactor time step underflow");
        const bool solved=backward_euler(network,state,temperature,step,full,config.newton_iterations)
            &&backward_euler(network,state,temperature,0.5*step,half,config.newton_iterations)
            &&backward_euler(network,half,temperature,0.5*step,next,config.newton_iterations);
        double error=0.0;
        if(solved) for(std::size_t i=0;i<state.size();++i)
            error=std::max(error,std::abs(full[i]-next[i])/(config.absolute_tolerance
                +config.relative_tolerance*std::max(std::abs(state[i]),std::abs(next[i]))));
        if(solved&&error<=1.0){
            state.swap(next); result.time_advanced+=step; ++result.accepted_steps;
            step*=std::clamp(0.9/std::sqrt(std::max(error,1.0e-12)),0.5,2.0);
        }else{
            ++result.rejected_steps;
            step*=solved?std::clamp(0.9/std::sqrt(error),0.1,0.5):0.5;
            if(step<config.minimum_step) throw std::runtime_error("implicit reactor failed at minimum step");
        }
    }
    std::copy(state.begin(),state.end(),concentrations.begin());
    return result;
}
} // namespace cfd::chemistry
