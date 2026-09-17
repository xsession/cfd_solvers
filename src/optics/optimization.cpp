#include "cfd/optics/optimization.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace cfd::optics {
OptimizationResult coordinate_descent(MeritFunction merit,std::vector<double>x,std::vector<double>step,std::size_t maxit,double tol){if(!merit||x.empty()||step.size()!=x.size()||maxit==0||!(tol>0))throw std::invalid_argument("invalid optical optimizer controls");double best=merit(x);for(std::size_t it=0;it<maxit;++it){bool improved=false;for(std::size_t j=0;j<x.size();++j){auto trial=x;trial[j]+=step[j];double m=merit(trial);if(m<best){x=std::move(trial);best=m;improved=true;continue;}trial=x;trial[j]-=step[j];m=merit(trial);if(m<best){x=std::move(trial);best=m;improved=true;}}if(!improved){double mx=0;for(auto&h:step){h*=0.5;mx=std::max(mx,std::abs(h));}if(mx<tol)return {x,best,it+1,true};}}return {x,best,maxit,false};}
std::vector<double> finite_difference_sensitivity(const MeritFunction& merit,std::span<const double> p,double rel){if(!merit||!(rel>0))throw std::invalid_argument("invalid sensitivity controls");std::vector<double>x(p.begin(),p.end()),g(p.size());for(std::size_t i=0;i<p.size();++i){const double h=rel*std::max(1.0,std::abs(p[i]));x[i]=p[i]+h;const double a=merit(x);x[i]=p[i]-h;const double b=merit(x);x[i]=p[i];g[i]=(a-b)/(2*h);}return g;}
MonteCarloStats monte_carlo_tolerance(const MeritFunction& merit,std::span<const double> nominal,std::span<const double> sigma,std::size_t samples,std::uint64_t seed){if(!merit||nominal.size()!=sigma.size()||samples==0)throw std::invalid_argument("invalid tolerance controls");std::mt19937_64 rng(seed);std::normal_distribution<double> normal;std::vector<double>x(nominal.begin(),nominal.end());MonteCarloStats s;double sum=0,sum2=0;s.minimum=std::numeric_limits<double>::infinity();s.maximum=-s.minimum;for(std::size_t k=0;k<samples;++k){for(std::size_t i=0;i<x.size();++i)x[i]=nominal[i]+sigma[i]*normal(rng);const double m=merit(x);sum+=m;sum2+=m*m;s.minimum=std::min(s.minimum,m);s.maximum=std::max(s.maximum,m);}s.samples=samples;s.mean=sum/samples;s.standard_deviation=std::sqrt(std::max(0.0,sum2/samples-s.mean*s.mean));return s;}

OptimizationResult differential_evolution(MeritFunction merit,std::span<const double> lower,std::span<const double> upper,std::size_t population,std::size_t generations,double mutation,double crossover,std::uint64_t seed){
    if(!merit||lower.empty()||lower.size()!=upper.size()||population<4||generations==0||!(mutation>0.0)||!(crossover>=0.0&&crossover<=1.0))throw std::invalid_argument("invalid differential-evolution controls");
    for(std::size_t j=0;j<lower.size();++j)if(!(lower[j]<upper[j])||!std::isfinite(lower[j])||!std::isfinite(upper[j]))throw std::invalid_argument("invalid differential-evolution bounds");
    std::mt19937_64 rng(seed);std::uniform_real_distribution<double> unit(0.0,1.0);std::uniform_int_distribution<std::size_t> pick(0,population-1);
    std::vector<std::vector<double>> pop(population,std::vector<double>(lower.size()));std::vector<double> score(population);
    for(std::size_t i=0;i<population;++i){for(std::size_t j=0;j<lower.size();++j)pop[i][j]=lower[j]+unit(rng)*(upper[j]-lower[j]);score[i]=merit(pop[i]);}
    std::size_t best=static_cast<std::size_t>(std::min_element(score.begin(),score.end())-score.begin());
    for(std::size_t gen=0;gen<generations;++gen){
        for(std::size_t i=0;i<population;++i){std::size_t a,b,c;do{a=pick(rng);}while(a==i);do{b=pick(rng);}while(b==i||b==a);do{c=pick(rng);}while(c==i||c==a||c==b);std::vector<double> trial=pop[i];const std::size_t forced=pick(rng)%lower.size();for(std::size_t j=0;j<lower.size();++j)if(unit(rng)<crossover||j==forced)trial[j]=std::clamp(pop[a][j]+mutation*(pop[b][j]-pop[c][j]),lower[j],upper[j]);const double value=merit(trial);if(value<score[i]){pop[i]=std::move(trial);score[i]=value;if(score[i]<score[best])best=i;}}
    }
    return {pop[best],score[best],generations,true};
}
} // namespace cfd::optics
