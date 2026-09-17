#include "cfd/workflow/runner.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>
namespace cfd::workflow {
std::vector<SweepResult> parameter_sweep(std::span<const double>v,const std::function<double(double)>&f){if(!f)throw std::invalid_argument("parameter sweep requires objective");std::vector<SweepResult>out;out.reserve(v.size());for(double x:v)out.push_back({x,f(x)});return out;}
InverseResult inverse_coordinate_search(const std::function<double(std::span<const double>)>&f,std::vector<double>x,std::vector<double>step,std::size_t maxit,double tol){if(!f||x.empty()||x.size()!=step.size()||maxit==0||!(tol>0))throw std::invalid_argument("invalid inverse search");double best=f(x);InverseResult out{x,best,0,false};for(std::size_t it=0;it<maxit;++it){bool improved=false;for(std::size_t j=0;j<x.size();++j){for(double sign:{-1.0,1.0}){auto trial=x;trial[j]+=sign*step[j];double val=f(trial);if(val<best){x=std::move(trial);best=val;improved=true;}}}if(!improved)for(double&s:step)s*=0.5;double maxs=*std::max_element(step.begin(),step.end());out={x,best,it+1,maxs<tol};if(out.converged)break;}return out;}
std::size_t RestartableWorkflow::add(std::string name,Action action,std::vector<std::size_t>deps){if(name.empty()||!action)throw std::invalid_argument("invalid workflow node");for(auto d:deps)if(d>=nodes_.size())throw std::invalid_argument("workflow dependency must precede node");nodes_.push_back({std::move(name),std::move(action),std::move(deps)});return nodes_.size()-1;}
void RestartableWorkflow::run(const std::filesystem::path&checkpoint){std::set<std::string>done;{std::ifstream in(checkpoint);std::string n;while(std::getline(in,n))if(!n.empty())done.insert(n);}completed_names_.assign(done.begin(),done.end());for(std::size_t i=0;i<nodes_.size();++i){auto&n=nodes_[i];if(done.contains(n.name))continue;for(auto d:n.deps)if(!done.contains(nodes_[d].name))throw std::runtime_error("workflow dependency incomplete");n.action();done.insert(n.name);completed_names_.push_back(n.name);std::ofstream out(checkpoint,std::ios::trunc);for(const auto&name:done)out<<name<<'\n';}}
}
