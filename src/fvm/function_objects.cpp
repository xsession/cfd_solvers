#include "cfd/fvm/function_objects.hpp"
#include <stdexcept>
namespace cfd::fvm {
void FunctionObjectManager::add(FunctionObject f){if(!f)throw std::invalid_argument("empty function object");objects_.push_back(std::move(f));}
void FunctionObjectManager::execute(const SolverReport&r)const{for(const auto&f:objects_)f(r);}
void ResidualHistory::observe(const SolverReport&r){records_.push_back(r);}
double ResidualHistory::latest(std::string_view field)const{if(records_.empty())throw std::runtime_error("empty residual history");auto it=records_.back().residuals.find(std::string(field));if(it==records_.back().residuals.end())throw std::out_of_range("residual field not found");return it->second;}
} // namespace cfd::fvm
