#include "cfd/solvers/fvm/reacting.hpp"
#include <stdexcept>
namespace cfd::fvm {
std::vector<std::vector<double>> reacting_cell_sources(const cfd::chemistry::ReactionNetwork& n,std::span<const std::vector<double>> c,std::span<const double> t){if(c.size()!=t.size())throw std::invalid_argument("reacting CFD cell-state size mismatch");std::vector<std::vector<double>>out(c.size());for(std::size_t i=0;i<c.size();++i)out[i]=n.source_terms(c[i],t[i]);return out;}
}
