#include "cfd/solvers/fvm/species_transport.hpp"
#include <stdexcept>
namespace cfd::fvm {
MultiSpeciesTransport::MultiSpeciesTransport(PolyMesh mesh,std::vector<double>d,double dt){if(d.empty())throw std::invalid_argument("species transport requires species");solver_.reserve(d.size());for(double x:d){if(!(x>=0))throw std::invalid_argument("invalid species diffusivity");ScalarTransportConfig c;c.dt=dt;c.diffusivity=x;solver_.emplace_back(mesh,c);}}
void MultiSpeciesTransport::initialize(std::size_t s,const std::function<double(Vec3)>&v){solver_.at(s).initialize(v);}void MultiSpeciesTransport::set_face_flux(std::vector<double>f){for(auto&s:solver_)s.set_face_flux(f);}void MultiSpeciesTransport::step(){for(auto&s:solver_)s.step();}void MultiSpeciesTransport::run(std::size_t n){for(std::size_t i=0;i<n;++i)step();}
} // namespace cfd::fvm
