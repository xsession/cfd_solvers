#pragma once
#include "cfd/solvers/fvm/scalar_transport.hpp"
#include <vector>
namespace cfd::fvm {
class MultiSpeciesTransport {
public:
    MultiSpeciesTransport(PolyMesh mesh,std::vector<double> diffusivity,double dt=1e-3);
    void initialize(std::size_t species,const std::function<double(Vec3)>& value);
    void set_face_flux(std::vector<double> flux);
    void step();void run(std::size_t steps);
    [[nodiscard]] std::size_t species_count()const noexcept{return solver_.size();}
    [[nodiscard]] const std::vector<double>& values(std::size_t s)const{return solver_.at(s).values();}
    [[nodiscard]] double inventory(std::size_t s)const{return solver_.at(s).volume_integral();}
private:std::vector<ScalarTransport> solver_;
};
} // namespace cfd::fvm
