#include "cfd/multiphysics/deformation_optics.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::multiphysics {
cfd::optics::SequentialOpticalSystem deform_optical_system_axially(const cfd::optics::SequentialOpticalSystem& system,const std::function<double(double)>& displacement){
    if(!displacement)throw std::invalid_argument("axial optical deformation requires displacement function");
    cfd::optics::SequentialOpticalSystem out(system.object_space_index());
    for(auto s:system.surfaces()){const double dz=displacement(s.vertex_z);if(!std::isfinite(dz))throw std::invalid_argument("non-finite optical deformation");s.vertex_z+=dz;out.add_surface(s);}return out;
}
} // namespace cfd::multiphysics
