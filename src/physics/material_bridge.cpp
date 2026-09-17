#include "cfd/physics/material_bridge.hpp"
#include <stdexcept>
namespace cfd::physics {
cfd::fdtd::Maxwell3DConfig apply_to_fdtd(cfd::fdtd::Maxwell3DConfig c,const Material&m){validate_material(m);c.epsilon_rx=m.relative_permittivity[0];c.epsilon_ry=m.relative_permittivity[1];c.epsilon_rz=m.relative_permittivity[2];return c;}
double fem_conductivity(const Material&m){validate_material(m);return m.electrical_conductivity;}
double optical_refractive_index(const Material&m){validate_material(m);if(!(m.refractive_index>0))throw std::invalid_argument("material has no optical refractive index");return m.refractive_index;}
}
