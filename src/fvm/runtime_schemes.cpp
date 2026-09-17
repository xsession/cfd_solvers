#include "cfd/fvm/runtime_schemes.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
FaceInterpolationScheme parse_face_interpolation_scheme(std::string_view s){if(s=="linear")return FaceInterpolationScheme::linear;if(s=="upwind")return FaceInterpolationScheme::upwind;if(s=="boundedLinear"||s=="bounded_linear")return FaceInterpolationScheme::bounded_linear;if(s=="minmod")return FaceInterpolationScheme::minmod;if(s=="vanLeer"||s=="van_leer")return FaceInterpolationScheme::van_leer;throw std::invalid_argument("unknown FVM interpolation scheme");}
const char* face_interpolation_scheme_name(FaceInterpolationScheme s) noexcept{switch(s){case FaceInterpolationScheme::linear:return "linear";case FaceInterpolationScheme::upwind:return "upwind";case FaceInterpolationScheme::bounded_linear:return "boundedLinear";case FaceInterpolationScheme::minmod:return "minmod";case FaceInterpolationScheme::van_leer:return "vanLeer";}return "unknown";}
double courant_limited_timestep(const PolyMesh& mesh,std::span<const double> flux,double target,double maxdt){if(flux.size()!=mesh.face_count()||!(target>0)||!(maxdt>0))throw std::invalid_argument("invalid Courant timestep controls");double dt=maxdt;for(std::size_t c=0;c<mesh.cell_count();++c){double outgoing=0;for(auto fi:mesh.cell_faces()[c]){const auto&f=mesh.faces()[fi];const double phi=f.owner==c?flux[fi]:-flux[fi];outgoing+=std::abs(phi);}if(outgoing>0)dt=std::min(dt,target*mesh.cells()[c].volume/outgoing);}return dt;}
} // namespace cfd::fvm
