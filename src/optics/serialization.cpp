#include "cfd/optics/serialization.hpp"
#include <array>
#include <charconv>
#include <iomanip>
#include <regex>
#include <sstream>
#include <stdexcept>
namespace cfd::optics {
namespace {
const char* type_name(SurfaceType t){switch(t){case SurfaceType::plane:return "plane";case SurfaceType::sphere:return "sphere";case SurfaceType::conic:return "conic";case SurfaceType::even_asphere:return "even_asphere";}return "plane";}
SurfaceType parse_type(const std::string&s){if(s=="plane")return SurfaceType::plane;if(s=="sphere")return SurfaceType::sphere;if(s=="conic")return SurfaceType::conic;if(s=="even_asphere")return SurfaceType::even_asphere;throw std::invalid_argument("unknown optical surface type");}
}
std::string sequential_system_to_json(const SequentialOpticalSystem& system,double object_index){std::ostringstream o;o<<std::setprecision(17)<<"{\"object_index\":"<<object_index<<",\"surfaces\":[";bool first=true;for(const auto&s:system.surfaces()){if(!first)o<<',';first=false;o<<"{\"type\":\""<<type_name(s.type)<<"\",\"z\":"<<s.vertex_z<<",\"radius\":"<<s.radius<<",\"aperture\":"<<s.aperture_radius<<",\"n_after\":"<<s.refractive_index_after<<",\"k\":"<<s.conic_constant<<",\"a\":[";for(std::size_t i=0;i<s.even_coefficients.size();++i){if(i)o<<',';o<<s.even_coefficients[i];}o<<"]}";}o<<"]}";return o.str();}
SequentialOpticalSystem sequential_system_from_json(const std::string&json){std::smatch m;std::regex idx("\\\"object_index\\\":([-+0-9.eE]+)");if(!std::regex_search(json,m,idx))throw std::invalid_argument("missing optical object index");double object_index=std::stod(m[1].str());SequentialOpticalSystem sys(object_index);std::regex surf("\\{\\\"type\\\":\\\"([^\\\"]+)\\\",\\\"z\\\":([-+0-9.eE]+),\\\"radius\\\":([-+0-9.eE]+),\\\"aperture\\\":([-+0-9.eE]+),\\\"n_after\\\":([-+0-9.eE]+),\\\"k\\\":([-+0-9.eE]+),\\\"a\\\":\\[([-+0-9.eE]+),([-+0-9.eE]+),([-+0-9.eE]+),([-+0-9.eE]+)\\]\\}");for(auto it=std::sregex_iterator(json.begin(),json.end(),surf);it!=std::sregex_iterator();++it){const auto&q=*it;SequentialSurface s;s.type=parse_type(q[1].str());s.vertex_z=std::stod(q[2].str());s.radius=std::stod(q[3].str());s.aperture_radius=std::stod(q[4].str());s.refractive_index_after=std::stod(q[5].str());s.conic_constant=std::stod(q[6].str());for(std::size_t i=0;i<4;++i)s.even_coefficients[i]=std::stod(q[7+static_cast<int>(i)].str());sys.add_surface(s);}return sys;}
} // namespace cfd::optics
