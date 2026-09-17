#include "cfd/solvers/optics/paraxial.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::optics {

ParaxialTraceResult trace_paraxial(const SequentialOpticalSystem& system,ParaxialRay ray,double n){
    if(!(n>0.0)) throw std::invalid_argument("paraxial object-space index must be positive");
    std::size_t count=0U;
    for(const auto&s:system.surfaces()){
        const double dz=s.vertex_z-ray.z;
        ray.height+=dz*ray.angle;
        ray.z=s.vertex_z;
        if(std::abs(ray.height)>s.aperture_radius) return {ray,false,count};
        const double n2=s.refractive_index_after;
        if(!(n2>0.0)) return {ray,false,count};
        const double curvature=s.type==SurfaceType::sphere?1.0/s.radius:0.0;
        ray.angle=(n/n2)*ray.angle-((n2-n)/n2)*curvature*ray.height;
        n=n2;++count;
    }
    return {ray,true,count};
}

double paraxial_back_focal_distance(const SequentialOpticalSystem& system,double object_space_index){
    if(system.surfaces().empty()) throw std::invalid_argument("paraxial focal distance requires surfaces");
    ParaxialRay ray{1.0e-3,0.0,system.surfaces().front().vertex_z};
    const auto result=trace_paraxial(system,ray,object_space_index);
    if(!result.valid||std::abs(result.ray.angle)<1.0e-18) return std::numeric_limits<double>::infinity();
    return -result.ray.height/result.ray.angle;
}

} // namespace cfd::optics
