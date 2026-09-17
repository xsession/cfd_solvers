#include "cfd/solvers/optics/sequential.hpp"
#include "cfd/core/parallel.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::optics { namespace {
Vec3 add(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};} Vec3 scale(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};} Vec3 sub(Vec3 a,Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
bool sag_and_slope(const SequentialSurface&s,double r,double&sag,double&slope){
    if(s.type==SurfaceType::plane){sag=0.0;slope=0.0;return true;} const double c=1.0/s.radius; const double k=s.type==SurfaceType::sphere?0.0:s.conic_constant;
    const double rad=1.0-(1.0+k)*c*c*r*r; if(!(rad>0.0))return false; const double root=std::sqrt(rad); sag=c*r*r/(1.0+root); slope=c*r/root;
    if(s.type==SurfaceType::even_asphere){for(std::size_t i=0;i<s.even_coefficients.size();++i){const int e=4+2*static_cast<int>(i);sag+=s.even_coefficients[i]*std::pow(r,e);slope+=static_cast<double>(e)*s.even_coefficients[i]*std::pow(r,e-1);}}
    return std::isfinite(sag)&&std::isfinite(slope);
}
bool intersect_surface(const SequentialSurface&s,const Ray&ray,Vec3&hit,Vec3&normal){constexpr double eps=1.0e-12;double t=0.0;
    if(s.type==SurfaceType::plane){if(std::abs(ray.direction.z)<eps)return false;t=(s.vertex_z-ray.origin.z)/ray.direction.z;if(!(t>eps))return false;hit=add(ray.origin,scale(ray.direction,t));normal={0,0,-1};}
    else if(s.type==SurfaceType::conic||s.type==SurfaceType::even_asphere){if(std::abs(ray.direction.z)<eps)return false;t=(s.vertex_z-ray.origin.z)/ray.direction.z;bool ok=false;
        for(std::size_t it=0;it<50;++it){hit=add(ray.origin,scale(ray.direction,t));const double r=std::hypot(hit.x,hit.y);double sag=0,slope=0;if(!sag_and_slope(s,r,sag,slope))return false;const double res=hit.z-s.vertex_z-sag;const double sx=r>0?slope*hit.x/r:0,sy=r>0?slope*hit.y/r:0;normal=normalized({sx,sy,-1});if(std::abs(res)<1e-12*(1+std::abs(sag))){ok=true;break;}const double der=ray.direction.z-sx*ray.direction.x-sy*ray.direction.y;if(std::abs(der)<eps)return false;t-=res/der;if(!std::isfinite(t))return false;}if(!ok||!(t>eps))return false;}
    else{if(std::abs(s.radius)<eps)return false;const Vec3 center{0,0,s.vertex_z+s.radius};const Vec3 oc=sub(ray.origin,center);const double b=dot(oc,ray.direction),c=dot(oc,oc)-s.radius*s.radius,disc=b*b-c;if(disc<0)return false;const double root=std::sqrt(disc),t0=-b-root,t1=-b+root;if(t0>eps)t=t0;else if(t1>eps)t=t1;else return false;hit=add(ray.origin,scale(ray.direction,t));normal=normalized(sub(hit,center));}
    return hit.x*hit.x+hit.y*hit.y<=s.aperture_radius*s.aperture_radius;
}
} // namespace

double surface_sag(const SequentialSurface&s,double r){if(!(r>=0)||!std::isfinite(r))throw std::invalid_argument("invalid sag radius");double sag=0,slope=0;if(!sag_and_slope(s,r,sag,slope))throw std::domain_error("surface sag outside real conic");return sag;}
SequentialOpticalSystem::SequentialOpticalSystem(double n):object_space_index_(n){if(!(n>0))throw std::invalid_argument("object-space index must be positive");}
void SequentialOpticalSystem::add_surface(SequentialSurface s){if(!(s.aperture_radius>0)||!(s.refractive_index_after>0)||!std::isfinite(s.vertex_z)||!std::isfinite(s.aperture_radius)||!std::isfinite(s.refractive_index_after)||!std::isfinite(s.conic_constant)||(s.type!=SurfaceType::plane&&(!std::isfinite(s.radius)||s.radius==0)))throw std::invalid_argument("invalid sequential optical surface");for(double a:s.even_coefficients)if(!std::isfinite(a))throw std::invalid_argument("invalid asphere coefficient");surfaces_.push_back(s);}
RayTraceResult SequentialOpticalSystem::trace(Ray ray) const {ray.direction=normalized(ray.direction);double n=object_space_index_;std::size_t count=0;for(const auto&s:surfaces_){Vec3 hit{},normal{};if(!intersect_surface(s,ray,hit,normal))return {ray,false,count};const auto transmitted=refract(ray.direction,normal,n,s.refractive_index_after);if(!transmitted)return {ray,false,count};ray.origin=hit;ray.direction=*transmitted;n=s.refractive_index_after;++count;}return {ray,true,count};}
std::vector<RayTraceResult> SequentialOpticalSystem::trace_many(const std::vector<Ray>&rays) const {std::vector<RayTraceResult> out(rays.size());cfd::core::parallel_for(rays.size(),[&](std::size_t i){out[i]=trace(rays[i]);});return out;}
double rms_spot_radius_at_plane(const std::vector<RayTraceResult>&rays,double z){double sx=0,sy=0;std::size_t n=0;for(const auto&r:rays){if(!r.valid||std::abs(r.ray.direction.z)<=1e-14)continue;const double t=(z-r.ray.origin.z)/r.ray.direction.z;if(t<0)continue;sx+=r.ray.origin.x+t*r.ray.direction.x;sy+=r.ray.origin.y+t*r.ray.direction.y;++n;}if(n==0)throw std::runtime_error("no valid rays reach image plane");const double cx=sx/static_cast<double>(n),cy=sy/static_cast<double>(n);double sum=0;for(const auto&r:rays){if(!r.valid||std::abs(r.ray.direction.z)<=1e-14)continue;const double t=(z-r.ray.origin.z)/r.ray.direction.z;if(t<0)continue;const double x=r.ray.origin.x+t*r.ray.direction.x,y=r.ray.origin.y+t*r.ray.direction.y;sum+=(x-cx)*(x-cx)+(y-cy)*(y-cy);}return std::sqrt(sum/static_cast<double>(n));}
} // namespace cfd::optics
