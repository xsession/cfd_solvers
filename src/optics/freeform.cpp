#include "cfd/optics/freeform.hpp"
#include <cmath>
#include <stdexcept>

namespace cfd::optics {
namespace {
double ipow(double x, unsigned p) { double v=1.0; for(unsigned i=0;i<p;++i)v*=x; return v; }
}
PolynomialFreeformSurface::PolynomialFreeformSurface(double z,std::vector<PolynomialTerm> terms,double aperture)
    :vertex_z_(z),terms_(std::move(terms)),aperture_radius_(aperture){
    if(!std::isfinite(z)||!(aperture>0.0)||!std::isfinite(aperture))throw std::invalid_argument("invalid freeform surface");
    for(const auto&t:terms_)if(!std::isfinite(t.coefficient))throw std::invalid_argument("invalid freeform coefficient");
}
double PolynomialFreeformSurface::sag(double x,double y)const{double z=0.0;for(const auto&t:terms_)z+=t.coefficient*ipow(x,t.x_power)*ipow(y,t.y_power);return z;}
Vec3 PolynomialFreeformSurface::normal(double x,double y)const{double sx=0,sy=0;for(const auto&t:terms_){if(t.x_power)sx+=t.coefficient*static_cast<double>(t.x_power)*ipow(x,t.x_power-1)*ipow(y,t.y_power);if(t.y_power)sy+=t.coefficient*static_cast<double>(t.y_power)*ipow(x,t.x_power)*ipow(y,t.y_power-1);}return normalized({sx,sy,-1.0});}
bool PolynomialFreeformSurface::intersect(const Ray& ray,Vec3& point,Vec3& n,std::size_t maxit,double tol)const{
    if(maxit==0||!(tol>0.0)||std::abs(ray.direction.z)<1e-15)return false;
    double t=(vertex_z_-ray.origin.z)/ray.direction.z;
    for(std::size_t it=0;it<maxit;++it){
        point={ray.origin.x+t*ray.direction.x,ray.origin.y+t*ray.direction.y,ray.origin.z+t*ray.direction.z};
        if(point.x*point.x+point.y*point.y>aperture_radius_*aperture_radius_)return false;
        n=normal(point.x,point.y);const double f=point.z-vertex_z_-sag(point.x,point.y);
        if(std::abs(f)<tol)return t>0.0;
        const double denom=-(n.x*ray.direction.x+n.y*ray.direction.y+n.z*ray.direction.z)/(-n.z);
        if(std::abs(denom)<1e-15)return false;t-=f/denom;if(!std::isfinite(t))return false;
    }
    return false;
}
} // namespace cfd::optics
