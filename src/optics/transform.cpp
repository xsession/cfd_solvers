#include "cfd/optics/transform.hpp"
#include <cmath>
namespace cfd::optics { namespace {
Vec3 mul(const std::array<double,9>&m,Vec3 v){return {m[0]*v.x+m[1]*v.y+m[2]*v.z,m[3]*v.x+m[4]*v.y+m[5]*v.z,m[6]*v.x+m[7]*v.y+m[8]*v.z};}
Vec3 tmul(const std::array<double,9>&m,Vec3 v){return {m[0]*v.x+m[3]*v.y+m[6]*v.z,m[1]*v.x+m[4]*v.y+m[7]*v.z,m[2]*v.x+m[5]*v.y+m[8]*v.z};}
}
RigidTransform RigidTransform::from_euler_xyz(Vec3 t,double rx,double ry,double rz){const double cx=std::cos(rx),sx=std::sin(rx),cy=std::cos(ry),sy=std::sin(ry),cz=std::cos(rz),sz=std::sin(rz);RigidTransform q;q.t_=t;q.r_={cz*cy,cz*sy*sx-sz*cx,cz*sy*cx+sz*sx,sz*cy,sz*sy*sx+cz*cx,sz*sy*cx-cz*sx,-sy,cy*sx,cy*cx};return q;}
Vec3 RigidTransform::apply_vector(Vec3 v)const noexcept{return mul(r_,v);} Vec3 RigidTransform::apply_point(Vec3 p)const noexcept{auto v=mul(r_,p);return {v.x+t_.x,v.y+t_.y,v.z+t_.z};}
Vec3 RigidTransform::inverse_vector(Vec3 v)const noexcept{return tmul(r_,v);} Vec3 RigidTransform::inverse_point(Vec3 p)const noexcept{return tmul(r_,{p.x-t_.x,p.y-t_.y,p.z-t_.z});}
Ray RigidTransform::apply(Ray r)const noexcept{r.origin=apply_point(r.origin);r.direction=apply_vector(r.direction);return r;} Ray RigidTransform::inverse(Ray r)const noexcept{r.origin=inverse_point(r.origin);r.direction=inverse_vector(r.direction);return r;}
} // namespace cfd::optics
