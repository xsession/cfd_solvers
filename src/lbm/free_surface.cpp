#include "cfd/solvers/lbm/free_surface.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace cfd::lbm {
namespace {
std::size_t wrap(long long i,std::size_t n){const auto m=static_cast<long long>(n);i%=m;if(i<0)i+=m;return static_cast<std::size_t>(i);}
std::size_t idx(std::size_t x,std::size_t y,std::size_t z,const FreeSurfaceField3D&s){return(z*s.ny+y)*s.nx+x;}
void classify(FreeSurfaceField3D&s,double eps){for(std::size_t i=0;i<s.fill.size();++i){s.fill[i]=std::clamp(s.fill[i],0.0,1.0);s.phase[i]=s.fill[i]<=eps?FreeSurfacePhase::gas:(s.fill[i]>=1.0-eps?FreeSurfacePhase::fluid:FreeSurfacePhase::interface_cell);}}
void validate_velocity(const FreeSurfaceField3D&s,std::span<const float>ux,std::span<const float>uy,std::span<const float>uz){if(ux.size()!=s.fill.size()||uy.size()!=s.fill.size()||uz.size()!=s.fill.size())throw std::invalid_argument("free-surface velocity field size mismatch");}
}

FreeSurfaceField3D make_free_surface_field(std::size_t nx,std::size_t ny,std::size_t nz,std::span<const double> fill,double eps){
    if(nx<2||ny<2||nz<2||fill.size()!=nx*ny*nz||!(eps>=0.0&&eps<0.5))throw std::invalid_argument("invalid free-surface field");
    FreeSurfaceField3D s{nx,ny,nz,std::vector<double>(fill.begin(),fill.end()),std::vector<FreeSurfacePhase>(fill.size())};classify(s,eps);return s;
}

double free_surface_volume(const FreeSurfaceField3D&s) noexcept{return std::accumulate(s.fill.begin(),s.fill.end(),0.0);}

void advect_free_surface_vof(FreeSurfaceField3D&s,std::span<const float>ux,std::span<const float>uy,std::span<const float>uz,double dt,double eps){
    validate_velocity(s,ux,uy,uz);if(!(dt>0.0))throw std::invalid_argument("free-surface dt must be positive");
    const double vmax=[&]{double m=0;for(std::size_t i=0;i<ux.size();++i)m=std::max(m,std::abs(static_cast<double>(ux[i]))+std::abs(static_cast<double>(uy[i]))+std::abs(static_cast<double>(uz[i])));return m;}();
    if(dt*vmax>1.0+1e-12)throw std::invalid_argument("free-surface CFL exceeds one cell per step");
    std::vector<double> delta(s.fill.size(),0.0);
    const auto flux=[&](std::size_t a,std::size_t b,double vel){const double face=0.5*vel;return face>=0.0?face*s.fill[a]:face*s.fill[b];};
    for(std::size_t z=0;z<s.nz;++z)for(std::size_t y=0;y<s.ny;++y)for(std::size_t x=0;x<s.nx;++x){
        const auto a=idx(x,y,z,s),xp=idx(wrap(static_cast<long long>(x)+1,s.nx),y,z,s),yp=idx(x,wrap(static_cast<long long>(y)+1,s.ny),z,s),zp=idx(x,y,wrap(static_cast<long long>(z)+1,s.nz),s);
        const double fx=dt*flux(a,xp,static_cast<double>(ux[a]+ux[xp]));
        const double fy=dt*flux(a,yp,static_cast<double>(uy[a]+uy[yp]));
        const double fz=dt*flux(a,zp,static_cast<double>(uz[a]+uz[zp]));
        delta[a]-=fx+fy+fz;delta[xp]+=fx;delta[yp]+=fy;delta[zp]+=fz;
    }
    for(std::size_t i=0;i<s.fill.size();++i)s.fill[i]+=delta[i];
    classify(s,eps);
}

std::vector<ParticleVector3> free_surface_normals(const FreeSurfaceField3D&s){
    std::vector<ParticleVector3> out(s.fill.size());
    for(std::size_t z=0;z<s.nz;++z)for(std::size_t y=0;y<s.ny;++y)for(std::size_t x=0;x<s.nx;++x){const auto n=idx(x,y,z,s);const double gx=0.5*(s.fill[idx(wrap(static_cast<long long>(x)+1,s.nx),y,z,s)]-s.fill[idx(wrap(static_cast<long long>(x)-1,s.nx),y,z,s)]),gy=0.5*(s.fill[idx(x,wrap(static_cast<long long>(y)+1,s.ny),z,s)]-s.fill[idx(x,wrap(static_cast<long long>(y)-1,s.ny),z,s)]),gz=0.5*(s.fill[idx(x,y,wrap(static_cast<long long>(z)+1,s.nz),s)]-s.fill[idx(x,y,wrap(static_cast<long long>(z)-1,s.nz),s)]);const double m=std::sqrt(gx*gx+gy*gy+gz*gz);if(m>1e-14)out[n]={gx/m,gy/m,gz/m};}
    return out;
}

std::vector<double> free_surface_curvature(const FreeSurfaceField3D&s){
    const auto normals=free_surface_normals(s);std::vector<double> k(s.fill.size(),0.0);
    for(std::size_t z=0;z<s.nz;++z)for(std::size_t y=0;y<s.ny;++y)for(std::size_t x=0;x<s.nx;++x){const auto xp=idx(wrap(static_cast<long long>(x)+1,s.nx),y,z,s),xm=idx(wrap(static_cast<long long>(x)-1,s.nx),y,z,s),yp=idx(x,wrap(static_cast<long long>(y)+1,s.ny),z,s),ym=idx(x,wrap(static_cast<long long>(y)-1,s.ny),z,s),zp=idx(x,y,wrap(static_cast<long long>(z)+1,s.nz),s),zm=idx(x,y,wrap(static_cast<long long>(z)-1,s.nz),s),n=idx(x,y,z,s);k[n]=-0.5*((normals[xp].x-normals[xm].x)+(normals[yp].y-normals[ym].y)+(normals[zp].z-normals[zm].z));}
    return k;
}

std::vector<ParticleVector3> free_surface_surface_tension_acceleration(const FreeSurfaceField3D&s,double sigma,double density){
    if(!(sigma>=0.0)||!(density>0.0)) throw std::invalid_argument("invalid free-surface material properties");
    const auto k=free_surface_curvature(s);
    std::vector<ParticleVector3>a(s.fill.size());
    for(std::size_t z=0;z<s.nz;++z)for(std::size_t y=0;y<s.ny;++y)for(std::size_t x=0;x<s.nx;++x){const auto n=idx(x,y,z,s);const double gx=0.5*(s.fill[idx(wrap(static_cast<long long>(x)+1,s.nx),y,z,s)]-s.fill[idx(wrap(static_cast<long long>(x)-1,s.nx),y,z,s)]),gy=0.5*(s.fill[idx(x,wrap(static_cast<long long>(y)+1,s.ny),z,s)]-s.fill[idx(x,wrap(static_cast<long long>(y)-1,s.ny),z,s)]),gz=0.5*(s.fill[idx(x,y,wrap(static_cast<long long>(z)+1,s.nz),s)]-s.fill[idx(x,y,wrap(static_cast<long long>(z)-1,s.nz),s)]);const double c=sigma*k[n]/density;a[n]={c*gx,c*gy,c*gz};}
    return a;
}
} // namespace cfd::lbm
