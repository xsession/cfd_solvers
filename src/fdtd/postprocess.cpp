#include "cfd/solvers/fdtd/postprocess.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
namespace cfd::fdtd {namespace {
using P=cfd::fem::Point3;P normed(P a){double n=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);if(!(n>0))throw std::invalid_argument("zero far-field direction");return {a.x/n,a.y/n,a.z/n};}
P cross(P a,P b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
ComplexVec3 cross(P a,ComplexVec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
ComplexVec3 cross(ComplexVec3 a,P b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
ComplexVec3 add(ComplexVec3 a,ComplexVec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}ComplexVec3 scale(ComplexVec3 a,std::complex<double>s){return {a.x*s,a.y*s,a.z*s};}
double abs2(ComplexVec3 a){return std::norm(a.x)+std::norm(a.y)+std::norm(a.z);}std::size_t id(std::size_t i,std::size_t j,std::size_t k,const SarGrid3D&g){return (k*g.ny+j)*g.nx+i;}
}
FarFieldSample nf2ff(std::span<const NearFieldSurfaceSample>s,P dir,double f,double eta){if(s.empty()||!(f>0)||!(eta>0))throw std::invalid_argument("invalid NF2FF inputs");const P r=normed(dir);constexpr double c0=299792458.0;const double k=2*std::numbers::pi*f/c0;ComplexVec3 sum{};for(const auto&q:s){if(!(q.area>0))throw std::invalid_argument("invalid NF2FF surface area");const P n=normed(q.normal);const ComplexVec3 J=cross(n,q.magnetic);const ComplexVec3 M=scale(cross(n,q.electric),-1.0);const auto transverseJ=cross(r,cross(J,r));const auto term=add(scale(transverseJ,eta),cross(M,r));const double phase=k*(r.x*q.position.x+r.y*q.position.y+r.z*q.position.z);sum=add(sum,scale(term,q.area*std::exp(std::complex<double>(0,phase))));}sum=scale(sum,std::complex<double>(0,k/(4*std::numbers::pi)));return {sum,abs2(sum)/(2*eta)};}
void SarGrid3D::validate()const{const auto n=nx*ny*nz;if(nx==0||ny==0||nz==0||!(dx>0)||!(dy>0)||!(dz>0)||ex.size()!=n||ey.size()!=n||ez.size()!=n||conductivity.size()!=n||density.size()!=n)throw std::invalid_argument("invalid SAR grid");for(std::size_t i=0;i<n;++i)if(conductivity[i]<0||!(density[i]>0))throw std::invalid_argument("invalid SAR material");}
std::vector<double> local_sar(const SarGrid3D&g){g.validate();std::vector<double> out(g.ex.size());for(std::size_t i=0;i<out.size();++i)out[i]=g.conductivity[i]*(g.ex[i]*g.ex[i]+g.ey[i]*g.ey[i]+g.ez[i]*g.ez[i])/g.density[i];return out;}
std::vector<double> mass_averaged_sar(const SarGrid3D&g,double target){g.validate();if(!(target>0))throw std::invalid_argument("SAR target mass must be positive");auto sar=local_sar(g);std::vector<double> out(sar.size());const double vol=g.dx*g.dy*g.dz;for(std::size_t k=0;k<g.nz;++k)for(std::size_t j=0;j<g.ny;++j)for(std::size_t i=0;i<g.nx;++i){double mass=0,weighted=0;const std::size_t maxr=std::max({g.nx,g.ny,g.nz});for(std::size_t r=0;r<maxr&&mass<target;++r){const auto imin=i>r?i-r:0,imax=std::min(g.nx-1,i+r),jmin=j>r?j-r:0,jmax=std::min(g.ny-1,j+r),kmin=k>r?k-r:0,kmax=std::min(g.nz-1,k+r);mass=0;weighted=0;for(std::size_t kk=kmin;kk<=kmax;++kk)for(std::size_t jj=jmin;jj<=jmax;++jj)for(std::size_t ii=imin;ii<=imax;++ii){auto n=id(ii,jj,kk,g);const double m=g.density[n]*vol;mass+=m;weighted+=m*sar[n];}}out[id(i,j,k,g)]=weighted/std::max(mass,1e-30);}return out;}
std::vector<double> sar_1g(const SarGrid3D&g){return mass_averaged_sar(g,0.001);}std::vector<double> sar_10g(const SarGrid3D&g){return mass_averaged_sar(g,0.010);}
} // namespace cfd::fdtd
