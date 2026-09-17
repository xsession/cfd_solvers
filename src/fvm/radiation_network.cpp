#include "cfd/fvm/radiation_network.hpp"
#include "cfd/fvm/radiation.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
RadiosityResult solve_gray_radiosity(std::span<const double>T,std::span<const double>e,std::span<const double>A,std::span<const double>F){const std::size_t n=T.size();if(n==0||e.size()!=n||A.size()!=n||F.size()!=n*n)throw std::invalid_argument("radiosity size mismatch");std::vector<double>M(n*n),b(n);for(std::size_t i=0;i<n;++i){if(T[i]<0||!(e[i]>0&&e[i]<=1)||!(A[i]>0))throw std::invalid_argument("invalid radiosity surface");double row=0;for(std::size_t j=0;j<n;++j){if(F[i*n+j]<0)throw std::invalid_argument("negative view factor");row+=F[i*n+j];M[i*n+j]=(i==j?1.0:0.0)-(1.0-e[i])*F[i*n+j];}if(std::abs(row-1.0)>1e-8)throw std::invalid_argument("view-factor row must sum to one");b[i]=e[i]*stefan_boltzmann*std::pow(T[i],4);}for(std::size_t k=0;k<n;++k){std::size_t p=k;for(std::size_t r=k+1;r<n;++r)if(std::abs(M[r*n+k])>std::abs(M[p*n+k]))p=r;if(std::abs(M[p*n+k])<1e-30)throw std::runtime_error("singular radiosity system");if(p!=k){for(std::size_t j=k;j<n;++j)std::swap(M[p*n+j],M[k*n+j]);std::swap(b[p],b[k]);}double d=M[k*n+k];for(std::size_t j=k;j<n;++j)M[k*n+j]/=d;b[k]/=d;for(std::size_t r=0;r<n;++r)if(r!=k){double q=M[r*n+k];for(std::size_t j=k;j<n;++j)M[r*n+j]-=q*M[k*n+j];b[r]-=q*b[k];}}
RadiosityResult out;out.radiosity=b;out.net_heat.assign(n,0);for(std::size_t i=0;i<n;++i){double G=0;for(std::size_t j=0;j<n;++j)G+=F[i*n+j]*b[j];out.net_heat[i]=A[i]*(b[i]-G);}return out;}
double OpticallyThinRadiation::volumetric_source(double t,double env)const{return optically_thin_radiation_source(t,env,absorption_);}
} // namespace cfd::fvm
