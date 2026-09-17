#include "cfd/optics/wavefront.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd::optics {namespace {
double factorial(int n){double v=1;for(int i=2;i<=n;++i)v*=i;return v;}
double radial(int n,int m,double r){m=std::abs(m);if(n<0||m>n||((n-m)&1))return 0.0;double sum=0;for(int k=0;k<=(n-m)/2;++k){double c=std::pow(-1.0,k)*factorial(n-k)/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k));sum+=c*std::pow(r,n-2*k);}return sum;}
std::vector<std::pair<int,int>> modes(unsigned maxn){std::vector<std::pair<int,int>> out;for(unsigned n=0;n<=maxn;++n)for(int m=-static_cast<int>(n);m<=static_cast<int>(n);m+=2)out.emplace_back(static_cast<int>(n),m);return out;}
std::vector<double> solve_dense(std::vector<double>a,std::vector<double>b,std::size_t n){for(std::size_t c=0;c<n;++c){std::size_t p=c;for(std::size_t r=c+1;r<n;++r)if(std::abs(a[r*n+c])>std::abs(a[p*n+c]))p=r;if(std::abs(a[p*n+c])<1e-14)throw std::runtime_error("singular Zernike fit");if(p!=c){for(std::size_t j=c;j<n;++j)std::swap(a[p*n+j],a[c*n+j]);std::swap(b[p],b[c]);}double d=a[c*n+c];for(std::size_t j=c;j<n;++j)a[c*n+j]/=d;b[c]/=d;for(std::size_t r=0;r<n;++r)if(r!=c){double f=a[r*n+c];for(std::size_t j=c;j<n;++j)a[r*n+j]-=f*a[c*n+j];b[r]-=f*b[c];}}return b;}
}
double zernike(int n,int m,double rho,double theta){if(rho<0||rho>1||!std::isfinite(rho)||!std::isfinite(theta))throw std::invalid_argument("invalid Zernike coordinate");double R=radial(n,m,rho);if(m>0)return R*std::cos(m*theta);if(m<0)return R*std::sin(-m*theta);return R;}
std::vector<ZernikeCoefficient> fit_zernike(std::span<const WavefrontSample>s,unsigned maxn){auto md=modes(maxn);if(s.size()<md.size())throw std::invalid_argument("insufficient wavefront samples");std::size_t n=md.size();std::vector<double>A(n*n),b(n);for(auto q:s){std::vector<double>v(n);for(std::size_t i=0;i<n;++i)v[i]=zernike(md[i].first,md[i].second,q.rho,q.theta);for(std::size_t i=0;i<n;++i){b[i]+=v[i]*q.opd;for(std::size_t j=0;j<n;++j)A[i*n+j]+=v[i]*v[j];}}auto x=solve_dense(std::move(A),std::move(b),n);std::vector<ZernikeCoefficient>out(n);for(std::size_t i=0;i<n;++i)out[i]={md[i].first,md[i].second,x[i]};return out;}
double wavefront_rms(std::span<const WavefrontSample>s,bool remove){if(s.empty())throw std::invalid_argument("empty wavefront");double mean=0;for(auto q:s)mean+=q.opd;mean/=static_cast<double>(s.size());if(!remove)mean=0;double v=0;for(auto q:s){double d=q.opd-mean;v+=d*d;}return std::sqrt(v/static_cast<double>(s.size()));}
double longitudinal_chromatic_shift(std::span<const double>f){if(f.empty())throw std::invalid_argument("empty chromatic focal set");auto [lo,hi]=std::minmax_element(f.begin(),f.end());return *hi-*lo;}
} // namespace cfd::optics
