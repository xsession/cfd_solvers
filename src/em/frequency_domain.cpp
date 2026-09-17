#include "cfd/em/frequency_domain.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace cfd::em {
namespace {
constexpr double epsilon0=8.8541878128e-12;
constexpr double mu0=1.25663706212e-6;

std::vector<Complex> solve_dense_complex(std::vector<Complex> matrix,std::vector<Complex> rhs,std::size_t n) {
    if(matrix.size()!=n*n||rhs.size()!=n)throw std::invalid_argument("complex dense solve size mismatch");
    for(std::size_t k=0;k<n;++k){
        std::size_t pivot=k;double best=std::abs(matrix[k*n+k]);
        for(std::size_t row=k+1;row<n;++row){const double value=std::abs(matrix[row*n+k]);if(value>best){best=value;pivot=row;}}
        if(!(best>std::numeric_limits<double>::epsilon()))throw std::runtime_error("singular frequency-domain Maxwell matrix");
        if(pivot!=k){for(std::size_t col=k;col<n;++col)std::swap(matrix[k*n+col],matrix[pivot*n+col]);std::swap(rhs[k],rhs[pivot]);}
        const Complex diag=matrix[k*n+k];
        for(std::size_t row=k+1;row<n;++row){
            const Complex factor=matrix[row*n+k]/diag;if(std::abs(factor)==0.0)continue;
            matrix[row*n+k]={};
            for(std::size_t col=k+1;col<n;++col)matrix[row*n+col]-=factor*matrix[k*n+col];
            rhs[row]-=factor*rhs[k];
        }
    }
    std::vector<Complex> x(n);
    for(std::size_t back=0;back<n;++back){const std::size_t row=n-1U-back;Complex sum=rhs[row];for(std::size_t col=row+1;col<n;++col)sum-=matrix[row*n+col]*x[col];x[row]=sum/matrix[row*n+row];}
    return x;
}

void jacobi_symmetric(std::vector<double>& a,std::vector<double>& vectors,std::size_t n) {
    vectors.assign(n*n,0.0);for(std::size_t i=0;i<n;++i)vectors[i*n+i]=1.0;
    const std::size_t max_iterations=80U*n*n;
    for(std::size_t iteration=0;iteration<max_iterations;++iteration){
        std::size_t p=0,q=0;double largest=0.0;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j){const double value=std::abs(a[i*n+j]);if(value>largest){largest=value;p=i;q=j;}}
        double diagonal_scale=1.0;for(std::size_t i=0;i<n;++i)diagonal_scale=std::max(diagonal_scale,std::abs(a[i*n+i]));
        if(largest<=1.0e-13*diagonal_scale)return;
        const double app=a[p*n+p],aqq=a[q*n+q],apq=a[p*n+q];
        const double angle=0.5*std::atan2(2.0*apq,aqq-app);const double c=std::cos(angle),s=std::sin(angle);
        for(std::size_t k=0;k<n;++k){if(k==p||k==q)continue;const double akp=a[k*n+p],akq=a[k*n+q];const double new_kp=c*akp-s*akq,new_kq=s*akp+c*akq;a[k*n+p]=a[p*n+k]=new_kp;a[k*n+q]=a[q*n+k]=new_kq;}
        a[p*n+p]=c*c*app-2.0*s*c*apq+s*s*aqq;
        a[q*n+q]=s*s*app+2.0*s*c*apq+c*c*aqq;
        a[p*n+q]=a[q*n+p]=0.0;
        for(std::size_t k=0;k<n;++k){const double vkp=vectors[k*n+p],vkq=vectors[k*n+q];vectors[k*n+p]=c*vkp-s*vkq;vectors[k*n+q]=s*vkp+c*vkq;}
    }
    throw std::runtime_error("cavity Jacobi eigensolver failed to converge");
}
}

FrequencyDomain1DResult solve_pec_driven_maxwell_1d(const FrequencyDomain1DConfig& config,
                                                      std::span<const Complex> impressed_current_density) {
    if(config.points<3U||!(config.length_m>0.0)||!(config.frequency_hz>0.0)
       ||!(config.relative_permittivity>0.0)||!(config.relative_permeability>0.0)
       ||!(config.conductivity_s_per_m>=0.0)||impressed_current_density.size()!=config.points)
        throw std::invalid_argument("invalid driven 1-D Maxwell configuration");
    const std::size_t interior=config.points-2U;const double dx=config.length_m/static_cast<double>(config.points-1U);
    const double omega=2.0*std::numbers::pi*config.frequency_hz;
    const double mu=mu0*config.relative_permeability,epsilon=epsilon0*config.relative_permittivity;
    const Complex k2{omega*omega*mu*epsilon,-omega*mu*config.conductivity_s_per_m};
    const double inv_dx2=1.0/(dx*dx);std::vector<Complex> matrix(interior*interior),rhs(interior);
    for(std::size_t row=0;row<interior;++row){
        matrix[row*interior+row]=Complex{-2.0*inv_dx2,0.0}+k2;
        if(row>0U)matrix[row*interior+row-1U]=inv_dx2;
        if(row+1U<interior)matrix[row*interior+row+1U]=inv_dx2;
        rhs[row]=Complex{0.0,omega*mu}*impressed_current_density[row+1U];
    }
    const auto interior_e=solve_dense_complex(std::move(matrix),std::move(rhs),interior);
    FrequencyDomain1DResult result;result.spacing_m=dx;result.electric_v_per_m.assign(config.points,{});result.magnetic_a_per_m.resize(config.points-1U);
    for(std::size_t i=0;i<interior;++i)result.electric_v_per_m[i+1U]=interior_e[i];
    const Complex h_factor{0.0,1.0/(omega*mu*dx)};
    for(std::size_t i=0;i+1U<config.points;++i)result.magnetic_a_per_m[i]=h_factor*(result.electric_v_per_m[i+1U]-result.electric_v_per_m[i]);
    return result;
}

std::vector<CavityMode1D> pec_cavity_eigenmodes_1d(std::size_t points,double length_m,
                                                    std::span<const double> relative_permittivity,
                                                    double relative_permeability,std::size_t mode_count) {
    if(points<3U||!(length_m>0.0)||relative_permittivity.size()!=points||!(relative_permeability>0.0)||mode_count==0U)
        throw std::invalid_argument("invalid PEC cavity eigenmode configuration");
    for(double value:relative_permittivity)if(!(value>0.0)||!std::isfinite(value))throw std::invalid_argument("invalid cavity permittivity");
    const std::size_t n=points-2U;if(mode_count>n)throw std::invalid_argument("requested too many cavity modes");
    const double dx=length_m/static_cast<double>(points-1U),inv_dx2=1.0/(dx*dx),mu=mu0*relative_permeability;
    std::vector<double> mass(n);for(std::size_t i=0;i<n;++i)mass[i]=mu*epsilon0*relative_permittivity[i+1U];
    std::vector<double> matrix(n*n,0.0);
    for(std::size_t row=0;row<n;++row){
        matrix[row*n+row]=2.0*inv_dx2/mass[row];
        if(row+1U<n){const double value=-inv_dx2/std::sqrt(mass[row]*mass[row+1U]);matrix[row*n+row+1U]=matrix[(row+1U)*n+row]=value;}
    }
    std::vector<double> eigenvectors;jacobi_symmetric(matrix,eigenvectors,n);
    std::vector<std::size_t> order(n);std::iota(order.begin(),order.end(),0U);
    std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return matrix[a*n+a]<matrix[b*n+b];});
    std::vector<CavityMode1D> modes;modes.reserve(mode_count);
    for(std::size_t m=0;m<mode_count;++m){
        const std::size_t column=order[m];const double lambda=matrix[column*n+column];if(!(lambda>0.0))throw std::runtime_error("non-positive Maxwell cavity eigenvalue");
        CavityMode1D mode;mode.angular_frequency_rad_per_s=std::sqrt(lambda);mode.frequency_hz=mode.angular_frequency_rad_per_s/(2.0*std::numbers::pi);mode.electric_shape.assign(points,0.0);
        double peak=0.0;for(std::size_t i=0;i<n;++i){mode.electric_shape[i+1U]=eigenvectors[i*n+column]/std::sqrt(mass[i]);peak=std::max(peak,std::abs(mode.electric_shape[i+1U]));}
        if(!(peak>0.0))throw std::runtime_error("empty cavity eigenvector");
        for(double& value:mode.electric_shape)value/=peak;
        modes.push_back(std::move(mode));
    }
    return modes;
}

std::vector<CavityMode1D> pec_cavity_eigenmodes_1d(std::size_t points,double length_m,double relative_permittivity,
                                                    double relative_permeability,std::size_t mode_count) {
    std::vector<double> epsilon(points,relative_permittivity);return pec_cavity_eigenmodes_1d(points,length_m,epsilon,relative_permeability,mode_count);
}

} // namespace cfd::em
