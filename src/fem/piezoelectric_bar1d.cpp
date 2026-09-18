#include "cfd/solvers/fem/piezoelectric_bar1d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cfd::fem {
namespace {
void dense_solve(std::vector<double>& a,std::vector<double>& b,std::size_t n){
    for(std::size_t k=0;k<n;++k){
        std::size_t pivot=k;double best=std::abs(a[k*n+k]);
        for(std::size_t r=k+1;r<n;++r){const double v=std::abs(a[r*n+k]);if(v>best){best=v;pivot=r;}}
        if(!(best>1e-18))throw std::runtime_error("singular piezoelectric FEM matrix");
        if(pivot!=k){for(std::size_t c=k;c<n;++c)std::swap(a[k*n+c],a[pivot*n+c]);std::swap(b[k],b[pivot]);}
        const double d=a[k*n+k];
        for(std::size_t r=k+1;r<n;++r){const double f=a[r*n+k]/d;if(f==0.0)continue;a[r*n+k]=0.0;for(std::size_t c=k+1;c<n;++c)a[r*n+c]-=f*a[k*n+c];b[r]-=f*b[k];}
    }
    for(std::size_t ii=0;ii<n;++ii){const std::size_t i=n-1U-ii;double s=b[i];for(std::size_t c=i+1;c<n;++c)s-=a[i*n+c]*b[c];b[i]=s/a[i*n+i];}
}
void prescribe(std::vector<double>& a,std::vector<double>& b,std::size_t n,std::size_t dof,double value){
    for(std::size_t r=0;r<n;++r){if(r==dof)continue;b[r]-=a[r*n+dof]*value;a[r*n+dof]=0.0;}
    for(std::size_t c=0;c<n;++c) a[dof*n+c]=0.0;
    a[dof*n+dof]=1.0;
    b[dof]=value;
}
}

PiezoelectricBar1D::PiezoelectricBar1D(PiezoelectricBar1DConfig c):config_(c){
    if(c.elements==0U||!(c.length_m>0.0)||!(c.area_m2>0.0)||!(c.elastic_modulus_pa>0.0)||
       !std::isfinite(c.piezoelectric_stress_c_per_m2)||!(c.permittivity_f_per_m>0.0))
        throw std::invalid_argument("invalid piezoelectric bar configuration");
    displacement_.assign(c.elements+1U,0.0);potential_.assign(c.elements+1U,0.0);
}

void PiezoelectricBar1D::solve_receive(double tip_force_n){solve(tip_force_n,false,0.0);}
void PiezoelectricBar1D::solve_transmit(double voltage,double tip_force_n){if(!std::isfinite(voltage))throw std::invalid_argument("invalid piezo voltage");solve(tip_force_n,true,voltage);}

void PiezoelectricBar1D::solve(double tip_force_n,bool prescribe_tip_voltage,double tip_voltage){
    if(!std::isfinite(tip_force_n))throw std::invalid_argument("invalid piezo tip force");
    const std::size_t nodes=config_.elements+1U, ndof=2U*nodes;
    std::vector<double> a(ndof*ndof,0.0),rhs(ndof,0.0);
    const double h=config_.length_m/static_cast<double>(config_.elements);
    const double kuu=config_.elastic_modulus_pa*config_.area_m2/h;
    const double kup=config_.piezoelectric_stress_c_per_m2*config_.area_m2/h;
    const double kpp=-config_.permittivity_f_per_m*config_.area_m2/h;
    const double pattern[2][2]={{1.0,-1.0},{-1.0,1.0}};
    for(std::size_t e=0;e<config_.elements;++e){
        const std::size_t n0=e,n1=e+1U;const std::size_t ns[2]={n0,n1};
        for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j){
            const double s=pattern[i][j];const auto ui=ns[i],uj=ns[j],pi=nodes+ns[i],pj=nodes+ns[j];
            a[ui*ndof+uj]+=kuu*s;
            a[ui*ndof+pj]+=kup*s;
            a[pi*ndof+uj]+=kup*s;
            a[pi*ndof+pj]+=kpp*s;
        }
    }
    rhs[nodes-1U]=tip_force_n;
    prescribe(a,rhs,ndof,0U,0.0);
    prescribe(a,rhs,ndof,nodes,0.0);
    if(prescribe_tip_voltage)prescribe(a,rhs,ndof,2U*nodes-1U,tip_voltage);
    dense_solve(a,rhs,ndof);
    for(std::size_t i=0;i<nodes;++i){displacement_[i]=rhs[i];potential_[i]=rhs[nodes+i];}
    // D = e*S - eps*dphi/dx, Q=A*D at the tip element.
    const double strain=(displacement_[nodes-1U]-displacement_[nodes-2U])/h;
    const double efield=-(potential_[nodes-1U]-potential_[nodes-2U])/h;
    const double d=config_.piezoelectric_stress_c_per_m2*strain+config_.permittivity_f_per_m*efield;
    electrode_charge_c_=config_.area_m2*d;
}

} // namespace cfd::fem
