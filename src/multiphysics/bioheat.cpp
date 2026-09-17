#include "cfd/multiphysics/bioheat.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::multiphysics {

double sar_from_rms_electric_field(double conductivity,double density,double electric_rms){
    if(!(conductivity>=0.0)||!(density>0.0)||!(electric_rms>=0.0)||!std::isfinite(conductivity)||!std::isfinite(density)||!std::isfinite(electric_rms))throw std::invalid_argument("invalid SAR material/field data");
    return conductivity*electric_rms*electric_rms/density;
}

PennesBioheat2D::PennesBioheat2D(PennesBioheat2DConfig config):config_(config){
    if(config_.nx<3U||config_.ny<3U||!(config_.width_m>0.0)||!(config_.height_m>0.0)
       ||!(config_.tissue_density_kg_per_m3>0.0)||!(config_.tissue_specific_heat_j_per_kg_k>0.0)
       ||!(config_.thermal_conductivity_w_per_m_k>=0.0)||!(config_.blood_density_kg_per_m3>0.0)
       ||!(config_.blood_specific_heat_j_per_kg_k>0.0)||!(config_.blood_perfusion_per_s>=0.0)
       ||!(config_.dt_s>0.0)||config_.max_iterations==0U||!(config_.relative_tolerance>0.0))
        throw std::invalid_argument("invalid Pennes bioheat configuration");
    temperature_.assign(config_.nx*config_.ny,config_.blood_temperature_k);sar_.assign(temperature_.size(),0.0);
}
void PennesBioheat2D::initialize(double temperature){if(!std::isfinite(temperature))throw std::invalid_argument("invalid initial tissue temperature");std::fill(temperature_.begin(),temperature_.end(),temperature);time_s_=0.0;}
void PennesBioheat2D::set_sar(double value){if(!(value>=0.0)||!std::isfinite(value))throw std::invalid_argument("invalid SAR");std::fill(sar_.begin(),sar_.end(),value);}
void PennesBioheat2D::set_sar(std::span<const double> values){if(values.size()!=sar_.size())throw std::invalid_argument("SAR map size mismatch");for(double value:values)if(!(value>=0.0)||!std::isfinite(value))throw std::invalid_argument("invalid SAR map value");std::copy(values.begin(),values.end(),sar_.begin());}
bool PennesBioheat2D::boundary_node(std::size_t x,std::size_t y) const noexcept{return x==0U||y==0U||x+1U==config_.nx||y+1U==config_.ny;}

void PennesBioheat2D::step_once(){
    const std::size_t n=temperature_.size();const bool periodic=config_.boundary==BioheatBoundary::periodic;
    const double dx=config_.width_m/static_cast<double>(periodic?config_.nx:config_.nx-1U),dy=config_.height_m/static_cast<double>(periodic?config_.ny:config_.ny-1U);
    const double kx=config_.thermal_conductivity_w_per_m_k/(dx*dx),ky=config_.thermal_conductivity_w_per_m_k/(dy*dy);
    const double capacity=config_.tissue_density_kg_per_m3*config_.tissue_specific_heat_j_per_kg_k;
    const double transient=capacity/config_.dt_s;
    const double perfusion=config_.blood_density_kg_per_m3*config_.blood_specific_heat_j_per_kg_k*config_.blood_perfusion_per_s;
    const double diagonal=transient+perfusion+2.0*(kx+ky);
    std::vector<double> rhs(n),next=temperature_;
    for(std::size_t y=0;y<config_.ny;++y)for(std::size_t x=0;x<config_.nx;++x){const std::size_t i=index(x,y);if(!periodic&&boundary_node(x,y)){rhs[i]=config_.fixed_boundary_temperature_k;next[i]=rhs[i];continue;}rhs[i]=transient*temperature_[i]+perfusion*config_.blood_temperature_k+config_.metabolic_heat_w_per_m3+config_.tissue_density_kg_per_m3*sar_[i];if(!periodic){if(x==1U)rhs[i]+=kx*config_.fixed_boundary_temperature_k;if(x+2U==config_.nx)rhs[i]+=kx*config_.fixed_boundary_temperature_k;if(y==1U)rhs[i]+=ky*config_.fixed_boundary_temperature_k;if(y+2U==config_.ny)rhs[i]+=ky*config_.fixed_boundary_temperature_k;}}
    const auto apply=[&](std::span<const double> in,std::span<double> out){
        for(std::size_t y=0;y<config_.ny;++y)for(std::size_t x=0;x<config_.nx;++x){const std::size_t i=index(x,y);if(!periodic&&boundary_node(x,y)){out[i]=in[i];continue;}double value=diagonal*in[i];
            if(periodic){const std::size_t xl=(x+config_.nx-1U)%config_.nx,xr=(x+1U)%config_.nx,yd=(y+config_.ny-1U)%config_.ny,yu=(y+1U)%config_.ny;value-=kx*(in[index(xl,y)]+in[index(xr,y)])+ky*(in[index(x,yd)]+in[index(x,yu)]);}
            else{if(x>1U)value-=kx*in[index(x-1U,y)];if(x+2U<config_.nx)value-=kx*in[index(x+1U,y)];if(y>1U)value-=ky*in[index(x,y-1U)];if(y+2U<config_.ny)value-=ky*in[index(x,y+1U)];}
            out[i]=value;
        }
    };
    linear_result_=cfd::core::conjugate_gradient(rhs,next,apply,workspace_,config_.max_iterations,config_.relative_tolerance);
    if(!linear_result_.converged)throw std::runtime_error("Pennes bioheat linear solve failed");
    temperature_=std::move(next);
    time_s_+=config_.dt_s;
}
void PennesBioheat2D::step(std::size_t steps){for(std::size_t i=0;i<steps;++i)step_once();}

} // namespace cfd::multiphysics
