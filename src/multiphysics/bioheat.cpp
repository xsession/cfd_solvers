#include "cfd/multiphysics/bioheat.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfd::multiphysics {

double sar_from_rms_electric_field(double conductivity,double density,double electric_rms){
    if(!(conductivity>=0.0)||!(density>0.0)||!(electric_rms>=0.0)||!std::isfinite(conductivity)||!std::isfinite(density)||!std::isfinite(electric_rms))throw std::invalid_argument("invalid SAR material/field data");
    return conductivity*electric_rms*electric_rms/density;
}


void VoxelTissueGrid3D::validate() const {
    if(nx==0U||ny==0U||nz==0U||!(dx_m>0.0)||!(dy_m>0.0)||!(dz_m>0.0)||materials.empty())
        throw std::invalid_argument("invalid voxel tissue grid geometry");
    const std::size_t n=nx*ny*nz;
    if(material_index.size()!=n||electric_rms_v_per_m.size()!=n)throw std::invalid_argument("voxel tissue field size mismatch");
    for(const auto& material:materials){
        if(!(material.density_kg_per_m3>0.0)||!(material.specific_heat_j_per_kg_k>0.0)||!(material.thermal_conductivity_w_per_m_k>=0.0)
           ||!(material.electrical_conductivity_s_per_m>=0.0)||!(material.blood_perfusion_per_s>=0.0))
            throw std::invalid_argument("invalid voxel tissue material");
    }
    for(std::size_t i=0;i<n;++i){
        if(material_index[i]>=materials.size()||!(electric_rms_v_per_m[i]>=0.0)||!std::isfinite(electric_rms_v_per_m[i]))
            throw std::invalid_argument("invalid voxel tissue index/field");
    }
}

std::vector<double> voxel_sar_w_per_kg(const VoxelTissueGrid3D& grid){
    grid.validate();const std::size_t n=grid.nx*grid.ny*grid.nz;std::vector<double> sar(n);
    for(std::size_t i=0;i<n;++i){const auto& material=grid.materials[grid.material_index[i]];sar[i]=sar_from_rms_electric_field(material.electrical_conductivity_s_per_m,material.density_kg_per_m3,grid.electric_rms_v_per_m[i]);}
    return sar;
}

double max_mass_averaged_sar_w_per_kg(const VoxelTissueGrid3D& grid,std::span<const double> sar,double averaging_mass_kg){
    grid.validate();const std::size_t n=grid.nx*grid.ny*grid.nz;if(sar.size()!=n||!(averaging_mass_kg>0.0))throw std::invalid_argument("invalid SAR averaging request");
    const double volume=grid.dx_m*grid.dy_m*grid.dz_m;double best=0.0;
    for(std::size_t center=0;center<n;++center){
        const std::size_t cz=center/(grid.nx*grid.ny),rem=center%(grid.nx*grid.ny),cy=rem/grid.nx,cx=rem%grid.nx;
        std::vector<std::pair<double,std::size_t>> order;order.reserve(n);
        for(std::size_t idx=0;idx<n;++idx){const std::size_t z=idx/(grid.nx*grid.ny),r=idx%(grid.nx*grid.ny),y=r/grid.nx,x=r%grid.nx;const double dx=(static_cast<double>(x)-static_cast<double>(cx))*grid.dx_m,dy=(static_cast<double>(y)-static_cast<double>(cy))*grid.dy_m,dz=(static_cast<double>(z)-static_cast<double>(cz))*grid.dz_m;order.push_back({dx*dx+dy*dy+dz*dz,idx});}
        std::sort(order.begin(),order.end(),[](const auto& a,const auto& b){return a.first<b.first;});
        double mass=0.0,weighted=0.0;
        for(const auto& item:order){const std::size_t idx=item.second;if(!(sar[idx]>=0.0)||!std::isfinite(sar[idx]))throw std::invalid_argument("invalid SAR value");const double voxel_mass=grid.materials[grid.material_index[idx]].density_kg_per_m3*volume;const double take=std::min(voxel_mass,averaging_mass_kg-mass);weighted+=take*sar[idx];mass+=take;if(mass>=averaging_mass_kg-1.0e-18)break;}
        if(mass>0.0)best=std::max(best,weighted/mass);
    }
    return best;
}

std::vector<double> project_voxel_sar_to_pennes2d(const VoxelTissueGrid3D& grid,std::span<const double> sar,std::size_t nx,std::size_t ny){
    grid.validate();const std::size_t n=grid.nx*grid.ny*grid.nz;if(sar.size()!=n||nx==0U||ny==0U)throw std::invalid_argument("invalid voxel-to-Pennes projection");
    std::vector<double> projected(nx*ny,0.0),mass(nx*ny,0.0);const double volume=grid.dx_m*grid.dy_m*grid.dz_m;
    for(std::size_t z=0;z<grid.nz;++z)for(std::size_t y=0;y<grid.ny;++y)for(std::size_t x=0;x<grid.nx;++x){
        const std::size_t src=(z*grid.ny+y)*grid.nx+x;const std::size_t tx=std::min(nx-1U,x*nx/grid.nx),ty=std::min(ny-1U,y*ny/grid.ny),dst=ty*nx+tx;
        const double m=grid.materials[grid.material_index[src]].density_kg_per_m3*volume;projected[dst]+=m*sar[src];mass[dst]+=m;
    }
    for(std::size_t i=0;i<projected.size();++i)if(mass[i]>0.0)projected[i]/=mass[i];
    return projected;
}


PennesBioheat2D::PennesBioheat2D(PennesBioheat2DConfig config):config_(config){
    if(config_.nx<3U||config_.ny<3U||!(config_.width_m>0.0)||!(config_.height_m>0.0)
       ||!(config_.tissue_density_kg_per_m3>0.0)||!(config_.tissue_specific_heat_j_per_kg_k>0.0)
       ||!(config_.thermal_conductivity_w_per_m_k>=0.0)||!(config_.blood_density_kg_per_m3>0.0)
       ||!(config_.blood_specific_heat_j_per_kg_k>0.0)||!(config_.blood_perfusion_per_s>=0.0)
       ||!(config_.dt_s>0.0)||config_.max_iterations==0U||!(config_.relative_tolerance>0.0))
        throw std::invalid_argument("invalid Pennes bioheat configuration");
    temperature_.assign(config_.nx*config_.ny,config_.blood_temperature_k);sar_.assign(temperature_.size(),0.0);volumetric_heating_.assign(temperature_.size(),0.0);
}
void PennesBioheat2D::initialize(double temperature){if(!std::isfinite(temperature))throw std::invalid_argument("invalid initial tissue temperature");std::fill(temperature_.begin(),temperature_.end(),temperature);time_s_=0.0;}
void PennesBioheat2D::set_sar(double value){if(!(value>=0.0)||!std::isfinite(value))throw std::invalid_argument("invalid SAR");std::fill(sar_.begin(),sar_.end(),value);}
void PennesBioheat2D::set_sar(std::span<const double> values){if(values.size()!=sar_.size())throw std::invalid_argument("SAR map size mismatch");for(double value:values)if(!(value>=0.0)||!std::isfinite(value))throw std::invalid_argument("invalid SAR map value");std::copy(values.begin(),values.end(),sar_.begin());}
void PennesBioheat2D::set_volumetric_heating(double value){if(!(value>=0.0)||!std::isfinite(value))throw std::invalid_argument("invalid volumetric heating");std::fill(volumetric_heating_.begin(),volumetric_heating_.end(),value);}
void PennesBioheat2D::set_volumetric_heating(std::span<const double> values){if(values.size()!=volumetric_heating_.size())throw std::invalid_argument("volumetric heating map size mismatch");for(double value:values)if(!(value>=0.0)||!std::isfinite(value))throw std::invalid_argument("invalid volumetric heating map value");std::copy(values.begin(),values.end(),volumetric_heating_.begin());}
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
    for(std::size_t y=0;y<config_.ny;++y)for(std::size_t x=0;x<config_.nx;++x){const std::size_t i=index(x,y);if(!periodic&&boundary_node(x,y)){rhs[i]=config_.fixed_boundary_temperature_k;next[i]=rhs[i];continue;}rhs[i]=transient*temperature_[i]+perfusion*config_.blood_temperature_k+config_.metabolic_heat_w_per_m3+config_.tissue_density_kg_per_m3*sar_[i]+volumetric_heating_[i];if(!periodic){if(x==1U)rhs[i]+=kx*config_.fixed_boundary_temperature_k;if(x+2U==config_.nx)rhs[i]+=kx*config_.fixed_boundary_temperature_k;if(y==1U)rhs[i]+=ky*config_.fixed_boundary_temperature_k;if(y+2U==config_.ny)rhs[i]+=ky*config_.fixed_boundary_temperature_k;}}
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
