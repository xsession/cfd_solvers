#include "cfd/solvers/fdtd/maxwell1d.hpp"

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fdtd {
namespace {
constexpr double epsilon0 = 8.8541878128e-12;
constexpr double mu0 = 1.25663706212e-6;
Maxwell1DConfig validated(Maxwell1DConfig c) {
    if(c.cells<4U || !(c.dx>0.0) || !std::isfinite(c.dx)
       || !(c.courant>0.0&&c.courant<1.0) || !(c.epsilon_r>0.0)
       || !std::isfinite(c.epsilon_r) || !(c.mu_r>0.0) || !std::isfinite(c.mu_r))
        throw std::invalid_argument("invalid Maxwell1D configuration");
    if(c.boundary==Boundary1D::pml && (c.pml_cells<2U || c.pml_cells>(c.cells-4U)/2U
       || !(c.pml_order>=1.0) || !std::isfinite(c.pml_order)
       || !(c.pml_target_reflection>0.0&&c.pml_target_reflection<1.0)))
        throw std::invalid_argument("invalid Maxwell1D PML configuration");
    return c;
}
}

Maxwell1D::Maxwell1D(Maxwell1DConfig config)
    : config_(validated(config)),mu_(mu0*config.mu_r),epsilon_(config.cells,epsilon0*config.epsilon_r),
      conductivity_(config.cells,0.0),e_decay_(config.cells,1.0),e_curl_(config.cells,0.0),
      pml_sigma_e_(config.cells,0.0),h_decay_(config.cells-1U,1.0),h_curl_(config.cells-1U,0.0),
      ez_(config.cells,0.0),hy_(config.cells-1U,0.0),
      dispersion_(config.cells,DispersionModel1D::none),polarization_(config.cells,0.0),
      polarization_previous_(config.cells,0.0),dispersion_a_(config.cells,0.0),
      dispersion_b_(config.cells,0.0),dispersion_c_(config.cells,0.0) {
    if(config.cells<4U||!(config.dx>0.0)||!(config.courant>0.0&&config.courant<1.0)||
       !(config.epsilon_r>0.0)||!(config.mu_r>0.0)) throw std::invalid_argument("invalid Maxwell1D configuration");
    if(config.boundary==Boundary1D::pml&&
       (config.pml_cells<2U||2U*config.pml_cells+4U>config.cells||!(config.pml_order>=1.0)||
        !(config.pml_target_reflection>0.0&&config.pml_target_reflection<1.0)))
        throw std::invalid_argument("invalid Maxwell1D PML configuration");
    const double wave_speed=1.0/std::sqrt(epsilon_.front()*mu_);
    dt_=config.courant*config.dx/wave_speed;
    initialize_pml_profile();
    preprocess_material_coefficients();
}

void Maxwell1D::initialize_pml_profile(){
    std::fill(pml_sigma_e_.begin(),pml_sigma_e_.end(),0.0);
    if(config_.boundary!=Boundary1D::pml) return;
    const double eta=std::sqrt(mu_/epsilon_.front());
    const double thickness=static_cast<double>(config_.pml_cells)*config_.dx;
    const double sigma_max=-(config_.pml_order+1.0)*std::log(config_.pml_target_reflection)
        /(2.0*eta*thickness);
    for(std::size_t i=0;i<config_.cells;++i){
        double depth=0.0;
        if(i<config_.pml_cells){
            depth=(static_cast<double>(config_.pml_cells)-static_cast<double>(i)-0.5)
                /static_cast<double>(config_.pml_cells);
        } else if(i+config_.pml_cells>=config_.cells){
            depth=(static_cast<double>(i)-static_cast<double>(config_.cells-config_.pml_cells)+0.5)
                /static_cast<double>(config_.pml_cells);
        }
        if(depth>0.0) pml_sigma_e_[i]=sigma_max*std::pow(std::min(depth,1.0),config_.pml_order);
    }
}

void Maxwell1D::preprocess_material_coefficients(){
    for(std::size_t i=0;i<config_.cells;++i){
        const double loss=(conductivity_[i]+pml_sigma_e_[i])*dt_/(2.0*epsilon_[i]);
        const double denom=1.0+loss;
        e_decay_[i]=(1.0-loss)/denom;
        e_curl_[i]=(dt_/(epsilon_[i]*config_.dx))/denom;
    }
    for(std::size_t i=0;i<hy_.size();++i){
        const double eps_edge=0.5*(epsilon_[i]+epsilon_[i+1U]);
        const double sigma_e_edge=0.5*(pml_sigma_e_[i]+pml_sigma_e_[i+1U]);
        const double sigma_m=sigma_e_edge*mu_/eps_edge;
        const double loss=sigma_m*dt_/(2.0*mu_);
        const double denom=1.0+loss;
        h_decay_[i]=(1.0-loss)/denom;
        h_curl_[i]=(dt_/(mu_*config_.dx))/denom;
    }
}

void Maxwell1D::clear_dispersion(std::size_t begin_cell,std::size_t end_cell){
    for(std::size_t i=begin_cell;i<end_cell;++i){
        dispersion_[i]=DispersionModel1D::none;
        polarization_[i]=0.0;
        polarization_previous_[i]=0.0;
        dispersion_a_[i]=0.0;
        dispersion_b_[i]=0.0;
        dispersion_c_[i]=0.0;
    }
}

void Maxwell1D::set_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_r,double conductivity){
    validate_material(begin_cell,end_cell,epsilon_r,conductivity,false);
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_r>0.0)||!(conductivity>=0.0)) throw std::invalid_argument("invalid Maxwell1D material range");
    clear_dispersion(begin_cell,end_cell);
    for(std::size_t i=begin_cell;i<end_cell;++i){epsilon_[i]=epsilon0*epsilon_r;conductivity_[i]=conductivity;}
    preprocess_material_coefficients();
}

void Maxwell1D::validate_material(std::size_t begin,std::size_t end,double epsilon_r,
                                  double conductivity,bool dispersive) const {
    if(begin>=end || end>config_.cells || !(epsilon_r>0.0) || !std::isfinite(epsilon_r)
       || !(conductivity>=0.0) || !std::isfinite(conductivity))
        throw std::invalid_argument("invalid Maxwell1D material");
    const double courant=dt_/(config_.dx*std::sqrt(epsilon0*epsilon_r*mu_));
    if(!(courant<1.0)) throw std::invalid_argument("material violates fixed-step Maxwell1D CFL limit");
    if(config_.boundary==Boundary1D::pml
       && (begin<config_.pml_cells || end>config_.cells-config_.pml_cells)
       && (dispersive || epsilon_r!=config_.epsilon_r || conductivity!=0.0))
        throw std::invalid_argument("1-D PML requires unchanged nondispersive background material");
}

void Maxwell1D::set_debye_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                                   double delta_epsilon_r,double relaxation_time_s,double conductivity){
    validate_material(begin_cell,end_cell,epsilon_infinity_r,conductivity,true);
    if(!std::isfinite(delta_epsilon_r)||!std::isfinite(relaxation_time_s))
        throw std::invalid_argument("non-finite Debye parameters");
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_infinity_r>0.0)||!(delta_epsilon_r>=0.0)
       ||!(relaxation_time_s>0.0)||!(conductivity>=0.0)) throw std::invalid_argument("invalid Debye material parameters");
    clear_dispersion(begin_cell,end_cell);
    const double q=dt_/(2.0*relaxation_time_s);
    const double a=(1.0-q)/(1.0+q);
    const double b=(epsilon0*delta_epsilon_r*q)/(1.0+q);
    for(std::size_t i=begin_cell;i<end_cell;++i){
        epsilon_[i]=epsilon0*epsilon_infinity_r;
        conductivity_[i]=conductivity;
        dispersion_[i]=DispersionModel1D::debye;
        dispersion_a_[i]=a;
        dispersion_b_[i]=b;
    }
    preprocess_material_coefficients();
}

void Maxwell1D::set_drude_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                                   double plasma_frequency_rad_s,double collision_frequency_rad_s,double conductivity){
    validate_material(begin_cell,end_cell,epsilon_infinity_r,conductivity,true);
    if(!std::isfinite(plasma_frequency_rad_s)||!std::isfinite(collision_frequency_rad_s))
        throw std::invalid_argument("non-finite Drude parameters");
    const double courant=dt_/(config_.dx*std::sqrt(epsilon0*epsilon_infinity_r*mu_));
    const double pole_dt=plasma_frequency_rad_s*dt_;
    if(!(4.0*courant*courant+pole_dt*pole_dt/epsilon_infinity_r<4.0))
        throw std::invalid_argument("Drude pole violates explicit ADE stability bound");
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_infinity_r>0.0)||!(plasma_frequency_rad_s>=0.0)
       ||!(collision_frequency_rad_s>=0.0)||!(conductivity>=0.0)) throw std::invalid_argument("invalid Drude material parameters");
    clear_dispersion(begin_cell,end_cell);
    const double denom=1.0+0.5*collision_frequency_rad_s*dt_;
    const double a=2.0/denom;
    const double b=-(1.0-0.5*collision_frequency_rad_s*dt_)/denom;
    const double c=epsilon0*plasma_frequency_rad_s*plasma_frequency_rad_s*dt_*dt_/denom;
    for(std::size_t i=begin_cell;i<end_cell;++i){
        epsilon_[i]=epsilon0*epsilon_infinity_r;
        conductivity_[i]=conductivity;
        dispersion_[i]=DispersionModel1D::drude;
        dispersion_a_[i]=a;
        dispersion_b_[i]=b;
        dispersion_c_[i]=c;
    }
    preprocess_material_coefficients();
}

void Maxwell1D::set_lorentz_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                                     double delta_epsilon_r,double resonance_frequency_rad_s,
                                     double damping_rad_s,double conductivity){
    validate_material(begin_cell,end_cell,epsilon_infinity_r,conductivity,true);
    if(!std::isfinite(delta_epsilon_r)||!std::isfinite(resonance_frequency_rad_s)||!std::isfinite(damping_rad_s))
        throw std::invalid_argument("non-finite Lorentz parameters");
    const double courant=dt_/(config_.dx*std::sqrt(epsilon0*epsilon_infinity_r*mu_));
    const double pole_dt=resonance_frequency_rad_s*dt_;
    if(!(4.0*courant*courant+pole_dt*pole_dt*(1.0+delta_epsilon_r/epsilon_infinity_r)<4.0))
        throw std::invalid_argument("Lorentz pole violates conservative explicit ADE stability bound");
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_infinity_r>0.0)||!(delta_epsilon_r>=0.0)
       ||!(resonance_frequency_rad_s>0.0)||!(damping_rad_s>=0.0)||!(conductivity>=0.0)) {
        throw std::invalid_argument("invalid Lorentz material parameters");
    }
    clear_dispersion(begin_cell,end_cell);
    const double omega_dt=resonance_frequency_rad_s*dt_;
    const double denom=1.0+0.5*damping_rad_s*dt_;
    const double a=(2.0-omega_dt*omega_dt)/denom;
    const double b=-(1.0-0.5*damping_rad_s*dt_)/denom;
    const double c=epsilon0*delta_epsilon_r*omega_dt*omega_dt/denom;
    for(std::size_t i=begin_cell;i<end_cell;++i){
        epsilon_[i]=epsilon0*epsilon_infinity_r;
        conductivity_[i]=conductivity;
        dispersion_[i]=DispersionModel1D::lorentz;
        dispersion_a_[i]=a;
        dispersion_b_[i]=b;
        dispersion_c_[i]=c;
    }
    preprocess_material_coefficients();
}

DispersionModel1D Maxwell1D::dispersion_model(std::size_t cell) const{
    if(cell>=config_.cells) throw std::out_of_range("Maxwell1D cell index out of range");
    return dispersion_[cell];
}

void Maxwell1D::initialize_gaussian(double center_fraction,double width_fraction){
    if(!(width_fraction>0.0)||!std::isfinite(width_fraction)||!std::isfinite(center_fraction))
        throw std::invalid_argument("invalid Maxwell1D Gaussian width/center");
    std::fill(hy_.begin(),hy_.end(),0.0);
    std::fill(polarization_.begin(),polarization_.end(),0.0);
    std::fill(polarization_previous_.begin(),polarization_previous_.end(),0.0);
    time_=0.0;
    const double center=center_fraction*static_cast<double>(config_.cells-1U);
    const double width=width_fraction*static_cast<double>(config_.cells);
    for(std::size_t i=0;i<config_.cells;++i){const double z=(static_cast<double>(i)-center)/width;ez_[i]=std::exp(-0.5*z*z);}
    if(config_.boundary!=Boundary1D::mur1){ez_.front()=0.0;ez_.back()=0.0;}
}

void Maxwell1D::add_soft_source(std::size_t cell,double value){
    if(cell>=config_.cells||!std::isfinite(value)) throw std::invalid_argument("invalid Maxwell1D soft source");
    ez_[cell]+=value;
}
void Maxwell1D::set_hard_source(std::size_t cell,double value){
    if(cell>=config_.cells||!std::isfinite(value)) throw std::invalid_argument("invalid Maxwell1D hard source");
    ez_[cell]=value;
}

double Maxwell1D::local_wave_speed(std::size_t cell) const{
    if(cell>=config_.cells) throw std::out_of_range("Maxwell1D cell index out of range");
    return 1.0/std::sqrt(epsilon_[cell]*mu_);
}

void Maxwell1D::apply_boundary(double old_left,double old_left_adj,double old_right,double old_right_adj){
    if(config_.boundary==Boundary1D::pec||config_.boundary==Boundary1D::pml){ez_.front()=0.0;ez_.back()=0.0;return;}
    const double cl=local_wave_speed(0U),cr=local_wave_speed(config_.cells-1U);
    const double kl=(cl*dt_-config_.dx)/(cl*dt_+config_.dx);
    const double kr=(cr*dt_-config_.dx)/(cr*dt_+config_.dx);
    ez_.front()=old_left_adj+kl*(ez_[1U]-old_left);
    ez_.back()=old_right_adj+kr*(ez_[config_.cells-2U]-old_right);
}

void Maxwell1D::step(std::size_t count){
    for(std::size_t step_i=0;step_i<count;++step_i){
        const double old_left=ez_[0U],old_left_adj=ez_[1U],old_right=ez_.back(),old_right_adj=ez_[config_.cells-2U];
        cfd::core::parallel_for(hy_.size(),[&](std::size_t i){
            hy_[i]=h_decay_[i]*hy_[i]+h_curl_[i]*(ez_[i+1U]-ez_[i]);
        });
        cfd::core::parallel_for(config_.cells-2U,[&](std::size_t j){
            const std::size_t i=j+1U;
            const double old_e=ez_[i];
            const double curl_h=(hy_[i]-hy_[i-1U])/config_.dx;
            switch(dispersion_[i]){
                case DispersionModel1D::none:
                    ez_[i]=e_decay_[i]*old_e+e_curl_[i]*(hy_[i]-hy_[i-1U]);
                    break;
                case DispersionModel1D::debye:{
                    const double p_old=polarization_[i];
                    const double a=dispersion_a_[i];
                    const double b=dispersion_b_[i];
                    const double lhs=epsilon_[i]/dt_+0.5*conductivity_[i]+b/dt_;
                    const double rhs=(epsilon_[i]/dt_-0.5*conductivity_[i])*old_e+curl_h
                        -((a-1.0)*p_old+b*old_e)/dt_;
                    const double new_e=rhs/lhs;
                    const double p_new=a*p_old+b*(new_e+old_e);
                    polarization_previous_[i]=p_old;
                    polarization_[i]=p_new;
                    ez_[i]=new_e;
                    break;
                }
                case DispersionModel1D::drude:
                case DispersionModel1D::lorentz:{
                    const double p_old=polarization_[i];
                    const double p_new=dispersion_a_[i]*p_old
                        +dispersion_b_[i]*polarization_previous_[i]+dispersion_c_[i]*old_e;
                    const double denominator=epsilon_[i]+0.5*conductivity_[i]*dt_;
                    const double numerator=(epsilon_[i]-0.5*conductivity_[i]*dt_)*old_e
                        +dt_*curl_h-(p_new-p_old);
                    polarization_previous_[i]=p_old;
                    polarization_[i]=p_new;
                    ez_[i]=numerator/denominator;
                    break;
                }
            }
        });
        apply_boundary(old_left,old_left_adj,old_right,old_right_adj);
        time_+=dt_;
    }
}

double Maxwell1D::energy() const{
    double e=0.0;for(std::size_t i=0;i<ez_.size();++i)e+=0.5*epsilon_[i]*ez_[i]*ez_[i]*config_.dx;
    for(double v:hy_) e+=0.5*mu_*v*v*config_.dx;
    return e;
}

} // namespace cfd::fdtd
