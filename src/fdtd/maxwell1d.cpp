#include "cfd/solvers/fdtd/maxwell1d.hpp"

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfd::fdtd {
namespace {
constexpr double epsilon0 = 8.8541878128e-12;
constexpr double mu0 = 1.25663706212e-6;

[[nodiscard]] double cpml_depth_e(std::size_t i,std::size_t cells,std::size_t layer) noexcept {
    if(i<layer) return (static_cast<double>(layer)-static_cast<double>(i)-0.5)/static_cast<double>(layer);
    if(i+layer>=cells) return (static_cast<double>(i)-static_cast<double>(cells-layer)+0.5)/static_cast<double>(layer);
    return 0.0;
}
[[nodiscard]] double cpml_depth_h(std::size_t i,std::size_t cells,std::size_t layer) noexcept {
    const double x=static_cast<double>(i)+0.5;
    if(x<static_cast<double>(layer)) return (static_cast<double>(layer)-x)/static_cast<double>(layer);
    if(x>static_cast<double>(cells-layer)) return (x-static_cast<double>(cells-layer))/static_cast<double>(layer);
    return 0.0;
}
}

Maxwell1D::Maxwell1D(Maxwell1DConfig config)
    : config_(config),mu_(mu0*config.mu_r),epsilon_(config.cells,epsilon0*config.epsilon_r),
      conductivity_(config.cells,0.0),e_decay_(config.cells,1.0),e_curl_(config.cells,0.0),
      pml_sigma_e_(config.cells,0.0),h_decay_(config.cells-1U,1.0),h_curl_(config.cells-1U,0.0),
      ez_(config.cells,0.0),hy_(config.cells-1U,0.0),
      dispersion_(config.cells,DispersionModel1D::none),polarization_(config.cells,0.0),
      polarization_previous_(config.cells,0.0),dispersion_a_(config.cells,0.0),
      dispersion_b_(config.cells,0.0),dispersion_c_(config.cells,0.0),
      cpml_kappa_e_(config.cells,1.0),cpml_b_e_(config.cells,1.0),cpml_c_e_(config.cells,0.0),
      cpml_psi_e_(config.cells,0.0),cpml_kappa_h_(config.cells-1U,1.0),
      cpml_b_h_(config.cells-1U,1.0),cpml_c_h_(config.cells-1U,0.0),cpml_psi_h_(config.cells-1U,0.0),
      lumped_epsilon_add_(config.cells,0.0),lumped_sigma_add_(config.cells,0.0),
      lumped_inductor_drive_(config.cells,0.0),lumped_inductor_current_density_(config.cells,0.0) {
    if(config.cells<4U||!(config.dx>0.0)||!(config.courant>0.0&&config.courant<1.0)||
       !(config.epsilon_r>0.0)||!(config.mu_r>0.0)) throw std::invalid_argument("invalid Maxwell1D configuration");
    if((config.boundary==Boundary1D::pml||config.boundary==Boundary1D::cpml)&&
       (config.pml_cells<2U||2U*config.pml_cells+4U>config.cells||!(config.pml_order>=1.0)||
        !(config.pml_target_reflection>0.0&&config.pml_target_reflection<1.0)))
        throw std::invalid_argument("invalid Maxwell1D PML configuration");
    if(config.boundary==Boundary1D::cpml&&
       (!(config.cpml_kappa_max>=1.0)||!(config.cpml_alpha_fraction>=0.0)))
        throw std::invalid_argument("invalid Maxwell1D CPML configuration");
    const double wave_speed=1.0/std::sqrt(epsilon_.front()*mu_);
    dt_=config.courant*config.dx/wave_speed;
    initialize_pml_profile();
    initialize_cpml_profile();
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
        const double depth=std::clamp(cpml_depth_e(i,config_.cells,config_.pml_cells),0.0,1.0);
        if(depth>0.0) pml_sigma_e_[i]=sigma_max*std::pow(depth,config_.pml_order);
    }
}

void Maxwell1D::initialize_cpml_profile(){
    std::fill(cpml_kappa_e_.begin(),cpml_kappa_e_.end(),1.0);
    std::fill(cpml_b_e_.begin(),cpml_b_e_.end(),1.0);
    std::fill(cpml_c_e_.begin(),cpml_c_e_.end(),0.0);
    std::fill(cpml_psi_e_.begin(),cpml_psi_e_.end(),0.0);
    std::fill(cpml_kappa_h_.begin(),cpml_kappa_h_.end(),1.0);
    std::fill(cpml_b_h_.begin(),cpml_b_h_.end(),1.0);
    std::fill(cpml_c_h_.begin(),cpml_c_h_.end(),0.0);
    std::fill(cpml_psi_h_.begin(),cpml_psi_h_.end(),0.0);
    if(config_.boundary!=Boundary1D::cpml) return;

    const double eps_ref=epsilon_.front();
    const double eta=std::sqrt(mu_/eps_ref);
    const double thickness=static_cast<double>(config_.pml_cells)*config_.dx;
    const double sigma_max=-(config_.pml_order+1.0)*std::log(config_.pml_target_reflection)
        /(2.0*eta*thickness);
    const double alpha_max=config_.cpml_alpha_fraction*sigma_max;

    for(std::size_t i=0;i<config_.cells;++i){
        const double depth=std::clamp(cpml_depth_e(i,config_.cells,config_.pml_cells),0.0,1.0);
        if(!(depth>0.0)) continue;
        const double grade=std::pow(depth,config_.pml_order);
        const double sigma=sigma_max*grade;
        const double kappa=1.0+(config_.cpml_kappa_max-1.0)*grade;
        const double alpha=alpha_max*(1.0-depth);
        const double b=std::exp(-(sigma/kappa+alpha)*dt_/epsilon_[i]);
        const double denom=sigma*kappa+kappa*kappa*alpha;
        cpml_kappa_e_[i]=kappa;
        cpml_b_e_[i]=b;
        cpml_c_e_[i]=(denom>0.0)?sigma*(b-1.0)/denom:0.0;
    }
    for(std::size_t i=0;i<hy_.size();++i){
        const double depth=std::clamp(cpml_depth_h(i,config_.cells,config_.pml_cells),0.0,1.0);
        if(!(depth>0.0)) continue;
        const double grade=std::pow(depth,config_.pml_order);
        const double sigma_e=sigma_max*grade;
        const double kappa=1.0+(config_.cpml_kappa_max-1.0)*grade;
        const double alpha_e=alpha_max*(1.0-depth);
        const double eps_edge=0.5*(epsilon_[i]+epsilon_[i+1U]);
        const double sigma_m=sigma_e*mu_/eps_edge;
        const double alpha_m=alpha_e*mu_/eps_edge;
        const double b=std::exp(-(sigma_m/kappa+alpha_m)*dt_/mu_);
        const double denom=sigma_m*kappa+kappa*kappa*alpha_m;
        cpml_kappa_h_[i]=kappa;
        cpml_b_h_[i]=b;
        cpml_c_h_[i]=(denom>0.0)?sigma_m*(b-1.0)/denom:0.0;
    }
}

void Maxwell1D::preprocess_material_coefficients(){
    for(std::size_t i=0;i<config_.cells;++i){
        const double effective_epsilon=epsilon_[i]+lumped_epsilon_add_[i];
        const double effective_sigma=conductivity_[i]+lumped_sigma_add_[i]+pml_sigma_e_[i];
        const double loss=effective_sigma*dt_/(2.0*effective_epsilon);
        const double denom=1.0+loss;
        e_decay_[i]=(1.0-loss)/denom;
        e_curl_[i]=(dt_/(effective_epsilon*config_.dx))/denom;
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
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_r>0.0)||!(conductivity>=0.0)) throw std::invalid_argument("invalid Maxwell1D material range");
    clear_dispersion(begin_cell,end_cell);
    for(std::size_t i=begin_cell;i<end_cell;++i){epsilon_[i]=epsilon0*epsilon_r;conductivity_[i]=conductivity;}
    initialize_pml_profile();
    initialize_cpml_profile();
    preprocess_material_coefficients();
}

void Maxwell1D::set_debye_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                                   double delta_epsilon_r,double relaxation_time_s,double conductivity){
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_infinity_r>0.0)||!(delta_epsilon_r>=0.0)
       ||!(relaxation_time_s>0.0)||!(conductivity>=0.0)) throw std::invalid_argument("invalid Debye material parameters");
    for(std::size_t i=begin_cell;i<end_cell;++i) if(lumped_epsilon_add_[i]!=0.0||lumped_sigma_add_[i]!=0.0||lumped_inductor_drive_[i]!=0.0)
        throw std::invalid_argument("dispersive and lumped material models cannot share a Maxwell1D cell");
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
    initialize_pml_profile(); initialize_cpml_profile(); preprocess_material_coefficients();
}

void Maxwell1D::set_drude_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                                   double plasma_frequency_rad_s,double collision_frequency_rad_s,double conductivity){
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_infinity_r>0.0)||!(plasma_frequency_rad_s>=0.0)
       ||!(collision_frequency_rad_s>=0.0)||!(conductivity>=0.0)) throw std::invalid_argument("invalid Drude material parameters");
    for(std::size_t i=begin_cell;i<end_cell;++i) if(lumped_epsilon_add_[i]!=0.0||lumped_sigma_add_[i]!=0.0||lumped_inductor_drive_[i]!=0.0)
        throw std::invalid_argument("dispersive and lumped material models cannot share a Maxwell1D cell");
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
    initialize_pml_profile(); initialize_cpml_profile(); preprocess_material_coefficients();
}

void Maxwell1D::set_lorentz_material(std::size_t begin_cell,std::size_t end_cell,double epsilon_infinity_r,
                                     double delta_epsilon_r,double resonance_frequency_rad_s,
                                     double damping_rad_s,double conductivity){
    if(begin_cell>=end_cell||end_cell>config_.cells||!(epsilon_infinity_r>0.0)||!(delta_epsilon_r>=0.0)
       ||!(resonance_frequency_rad_s>0.0)||!(damping_rad_s>=0.0)||!(conductivity>=0.0)) {
        throw std::invalid_argument("invalid Lorentz material parameters");
    }
    for(std::size_t i=begin_cell;i<end_cell;++i) if(lumped_epsilon_add_[i]!=0.0||lumped_sigma_add_[i]!=0.0||lumped_inductor_drive_[i]!=0.0)
        throw std::invalid_argument("dispersive and lumped material models cannot share a Maxwell1D cell");
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
    initialize_pml_profile(); initialize_cpml_profile(); preprocess_material_coefficients();
}

void Maxwell1D::set_tfsf_source(std::size_t first_total_cell,std::function<double(double)> incident_e,double impedance){
    if(first_total_cell<1U||first_total_cell+2U>=config_.cells||!incident_e) throw std::invalid_argument("invalid Maxwell1D TFSF interface");
    if((config_.boundary==Boundary1D::pml||config_.boundary==Boundary1D::cpml)&&
       (first_total_cell<=config_.pml_cells||first_total_cell+config_.pml_cells>=config_.cells))
        throw std::invalid_argument("Maxwell1D TFSF interface must be outside the absorbing layer");
    if(dispersion_[first_total_cell]!=DispersionModel1D::none||
       lumped_epsilon_add_[first_total_cell]!=0.0||lumped_sigma_add_[first_total_cell]!=0.0||
       lumped_inductor_drive_[first_total_cell]!=0.0)
        throw std::invalid_argument("Maxwell1D TFSF interface must be in a nondispersive field cell");
    if(!(impedance>0.0)) impedance=std::sqrt(mu_/epsilon_[first_total_cell]);
    tfsf_enabled_=true;
    tfsf_cell_=first_total_cell;
    tfsf_incident_e_=std::move(incident_e);
    tfsf_impedance_=impedance;
}
void Maxwell1D::clear_tfsf_source() noexcept {
    tfsf_enabled_=false;
    tfsf_incident_e_={};
    tfsf_impedance_=0.0;
}

void Maxwell1D::set_parallel_lumped_rlc(std::size_t cell,double length,double area,double resistance,double inductance,double capacitance){
    if(cell==0U||cell+1U>=config_.cells||!(length>0.0)||!(area>0.0)||!(capacitance>=0.0)||
       (!(resistance>0.0)&&!std::isinf(resistance))||(!(inductance>0.0)&&!std::isinf(inductance)))
        throw std::invalid_argument("invalid Maxwell1D lumped RLC parameters");
    if(dispersion_[cell]!=DispersionModel1D::none) throw std::invalid_argument("lumped RLC cannot overlap a dispersive Maxwell1D cell");
    lumped_epsilon_add_[cell]=capacitance*length/area;
    lumped_sigma_add_[cell]=std::isfinite(resistance)?length/(resistance*area):0.0;
    lumped_inductor_drive_[cell]=std::isfinite(inductance)?length/(inductance*area):0.0;
    lumped_inductor_current_density_[cell]=0.0;
    preprocess_material_coefficients();
}
void Maxwell1D::clear_lumped_rlc(std::size_t cell){
    if(cell>=config_.cells) throw std::out_of_range("Maxwell1D lumped cell out of range");
    lumped_epsilon_add_[cell]=0.0;
    lumped_sigma_add_[cell]=0.0;
    lumped_inductor_drive_[cell]=0.0;
    lumped_inductor_current_density_[cell]=0.0;
    preprocess_material_coefficients();
}
double Maxwell1D::lumped_inductor_current_density(std::size_t cell) const{
    if(cell>=config_.cells) throw std::out_of_range("Maxwell1D lumped cell out of range");
    return lumped_inductor_current_density_[cell];
}

DispersionModel1D Maxwell1D::dispersion_model(std::size_t cell) const{
    if(cell>=config_.cells) throw std::out_of_range("Maxwell1D cell index out of range");
    return dispersion_[cell];
}

void Maxwell1D::initialize_gaussian(double center_fraction,double width_fraction){
    if(!(width_fraction>0.0)) throw std::invalid_argument("invalid Maxwell1D Gaussian width");
    const double center=center_fraction*static_cast<double>(config_.cells-1U);
    const double width=width_fraction*static_cast<double>(config_.cells);
    for(std::size_t i=0;i<config_.cells;++i){const double z=(static_cast<double>(i)-center)/width;ez_[i]=std::exp(-0.5*z*z);}
    if(config_.boundary==Boundary1D::pec){ez_.front()=0.0;ez_.back()=0.0;}
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
    if(config_.boundary==Boundary1D::pec||config_.boundary==Boundary1D::pml||config_.boundary==Boundary1D::cpml){
        ez_.front()=0.0; ez_.back()=0.0; return;
    }
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
            double derivative=(ez_[i]-ez_[i+1U])/config_.dx;
            if(config_.boundary==Boundary1D::cpml){
                cpml_psi_h_[i]=cpml_b_h_[i]*cpml_psi_h_[i]+cpml_c_h_[i]*derivative;
                derivative=derivative/cpml_kappa_h_[i]+cpml_psi_h_[i];
            }
            hy_[i]=h_decay_[i]*hy_[i]+h_curl_[i]*config_.dx*derivative;
        });
        if(tfsf_enabled_){
            const std::size_t edge=tfsf_cell_-1U;
            const double incident_e=tfsf_incident_e_(time_);
            hy_[edge]-=dt_/(mu_*config_.dx)*incident_e;
        }

        cfd::core::parallel_for(config_.cells-2U,[&](std::size_t j){
            const std::size_t i=j+1U;
            const double old_e=ez_[i];
            if(lumped_inductor_drive_[i]!=0.0)
                lumped_inductor_current_density_[i]+=dt_*lumped_inductor_drive_[i]*old_e;
            double curl_h=(hy_[i-1U]-hy_[i])/config_.dx;
            if(config_.boundary==Boundary1D::cpml){
                cpml_psi_e_[i]=cpml_b_e_[i]*cpml_psi_e_[i]+cpml_c_e_[i]*curl_h;
                curl_h=curl_h/cpml_kappa_e_[i]+cpml_psi_e_[i];
            }
            switch(dispersion_[i]){
                case DispersionModel1D::none:
                    ez_[i]=e_decay_[i]*old_e+e_curl_[i]*config_.dx*
                        (curl_h-lumped_inductor_current_density_[i]);
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
        if(tfsf_enabled_){
            const double incident_e_half=tfsf_incident_e_(time_+0.5*dt_+0.5*config_.dx*std::sqrt(epsilon_[tfsf_cell_]*mu_));
            const double incident_h=-incident_e_half/tfsf_impedance_;
            ez_[tfsf_cell_]+=e_curl_[tfsf_cell_]*incident_h;
        }
        apply_boundary(old_left,old_left_adj,old_right,old_right_adj);
        time_+=dt_;
    }
}

double Maxwell1D::energy() const{
    double e=0.0;
    for(std::size_t i=0;i<ez_.size();++i)
        e+=0.5*(epsilon_[i]+lumped_epsilon_add_[i])*ez_[i]*ez_[i]*config_.dx;
    for(double v:hy_) e+=0.5*mu_*v*v*config_.dx;
    return e;
}

} // namespace cfd::fdtd
