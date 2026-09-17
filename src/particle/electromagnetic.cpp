#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace cfd::particle {
namespace {
constexpr double epsilon0=8.8541878128e-12;
constexpr double mu0=1.25663706212e-6;

constexpr double speed_of_light=299792458.0;
constexpr double elementary_charge=1.602176634e-19;

Vec3 add(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Vec3 scale(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
Vec3 cross(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double norm2(Vec3 a){return a.x*a.x+a.y*a.y+a.z*a.z;}
double dot(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}

}

void boris_push(ChargedParticle& particle,const ElectromagneticField& field,double dt_s){
    if(!(particle.mass_kg>0.0)||!std::isfinite(particle.mass_kg)||!std::isfinite(particle.charge_c)||!(dt_s>0.0)||!std::isfinite(dt_s))
        throw std::invalid_argument("invalid Boris particle parameters");
    const double half_qdt_m=0.5*particle.charge_c*dt_s/particle.mass_kg;
    const Vec3 v_minus=add(particle.velocity_m_per_s,scale(field.electric_v_per_m,half_qdt_m));
    const Vec3 t=scale(field.magnetic_t,half_qdt_m);const Vec3 s=scale(t,2.0/(1.0+norm2(t)));
    const Vec3 v_prime=add(v_minus,cross(v_minus,t));
    const Vec3 v_plus=add(v_minus,cross(v_prime,s));
    particle.velocity_m_per_s=add(v_plus,scale(field.electric_v_per_m,half_qdt_m));
    particle.position_m=add(particle.position_m,scale(particle.velocity_m_per_s,dt_s));
}

void track_particle(ChargedParticle& particle,const FieldSampler& field,double start_time_s,double dt_s,std::size_t steps){
    if(!field)throw std::invalid_argument("particle tracker requires a field sampler");
    double time=start_time_s;for(std::size_t i=0;i<steps;++i){boris_push(particle,field(particle.position_m,time+0.5*dt_s),dt_s);time+=dt_s;}
}


double lorentz_gamma(const ChargedParticle& particle){
    const double v2=norm2(particle.velocity_m_per_s),c2=speed_of_light*speed_of_light;if(!(v2<c2))throw std::invalid_argument("particle speed must be below c");return 1.0/std::sqrt(1.0-v2/c2);
}

double kinetic_energy_j(const ChargedParticle& particle){if(!(particle.mass_kg>0.0))throw std::invalid_argument("particle mass must be positive");return (lorentz_gamma(particle)-1.0)*particle.mass_kg*speed_of_light*speed_of_light;}

void relativistic_boris_push(ChargedParticle& particle,const ElectromagneticField& field,double dt_s){
    if(!(particle.mass_kg>0.0)||!std::isfinite(particle.mass_kg)||!std::isfinite(particle.charge_c)||!(dt_s>0.0)||!std::isfinite(dt_s))throw std::invalid_argument("invalid relativistic Boris parameters");
    const double gamma0=lorentz_gamma(particle);Vec3 u=scale(particle.velocity_m_per_s,gamma0);const double half_qdt_m=0.5*particle.charge_c*dt_s/particle.mass_kg;const Vec3 u_minus=add(u,scale(field.electric_v_per_m,half_qdt_m));const double gamma_minus=std::sqrt(1.0+norm2(u_minus)/(speed_of_light*speed_of_light));const Vec3 t=scale(field.magnetic_t,half_qdt_m/gamma_minus);const Vec3 svec=scale(t,2.0/(1.0+norm2(t)));const Vec3 u_prime=add(u_minus,cross(u_minus,t));const Vec3 u_plus=add(u_minus,cross(u_prime,svec));u=add(u_plus,scale(field.electric_v_per_m,half_qdt_m));const double gamma_new=std::sqrt(1.0+norm2(u)/(speed_of_light*speed_of_light));particle.velocity_m_per_s=scale(u,1.0/gamma_new);particle.position_m=add(particle.position_m,scale(particle.velocity_m_per_s,dt_s));
}

void vay_push(ChargedParticle& particle,const ElectromagneticField& field,double dt_s){
    if(!(particle.mass_kg>0.0)||!std::isfinite(particle.mass_kg)||!std::isfinite(particle.charge_c)||!(dt_s>0.0)||!std::isfinite(dt_s))
        throw std::invalid_argument("invalid Vay particle parameters");
    const double gamma0=lorentz_gamma(particle);
    const Vec3 u0=scale(particle.velocity_m_per_s,gamma0);
    const double half_qdt_m=0.5*particle.charge_c*dt_s/particle.mass_kg;
    const Vec3 tau=scale(field.magnetic_t,half_qdt_m);
    const Vec3 u_minus=add(u0,scale(field.electric_v_per_m,half_qdt_m));
    const Vec3 v0=particle.velocity_m_per_s;
    const Vec3 u_prime=add(add(u_minus,scale(field.electric_v_per_m,half_qdt_m)),scale(cross(v0,field.magnetic_t),half_qdt_m));
    const double gamma_prime=std::sqrt(1.0+norm2(u_prime)/(speed_of_light*speed_of_light));
    const double tau2=norm2(tau);
    const double sigma=gamma_prime*gamma_prime-tau2;
    const double u_star=dot(u_prime,tau)/speed_of_light;
    const double gamma_new=std::sqrt(0.5*(sigma+std::sqrt(sigma*sigma+4.0*(tau2+u_star*u_star))));
    const Vec3 t=scale(tau,1.0/gamma_new);
    const double t2=norm2(t);
    const Vec3 numerator=add(add(u_prime,cross(u_prime,t)),scale(t,dot(u_prime,t)));
    const Vec3 u_new=scale(numerator,1.0/(1.0+t2));
    const double gamma_final=std::sqrt(1.0+norm2(u_new)/(speed_of_light*speed_of_light));
    particle.velocity_m_per_s=scale(u_new,1.0/gamma_final);
    particle.position_m=add(particle.position_m,scale(particle.velocity_m_per_s,dt_s));
}

CollisionStatistics apply_monte_carlo_collisions(std::vector<ChargedParticle>& particles,double dt_s,const NeutralCollisionModel& model){
    if(!(dt_s>0.0)||!std::isfinite(dt_s)||!(model.neutral_density_per_m3>=0.0)||!(model.elastic_cross_section_m2>=0.0)
       ||!(model.ionization_cross_section_m2>=0.0)||!(model.ionization_energy_ev>=0.0))
        throw std::invalid_argument("invalid Monte-Carlo collision model");
    std::mt19937_64 rng(model.random_seed);
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    CollisionStatistics stats;
    const std::size_t initial_count=particles.size();
    for(std::size_t i=0;i<initial_count;++i){
        ChargedParticle& p=particles[i];
        if(!(p.mass_kg>0.0))throw std::invalid_argument("invalid collision particle mass");
        const double speed=std::sqrt(norm2(p.velocity_m_per_s));
        const double sigma_total=model.elastic_cross_section_m2+model.ionization_cross_section_m2;
        if(speed==0.0||sigma_total==0.0||model.neutral_density_per_m3==0.0)continue;
        const double probability=1.0-std::exp(-model.neutral_density_per_m3*sigma_total*speed*dt_s);
        if(uniform(rng)>probability)continue;
        const bool ionize=uniform(rng)<model.ionization_cross_section_m2/sigma_total;
        if(ionize){
            const double energy=kinetic_energy_j(p);
            const double loss=model.ionization_energy_ev*elementary_charge;
            if(energy>loss){
                const double new_energy=energy-loss;
                const double gamma_new=1.0+new_energy/(p.mass_kg*speed_of_light*speed_of_light);
                const double beta2=std::max(0.0,1.0-1.0/(gamma_new*gamma_new));
                const double new_speed=speed_of_light*std::sqrt(beta2);
                const double factor=new_speed/speed;
                p.velocity_m_per_s=scale(p.velocity_m_per_s,factor);
                const double macro_weight=p.weight;
                ChargedParticle secondary=p;
                secondary.velocity_m_per_s=scale(p.velocity_m_per_s,-0.05);
                stats.energy_loss_j+=loss*macro_weight;
                particles.push_back(secondary);
                ++stats.ionization_events;
                continue;
            }
        }
        const double mu=2.0*uniform(rng)-1.0;
        const double phi=2.0*std::numbers::pi*uniform(rng);
        const double sint=std::sqrt(std::max(0.0,1.0-mu*mu));
        p.velocity_m_per_s={speed*sint*std::cos(phi),speed*sint*std::sin(phi),speed*mu};
        ++stats.elastic_events;
    }
    return stats;
}

double secondary_electron_yield(double energy_ev,const SecondaryEmissionModel& model){
    if(!(model.maximum_yield>=0.0)||!(model.energy_at_maximum_ev>model.threshold_energy_ev)||!(model.threshold_energy_ev>=0.0)||!std::isfinite(energy_ev))
        throw std::invalid_argument("invalid secondary-emission model/energy");
    if(energy_ev<=model.threshold_energy_ev)return 0.0;
    const double x=(energy_ev-model.threshold_energy_ev)/(model.energy_at_maximum_ev-model.threshold_energy_ev);
    return model.maximum_yield*x*std::exp(1.0-x);
}

ParticleWallInteraction apply_particle_box_boundary(ChargedParticle& particle,const AxisAlignedParticleBox& box,ParticleWallMode mode,const SecondaryEmissionModel* secondary){
    if(!(box.maximum_m.x>box.minimum_m.x&&box.maximum_m.y>box.minimum_m.y&&box.maximum_m.z>box.minimum_m.z))
        throw std::invalid_argument("invalid particle boundary box");
    ParticleWallInteraction result;
    const double incident_ev=kinetic_energy_j(particle)/elementary_charge;
    auto outside=[&](double value,double lo,double hi){return value<lo||value>hi;};
    if(!outside(particle.position_m.x,box.minimum_m.x,box.maximum_m.x)&&!outside(particle.position_m.y,box.minimum_m.y,box.maximum_m.y)&&!outside(particle.position_m.z,box.minimum_m.z,box.maximum_m.z))return result;
    result.impacted=true;
    auto absorb_or_reflect=[&](double& position,double& velocity,double lo,double hi,Vec3 low_normal,Vec3 high_normal){if(position>=lo&&position<=hi)return;result.impact_normal=position<lo?low_normal:high_normal;if(mode==ParticleWallMode::absorb){result.alive=false;return;}while(position<lo||position>hi){if(position<lo){position=2.0*lo-position;velocity=-velocity;}else if(position>hi){position=2.0*hi-position;velocity=-velocity;}}};
    absorb_or_reflect(particle.position_m.x,particle.velocity_m_per_s.x,box.minimum_m.x,box.maximum_m.x,{-1,0,0},{1,0,0});if(result.alive)absorb_or_reflect(particle.position_m.y,particle.velocity_m_per_s.y,box.minimum_m.y,box.maximum_m.y,{0,-1,0},{0,1,0});if(result.alive)absorb_or_reflect(particle.position_m.z,particle.velocity_m_per_s.z,box.minimum_m.z,box.maximum_m.z,{0,0,-1},{0,0,1});if(!result.alive&&secondary)result.secondary_macro_weight=particle.weight*secondary_electron_yield(incident_ev,*secondary);return result;
}

std::vector<double> periodic_electric_field_from_charge_density(std::span<const double> rho,double length_m,
                                                                 double relative_permittivity,bool remove_mean){
    const std::size_t n=rho.size();if(n<2U||!(length_m>0.0)||!(relative_permittivity>0.0))throw std::invalid_argument("invalid periodic Poisson grid");
    double mean=0.0;for(double value:rho){if(!std::isfinite(value))throw std::invalid_argument("non-finite charge density");mean+=value;}mean/=static_cast<double>(n);
    if(!remove_mean&&std::abs(mean)>1.0e-12*std::max(1.0,*std::max_element(rho.begin(),rho.end(),[](double a,double b){return std::abs(a)<std::abs(b);})))
        throw std::invalid_argument("periodic Poisson requires mean-zero charge density");
    using Complex=std::complex<double>;std::vector<Complex> field_hat(n,Complex{});const double epsilon=epsilon0*relative_permittivity;
    for(std::size_t mode=1U;mode<n;++mode){
        const long signed_mode=(mode<=n/2U)?static_cast<long>(mode):static_cast<long>(mode)-static_cast<long>(n);
        const double k=2.0*std::numbers::pi*static_cast<double>(signed_mode)/length_m;Complex rho_hat{};
        for(std::size_t j=0;j<n;++j){const double angle=-2.0*std::numbers::pi*static_cast<double>(mode*j)/static_cast<double>(n);rho_hat+=(rho[j]-(remove_mean?mean:0.0))*Complex{std::cos(angle),std::sin(angle)};}
        rho_hat/=static_cast<double>(n);field_hat[mode]=Complex{0.0,-1.0}*rho_hat/(k*epsilon);
    }
    std::vector<double> field(n,0.0);
    for(std::size_t j=0;j<n;++j){Complex value{};for(std::size_t mode=1U;mode<n;++mode){const double angle=2.0*std::numbers::pi*static_cast<double>(mode*j)/static_cast<double>(n);value+=field_hat[mode]*Complex{std::cos(angle),std::sin(angle)};}field[j]=value.real();}
    return field;
}

ElectrostaticPic1D::ElectrostaticPic1D(ElectrostaticPic1DConfig config):config_(config){
    if(config_.grid_points<2U||!(config_.length_m>0.0)||!(config_.relative_permittivity>0.0)||!(config_.dt_s>0.0))throw std::invalid_argument("invalid electrostatic PIC configuration");
    charge_density_.assign(config_.grid_points,0.0);electric_field_.assign(config_.grid_points,0.0);
}
void ElectrostaticPic1D::set_particles(std::vector<PicParticle1D> particles){for(auto& p:particles){if(!(p.mass_kg>0.0)||!std::isfinite(p.charge_c)||!std::isfinite(p.weight))throw std::invalid_argument("invalid PIC particle");p.position_m=wrap_position(p.position_m);}particles_=std::move(particles);}
double ElectrostaticPic1D::wrap_position(double x) const {x=std::fmod(x,config_.length_m);if(x<0.0)x+=config_.length_m;return x;}
void ElectrostaticPic1D::deposit_and_solve(){
    std::fill(charge_density_.begin(),charge_density_.end(),0.0);const double dx=config_.length_m/static_cast<double>(config_.grid_points);
    for(const auto& p:particles_){const double u=wrap_position(p.position_m)/dx;const std::size_t left=static_cast<std::size_t>(std::floor(u))%config_.grid_points;const std::size_t right=(left+1U)%config_.grid_points;const double fraction=u-std::floor(u);const double density=p.charge_c*p.weight/dx;charge_density_[left]+=density*(1.0-fraction);charge_density_[right]+=density*fraction;}
    electric_field_=periodic_electric_field_from_charge_density(charge_density_,config_.length_m,config_.relative_permittivity,config_.neutralize_mean_charge);
}
double ElectrostaticPic1D::interpolate_field(double position_m) const {const double dx=config_.length_m/static_cast<double>(config_.grid_points),u=wrap_position(position_m)/dx;const std::size_t left=static_cast<std::size_t>(std::floor(u))%config_.grid_points,right=(left+1U)%config_.grid_points;const double fraction=u-std::floor(u);return electric_field_[left]*(1.0-fraction)+electric_field_[right]*fraction;}
void ElectrostaticPic1D::step(std::size_t steps){for(std::size_t step_index=0;step_index<steps;++step_index){deposit_and_solve();for(auto& p:particles_){p.velocity_m_per_s+=(p.charge_c/p.mass_kg)*interpolate_field(p.position_m)*config_.dt_s;p.position_m=wrap_position(p.position_m+p.velocity_m_per_s*config_.dt_s);}time_s_+=config_.dt_s;}deposit_and_solve();}


namespace {
std::vector<double> deposit_cic_charge_density(std::span<const ChargedParticle> particles,std::size_t grid_points,double length_m){
    if(grid_points<2U||!(length_m>0.0))throw std::invalid_argument("invalid CIC charge grid");
    std::vector<double> rho(grid_points,0.0);const double cell=length_m/static_cast<double>(grid_points);
    auto wrap=[&](double value){value=std::fmod(value,length_m);if(value<0.0)value+=length_m;return value;};
    for(const auto& particle:particles){
        if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid charged particle for deposition");
        const double u=wrap(particle.position_m.x)/cell;const auto left_index=static_cast<std::size_t>(std::floor(u))%grid_points;const auto right_index=(left_index+1U)%grid_points;const double fraction=u-std::floor(u);const double density=particle.charge_c*particle.weight/cell;
        rho[left_index]+=density*(1.0-fraction);rho[right_index]+=density*fraction;
    }
    return rho;
}

std::complex<double> fourier_coefficient(std::span<const double> values,std::size_t mode){
    using Complex=std::complex<double>;Complex coefficient{};const std::size_t count=values.size();
    for(std::size_t index=0;index<count;++index){const double angle=-2.0*std::numbers::pi*static_cast<double>(mode*index)/static_cast<double>(count);coefficient+=values[index]*Complex{std::cos(angle),std::sin(angle)};}
    return coefficient/static_cast<double>(count);
}
}

ChargeConservingCurrent1D deposit_charge_conserving_current_1d(std::span<const ChargedParticle> particles_old,
                                                               std::span<const ChargedParticle> particles_new,
                                                               std::size_t grid_points,double length_m,double dt_s){
    if(particles_old.size()!=particles_new.size())throw std::invalid_argument("old/new particle arrays must have equal size");
    if(grid_points<2U||!(length_m>0.0)||!(dt_s>0.0)||!std::isfinite(dt_s))throw std::invalid_argument("invalid charge-conserving deposition parameters");
    ChargeConservingCurrent1D result;result.charge_density_old_c_per_m3=deposit_cic_charge_density(particles_old,grid_points,length_m);result.charge_density_new_c_per_m3=deposit_cic_charge_density(particles_new,grid_points,length_m);result.current_x_a_per_m2.assign(grid_points,0.0);
    using Complex=std::complex<double>;std::vector<Complex> current_hat(grid_points,Complex{}),divergence_hat(grid_points,Complex{});const double two_pi=2.0*std::numbers::pi;
    for(std::size_t mode=1U;mode<grid_points;++mode){
        const long signed_mode=(mode<=grid_points/2U)?static_cast<long>(mode):static_cast<long>(mode)-static_cast<long>(grid_points);
        const double wavenumber=two_pi*static_cast<double>(signed_mode)/length_m;const Complex delta_hat=fourier_coefficient(result.charge_density_new_c_per_m3,mode)-fourier_coefficient(result.charge_density_old_c_per_m3,mode);
        current_hat[mode]=Complex{0.0,1.0}*delta_hat/(wavenumber*dt_s);divergence_hat[mode]=Complex{0.0,1.0}*wavenumber*current_hat[mode];
    }
    std::vector<double> divergence(grid_points,0.0);
    for(std::size_t index=0;index<grid_points;++index){Complex current_value{},divergence_value{};for(std::size_t mode=1U;mode<grid_points;++mode){const double angle=two_pi*static_cast<double>(mode*index)/static_cast<double>(grid_points);const Complex phase{std::cos(angle),std::sin(angle)};current_value+=current_hat[mode]*phase;divergence_value+=divergence_hat[mode]*phase;}result.current_x_a_per_m2[index]=current_value.real();divergence[index]=divergence_value.real();}
    double residual=0.0;for(std::size_t index=0;index<grid_points;++index){const double continuity=(result.charge_density_new_c_per_m3[index]-result.charge_density_old_c_per_m3[index])/dt_s+divergence[index];residual=std::max(residual,std::abs(continuity));}
    result.continuity_linf_residual=residual;return result;
}

ElectromagneticPic1D::ElectromagneticPic1D(ElectromagneticPic1DConfig config):config_(config){
    if(config_.grid_points<2U||!(config_.length_m>0.0)||!(config_.relative_permittivity>0.0)||!(config_.relative_permeability>0.0)||!(config_.dt_s>0.0))throw std::invalid_argument("invalid electromagnetic PIC configuration");
    const double courant=config_.dt_s*medium_light_speed()/dx();if(courant>1.0)throw std::invalid_argument("1-D electromagnetic PIC CFL must be <= 1");
    charge_density_.assign(config_.grid_points,0.0);current_x_.assign(config_.grid_points,0.0);current_y_.assign(config_.grid_points,0.0);current_z_.assign(config_.grid_points,0.0);electric_x_.assign(config_.grid_points,0.0);electric_y_.assign(config_.grid_points,0.0);electric_z_.assign(config_.grid_points,0.0);magnetic_y_.assign(config_.grid_points,0.0);magnetic_z_.assign(config_.grid_points,0.0);
}

double ElectromagneticPic1D::epsilon() const noexcept{return epsilon0*config_.relative_permittivity;}
double ElectromagneticPic1D::mu() const noexcept{return (4.0e-7*std::numbers::pi)*config_.relative_permeability;}
double ElectromagneticPic1D::medium_light_speed() const noexcept{return 1.0/std::sqrt(epsilon()*mu());}

double ElectromagneticPic1D::wrap_position(double position_m) const {position_m=std::fmod(position_m,config_.length_m);if(position_m<0.0)position_m+=config_.length_m;return position_m;}

void ElectromagneticPic1D::set_particles(std::vector<ChargedParticle> particles){
    for(auto& particle:particles){if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid electromagnetic PIC particle");particle.position_m.x=wrap_position(particle.position_m.x);}particles_=std::move(particles);deposit_sources();}

void ElectromagneticPic1D::set_transverse_fields(std::span<const double> electric_y_v_per_m,std::span<const double> electric_z_v_per_m,std::span<const double> magnetic_y_t,std::span<const double> magnetic_z_t){
    const std::size_t count=config_.grid_points;if(electric_y_v_per_m.size()!=count||electric_z_v_per_m.size()!=count||magnetic_y_t.size()!=count||magnetic_z_t.size()!=count)throw std::invalid_argument("transverse PIC field arrays must match grid size");
    electric_y_.assign(electric_y_v_per_m.begin(),electric_y_v_per_m.end());electric_z_.assign(electric_z_v_per_m.begin(),electric_z_v_per_m.end());magnetic_y_.assign(magnetic_y_t.begin(),magnetic_y_t.end());magnetic_z_.assign(magnetic_z_t.begin(),magnetic_z_t.end());
}

void ElectromagneticPic1D::initialize_right_traveling_mode(double electric_amplitude_v_per_m,std::size_t mode_index){
    if(mode_index==0U)throw std::invalid_argument("mode index must be non-zero");
    const double cell=dx();const double k=2.0*std::numbers::pi*static_cast<double>(mode_index)/config_.length_m;const double c=medium_light_speed();
    for(std::size_t index=0;index<config_.grid_points;++index){const double x_node=static_cast<double>(index)*cell;const double x_half=(static_cast<double>(index)+0.5)*cell;electric_y_[index]=electric_amplitude_v_per_m*std::sin(k*x_node);electric_z_[index]=0.0;magnetic_y_[index]=0.0;magnetic_z_[index]=(electric_amplitude_v_per_m/c)*std::sin(k*x_half);}
}

void ElectromagneticPic1D::deposit_sources(){
    charge_density_=deposit_cic_charge_density(particles_,config_.grid_points,config_.length_m);std::fill(current_x_.begin(),current_x_.end(),0.0);std::fill(current_y_.begin(),current_y_.end(),0.0);std::fill(current_z_.begin(),current_z_.end(),0.0);
    if(config_.solve_longitudinal_poisson)electric_x_=periodic_electric_field_from_charge_density(charge_density_,config_.length_m,config_.relative_permittivity,config_.neutralize_mean_charge);else std::fill(electric_x_.begin(),electric_x_.end(),0.0);
}

double ElectromagneticPic1D::interpolate_node_field(std::span<const double> nodal,double position_m) const {const double cell=dx();const double u=wrap_position(position_m)/cell;const auto left_index=static_cast<std::size_t>(std::floor(u))%config_.grid_points;const auto right_index=(left_index+1U)%config_.grid_points;const double fraction=u-std::floor(u);return nodal[left_index]*(1.0-fraction)+nodal[right_index]*fraction;}

double ElectromagneticPic1D::interpolate_half_field(std::span<const double> half_cell,double position_m) const {const double cell=dx();double shifted=wrap_position(position_m-0.5*cell);const double u=shifted/cell;const auto left_index=static_cast<std::size_t>(std::floor(u))%config_.grid_points;const auto right_index=(left_index+1U)%config_.grid_points;const double fraction=u-std::floor(u);return half_cell[left_index]*(1.0-fraction)+half_cell[right_index]*fraction;}

ElectromagneticField ElectromagneticPic1D::gather_field(double position_m) const {return {{interpolate_node_field(electric_x_,position_m),interpolate_node_field(electric_y_,position_m),interpolate_node_field(electric_z_,position_m)},{0.0,interpolate_half_field(magnetic_y_,position_m),interpolate_half_field(magnetic_z_,position_m)}};}

void ElectromagneticPic1D::advance_magnetic(double dt_s){const double cell=dx();for(std::size_t index=0;index<config_.grid_points;++index){const auto right_index=(index+1U)%config_.grid_points;magnetic_y_[index]+=dt_s*(electric_z_[right_index]-electric_z_[index])/cell;magnetic_z_[index]-=dt_s*(electric_y_[right_index]-electric_y_[index])/cell;}}

void ElectromagneticPic1D::deposit_transverse_current(std::span<const ChargedParticle> particles_old){
    std::fill(current_y_.begin(),current_y_.end(),0.0);std::fill(current_z_.begin(),current_z_.end(),0.0);const double cell=dx();
    for(std::size_t particle_index=0;particle_index<particles_.size();++particle_index){
        double x_old=particles_old[particle_index].position_m.x;double x_new=particles_[particle_index].position_m.x;double delta=x_new-x_old;if(delta>0.5*config_.length_m)x_old+=config_.length_m;else if(delta<-0.5*config_.length_m)x_new+=config_.length_m;const double x_mid=wrap_position(0.5*(x_old+x_new));
        const double u=x_mid/cell;const auto left_index=static_cast<std::size_t>(std::floor(u))%config_.grid_points;const auto right_index=(left_index+1U)%config_.grid_points;const double fraction=u-std::floor(u);const double scale_current=particles_[particle_index].charge_c*particles_[particle_index].weight/cell;const double vy=0.5*(particles_old[particle_index].velocity_m_per_s.y+particles_[particle_index].velocity_m_per_s.y);const double vz=0.5*(particles_old[particle_index].velocity_m_per_s.z+particles_[particle_index].velocity_m_per_s.z);
        current_y_[left_index]+=scale_current*vy*(1.0-fraction);current_y_[right_index]+=scale_current*vy*fraction;current_z_[left_index]+=scale_current*vz*(1.0-fraction);current_z_[right_index]+=scale_current*vz*fraction;
    }
}

void ElectromagneticPic1D::step(std::size_t steps,const NeutralCollisionModel* collision_model){
    for(std::size_t step_index=0;step_index<steps;++step_index){
        const auto particles_old=particles_;advance_magnetic(0.5*config_.dt_s);
        for(auto& particle:particles_){const auto field=gather_field(particle.position_m.x);vay_push(particle,field,config_.dt_s);particle.position_m.x=wrap_position(particle.position_m.x);}
        const auto cc=deposit_charge_conserving_current_1d(particles_old,particles_,config_.grid_points,config_.length_m,config_.dt_s);charge_density_=cc.charge_density_new_c_per_m3;current_x_=cc.current_x_a_per_m2;last_continuity_residual_=cc.continuity_linf_residual;deposit_transverse_current(particles_old);
        const double cell=dx();const double eps=epsilon();const double permeability=mu();
        for(std::size_t index=0;index<config_.grid_points;++index){const auto left_index=(index+config_.grid_points-1U)%config_.grid_points;const double curl_h_y=-(magnetic_z_[index]-magnetic_z_[left_index])/(permeability*cell);const double curl_h_z=(magnetic_y_[index]-magnetic_y_[left_index])/(permeability*cell);electric_y_[index]+=config_.dt_s*(curl_h_y-current_y_[index])/eps;electric_z_[index]+=config_.dt_s*(curl_h_z-current_z_[index])/eps;}
        advance_magnetic(0.5*config_.dt_s);
        if(config_.solve_longitudinal_poisson)electric_x_=periodic_electric_field_from_charge_density(charge_density_,config_.length_m,config_.relative_permittivity,config_.neutralize_mean_charge);else std::fill(electric_x_.begin(),electric_x_.end(),0.0);
        if(collision_model){const auto stats=apply_monte_carlo_collisions(particles_,config_.dt_s,*collision_model);accumulated_collisions_.elastic_events+=stats.elastic_events;accumulated_collisions_.ionization_events+=stats.ionization_events;accumulated_collisions_.energy_loss_j+=stats.energy_loss_j;for(auto& particle:particles_)particle.position_m.x=wrap_position(particle.position_m.x);}
        time_s_+=config_.dt_s;
    }
}

ElectromagneticPicDiagnostics ElectromagneticPic1D::diagnostics() const {
    ElectromagneticPicDiagnostics result;const double cell=dx();const double eps=epsilon();const double permeability=mu();
    for(std::size_t index=0;index<config_.grid_points;++index){const double e2=electric_x_[index]*electric_x_[index]+electric_y_[index]*electric_y_[index]+electric_z_[index]*electric_z_[index];const double b2=magnetic_y_[index]*magnetic_y_[index]+magnetic_z_[index]*magnetic_z_[index];result.field_energy_j+=0.5*eps*e2*cell+0.5*b2*cell/permeability;}
    for(const auto& particle:particles_)result.particle_kinetic_energy_j+=particle.weight*kinetic_energy_j(particle);
    result.total_energy_j=result.field_energy_j+result.particle_kinetic_energy_j;result.charge_continuity_linf_residual=last_continuity_residual_;result.collision_statistics=accumulated_collisions_;result.particle_count=particles_.size();return result;
}

namespace {
std::size_t flat_index_2d(std::size_t ix,std::size_t iy,std::size_t nx){return iy*nx+ix;}

double wrap_periodic_value(double value,double length){value=std::fmod(value,length);if(value<0.0)value+=length;return value;}

std::complex<double> fourier_coefficient_2d(std::span<const double> values,std::size_t nx,std::size_t ny,
                                            std::size_t mode_x,std::size_t mode_y){
    using Complex=std::complex<double>;Complex coefficient{};const double two_pi=2.0*std::numbers::pi;
    for(std::size_t iy=0;iy<ny;++iy){
        for(std::size_t ix=0;ix<nx;++ix){
            const double angle=-two_pi*(static_cast<double>(mode_x*ix)/static_cast<double>(nx)+static_cast<double>(mode_y*iy)/static_cast<double>(ny));
            coefficient+=values[flat_index_2d(ix,iy,nx)]*Complex{std::cos(angle),std::sin(angle)};
        }
    }
    return coefficient/static_cast<double>(nx*ny);
}

std::pair<long,long> signed_modes_2d(std::size_t mode_x,std::size_t mode_y,std::size_t nx,std::size_t ny){
    const long sx=(mode_x<=nx/2U)?static_cast<long>(mode_x):static_cast<long>(mode_x)-static_cast<long>(nx);
    const long sy=(mode_y<=ny/2U)?static_cast<long>(mode_y):static_cast<long>(mode_y)-static_cast<long>(ny);
    return {sx,sy};
}
}

std::vector<double> deposit_cic_charge_density_2d(std::span<const PicParticle2D> particles,std::size_t nx,std::size_t ny,
                                                   double length_x_m,double length_y_m){
    if(nx<2U||ny<2U||!(length_x_m>0.0)||!(length_y_m>0.0))throw std::invalid_argument("invalid 2-D CIC grid");
    std::vector<double> rho(nx*ny,0.0);const double dx_value=length_x_m/static_cast<double>(nx);const double dy_value=length_y_m/static_cast<double>(ny);
    for(const auto& particle:particles){
        if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid 2-D PIC particle");
        const double ux=wrap_periodic_value(particle.position_m.x,length_x_m)/dx_value;const double uy=wrap_periodic_value(particle.position_m.y,length_y_m)/dy_value;
        const auto ix0=static_cast<std::size_t>(std::floor(ux))%nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%ny;const auto ix1=(ix0+1U)%nx;const auto iy1=(iy0+1U)%ny;
        const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);const double density=particle.charge_c*particle.weight/(dx_value*dy_value);
        rho[flat_index_2d(ix0,iy0,nx)]+=density*(1.0-fx)*(1.0-fy);
        rho[flat_index_2d(ix1,iy0,nx)]+=density*fx*(1.0-fy);
        rho[flat_index_2d(ix0,iy1,nx)]+=density*(1.0-fx)*fy;
        rho[flat_index_2d(ix1,iy1,nx)]+=density*fx*fy;
    }
    return rho;
}

ElectrostaticField2D periodic_electric_field_from_charge_density_2d(std::span<const double> rho,std::size_t nx,std::size_t ny,
                                                                     double length_x_m,double length_y_m,double relative_permittivity,
                                                                     bool remove_mean){
    if(rho.size()!=nx*ny||nx<2U||ny<2U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(relative_permittivity>0.0))throw std::invalid_argument("invalid 2-D periodic Poisson grid");
    double mean_value=0.0;double max_abs_value=0.0;for(double value:rho){if(!std::isfinite(value))throw std::invalid_argument("non-finite 2-D charge density");mean_value+=value;max_abs_value=std::max(max_abs_value,std::abs(value));}mean_value/=static_cast<double>(rho.size());
    if(!remove_mean&&std::abs(mean_value)>1.0e-12*std::max(1.0,max_abs_value))throw std::invalid_argument("2-D periodic Poisson requires mean-zero charge density");
    std::vector<double> centered(rho.begin(),rho.end());if(remove_mean)for(double& value:centered)value-=mean_value;
    using Complex=std::complex<double>;std::vector<Complex> ex_hat(nx*ny),ey_hat(nx*ny);const double epsilon=epsilon0*relative_permittivity;const double two_pi=2.0*std::numbers::pi;
    for(std::size_t my=0;my<ny;++my){
        for(std::size_t mx=0;mx<nx;++mx){
            if(mx==0U&&my==0U){continue;}
            const auto [sx,sy]=signed_modes_2d(mx,my,nx,ny);const double kx=two_pi*static_cast<double>(sx)/length_x_m;const double ky=two_pi*static_cast<double>(sy)/length_y_m;const double k2=kx*kx+ky*ky;if(k2==0.0){continue;}
            const Complex rho_hat=fourier_coefficient_2d(centered,nx,ny,mx,my);const Complex factor=Complex{0.0,-1.0}*rho_hat/(epsilon*k2);
            ex_hat[flat_index_2d(mx,my,nx)]=factor*kx;ey_hat[flat_index_2d(mx,my,nx)]=factor*ky;
        }
    }
    ElectrostaticField2D field;field.electric_x_v_per_m.assign(nx*ny,0.0);field.electric_y_v_per_m.assign(nx*ny,0.0);
    for(std::size_t iy=0;iy<ny;++iy){
        for(std::size_t ix=0;ix<nx;++ix){
            Complex ex_value{},ey_value{};
            for(std::size_t my=0;my<ny;++my){
                for(std::size_t mx=0;mx<nx;++mx){
                    const double angle=two_pi*(static_cast<double>(mx*ix)/static_cast<double>(nx)+static_cast<double>(my*iy)/static_cast<double>(ny));const Complex phase{std::cos(angle),std::sin(angle)};
                    ex_value+=ex_hat[flat_index_2d(mx,my,nx)]*phase;ey_value+=ey_hat[flat_index_2d(mx,my,nx)]*phase;
                }
            }
            field.electric_x_v_per_m[flat_index_2d(ix,iy,nx)]=ex_value.real();field.electric_y_v_per_m[flat_index_2d(ix,iy,nx)]=ey_value.real();
        }
    }
    return field;
}

ChargeConservingCurrent2D deposit_charge_conserving_current_2d(std::span<const PicParticle2D> particles_old,
                                                                std::span<const PicParticle2D> particles_new,
                                                                std::size_t nx,std::size_t ny,double length_x_m,
                                                                double length_y_m,double dt_s){
    if(particles_old.size()!=particles_new.size())throw std::invalid_argument("old/new 2-D particle arrays must have equal size");
    if(nx<2U||ny<2U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(dt_s>0.0)||!std::isfinite(dt_s))throw std::invalid_argument("invalid 2-D charge-conserving deposition parameters");
    ChargeConservingCurrent2D result;result.charge_density_old_c_per_m3=deposit_cic_charge_density_2d(particles_old,nx,ny,length_x_m,length_y_m);result.charge_density_new_c_per_m3=deposit_cic_charge_density_2d(particles_new,nx,ny,length_x_m,length_y_m);result.current_x_a_per_m2.assign(nx*ny,0.0);result.current_y_a_per_m2.assign(nx*ny,0.0);
    std::vector<double> delta(nx*ny);for(std::size_t idx=0;idx<delta.size();++idx)delta[idx]=result.charge_density_new_c_per_m3[idx]-result.charge_density_old_c_per_m3[idx];
    using Complex=std::complex<double>;std::vector<Complex> jx_hat(nx*ny),jy_hat(nx*ny),div_hat(nx*ny);const double two_pi=2.0*std::numbers::pi;
    for(std::size_t my=0;my<ny;++my){
        for(std::size_t mx=0;mx<nx;++mx){
            if(mx==0U&&my==0U){continue;}
            const auto [sx,sy]=signed_modes_2d(mx,my,nx,ny);const double kx=two_pi*static_cast<double>(sx)/length_x_m;const double ky=two_pi*static_cast<double>(sy)/length_y_m;const double k2=kx*kx+ky*ky;if(k2==0.0){continue;}
            const Complex delta_hat=fourier_coefficient_2d(delta,nx,ny,mx,my);const Complex factor=Complex{0.0,1.0}*delta_hat/(dt_s*k2);
            jx_hat[flat_index_2d(mx,my,nx)]=factor*kx;jy_hat[flat_index_2d(mx,my,nx)]=factor*ky;div_hat[flat_index_2d(mx,my,nx)]=Complex{0.0,1.0}*(kx*jx_hat[flat_index_2d(mx,my,nx)]+ky*jy_hat[flat_index_2d(mx,my,nx)]);
        }
    }
    std::vector<double> divergence(nx*ny,0.0);
    for(std::size_t iy=0;iy<ny;++iy){
        for(std::size_t ix=0;ix<nx;++ix){
            Complex jx_value{},jy_value{},div_value{};
            for(std::size_t my=0;my<ny;++my){
                for(std::size_t mx=0;mx<nx;++mx){
                    const double angle=two_pi*(static_cast<double>(mx*ix)/static_cast<double>(nx)+static_cast<double>(my*iy)/static_cast<double>(ny));const Complex phase{std::cos(angle),std::sin(angle)};const std::size_t k_index=flat_index_2d(mx,my,nx);
                    jx_value+=jx_hat[k_index]*phase;jy_value+=jy_hat[k_index]*phase;div_value+=div_hat[k_index]*phase;
                }
            }
            const std::size_t real_index=flat_index_2d(ix,iy,nx);result.current_x_a_per_m2[real_index]=jx_value.real();result.current_y_a_per_m2[real_index]=jy_value.real();divergence[real_index]=div_value.real();
        }
    }
    double residual=0.0;for(std::size_t idx=0;idx<delta.size();++idx){const double continuity=delta[idx]/dt_s+divergence[idx];residual=std::max(residual,std::abs(continuity));}result.continuity_linf_residual=residual;return result;
}


ElectrostaticPic2D::ElectrostaticPic2D(ElectrostaticPic2DConfig config):config_(config){
    if(config_.nx<2U||config_.ny<2U||!(config_.length_x_m>0.0)||!(config_.length_y_m>0.0)||!(config_.relative_permittivity>0.0)||!(config_.dt_s>0.0))throw std::invalid_argument("invalid electrostatic PIC 2-D configuration");
    charge_density_.assign(config_.nx*config_.ny,0.0);electric_x_.assign(config_.nx*config_.ny,0.0);electric_y_.assign(config_.nx*config_.ny,0.0);
}

Vec2 ElectrostaticPic2D::wrap_position(Vec2 position_m) const {position_m.x=wrap_periodic_value(position_m.x,config_.length_x_m);position_m.y=wrap_periodic_value(position_m.y,config_.length_y_m);return position_m;}

void ElectrostaticPic2D::set_particles(std::vector<PicParticle2D> particles){for(auto& particle:particles){if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid electrostatic PIC 2-D particle");particle.position_m=wrap_position(particle.position_m);}particles_=std::move(particles);}

void ElectrostaticPic2D::deposit_and_solve(){charge_density_=deposit_cic_charge_density_2d(particles_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m);const auto field=periodic_electric_field_from_charge_density_2d(charge_density_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m,config_.relative_permittivity,config_.neutralize_mean_charge);electric_x_=field.electric_x_v_per_m;electric_y_=field.electric_y_v_per_m;}

double ElectrostaticPic2D::bilinear(std::span<const double> nodal,Vec2 position_m) const {const Vec2 wrapped=wrap_position(position_m);const double ux=wrapped.x/dx();const double uy=wrapped.y/dy();const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;const auto ix1=(ix0+1U)%config_.nx;const auto iy1=(iy0+1U)%config_.ny;const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);return nodal[flat_index_2d(ix0,iy0,config_.nx)]*(1.0-fx)*(1.0-fy)+nodal[flat_index_2d(ix1,iy0,config_.nx)]*fx*(1.0-fy)+nodal[flat_index_2d(ix0,iy1,config_.nx)]*(1.0-fx)*fy+nodal[flat_index_2d(ix1,iy1,config_.nx)]*fx*fy;}

Vec2 ElectrostaticPic2D::gather_electric_field(Vec2 position_m) const {return {bilinear(electric_x_,position_m),bilinear(electric_y_,position_m)};}

void ElectrostaticPic2D::step(std::size_t steps){for(std::size_t step_index=0;step_index<steps;++step_index){deposit_and_solve();for(auto& particle:particles_){const Vec2 field=gather_electric_field(particle.position_m);particle.velocity_m_per_s.x+=(particle.charge_c/particle.mass_kg)*field.x*config_.dt_s;particle.velocity_m_per_s.y+=(particle.charge_c/particle.mass_kg)*field.y*config_.dt_s;particle.position_m.x+=particle.velocity_m_per_s.x*config_.dt_s;particle.position_m.y+=particle.velocity_m_per_s.y*config_.dt_s;particle.position_m=wrap_position(particle.position_m);}time_s_+=config_.dt_s;}deposit_and_solve();}



ElectromagneticPic2D::ElectromagneticPic2D(ElectromagneticPic2DConfig config):config_(config){
    if(config_.nx<2U||config_.ny<2U||!(config_.length_x_m>0.0)||!(config_.length_y_m>0.0)
       ||!(config_.relative_permittivity>0.0)||!(config_.relative_permeability>0.0)||!(config_.dt_s>0.0)||!std::isfinite(config_.dt_s))
        throw std::invalid_argument("invalid electromagnetic PIC 2-D configuration");
    const std::size_t n=config_.nx*config_.ny;
    charge_density_.assign(n,0.0);current_x_.assign(n,0.0);current_y_.assign(n,0.0);current_z_.assign(n,0.0);
    electric_x_.assign(n,0.0);electric_y_.assign(n,0.0);electric_z_.assign(n,0.0);
    magnetic_x_.assign(n,0.0);magnetic_y_.assign(n,0.0);magnetic_z_.assign(n,0.0);
}

double ElectromagneticPic2D::epsilon() const noexcept {return epsilon0*config_.relative_permittivity;}
double ElectromagneticPic2D::mu() const noexcept {return mu0*config_.relative_permeability;}
double ElectromagneticPic2D::medium_light_speed() const noexcept {return 1.0/std::sqrt(epsilon()*mu());}

Vec2 ElectromagneticPic2D::wrap_position(Vec2 position_m) const {position_m.x=wrap_periodic_value(position_m.x,config_.length_x_m);position_m.y=wrap_periodic_value(position_m.y,config_.length_y_m);return position_m;}

void ElectromagneticPic2D::set_particles(std::vector<PicParticle2D> particles){
    for(auto& particle:particles){
        if(!(particle.mass_kg>0.0)||!std::isfinite(particle.mass_kg)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))
            throw std::invalid_argument("invalid electromagnetic PIC 2-D particle");
        particle.position_m=wrap_position(particle.position_m);
    }
    particles_=std::move(particles);
}

void ElectromagneticPic2D::set_fields(std::span<const double> electric_x_v_per_m,
                                      std::span<const double> electric_y_v_per_m,
                                      std::span<const double> electric_z_v_per_m,
                                      std::span<const double> magnetic_x_t,
                                      std::span<const double> magnetic_y_t,
                                      std::span<const double> magnetic_z_t){
    const std::size_t n=config_.nx*config_.ny;
    if(electric_x_v_per_m.size()!=n||electric_y_v_per_m.size()!=n||electric_z_v_per_m.size()!=n||magnetic_x_t.size()!=n||magnetic_y_t.size()!=n||magnetic_z_t.size()!=n)
        throw std::invalid_argument("2-D EM PIC field arrays must match grid size");
    electric_x_.assign(electric_x_v_per_m.begin(),electric_x_v_per_m.end());electric_y_.assign(electric_y_v_per_m.begin(),electric_y_v_per_m.end());electric_z_.assign(electric_z_v_per_m.begin(),electric_z_v_per_m.end());
    magnetic_x_.assign(magnetic_x_t.begin(),magnetic_x_t.end());magnetic_y_.assign(magnetic_y_t.begin(),magnetic_y_t.end());magnetic_z_.assign(magnetic_z_t.begin(),magnetic_z_t.end());
}

void ElectromagneticPic2D::initialize_tm_z_mode(double electric_z_amplitude_v_per_m,std::size_t mode_x,std::size_t mode_y){
    if(mode_x==0U&&mode_y==0U)throw std::invalid_argument("2-D EM PIC mode index must be non-zero");
    const double kx=2.0*std::numbers::pi*static_cast<double>(mode_x)/config_.length_x_m;
    const double ky=2.0*std::numbers::pi*static_cast<double>(mode_y)/config_.length_y_m;
    const double k=std::hypot(kx,ky);const double omega=medium_light_speed()*k;
    for(std::size_t iy=0;iy<config_.ny;++iy){
        const double y=static_cast<double>(iy)*dy();
        for(std::size_t ix=0;ix<config_.nx;++ix){
            const double x=static_cast<double>(ix)*dx();const double phase=kx*x+ky*y;const std::size_t idx=flat_index_2d(ix,iy,config_.nx);
            electric_x_[idx]=0.0;electric_y_[idx]=0.0;electric_z_[idx]=electric_z_amplitude_v_per_m*std::sin(phase);
            // Right-travelling TM_z plane wave: B=(1/omega) k x E.
            magnetic_x_[idx]=(ky/omega)*electric_z_[idx];magnetic_y_[idx]=-(kx/omega)*electric_z_[idx];magnetic_z_[idx]=0.0;
        }
    }
}

double ElectromagneticPic2D::bilinear(std::span<const double> nodal,Vec2 position_m) const {
    const Vec2 wrapped=wrap_position(position_m);const double ux=wrapped.x/dx();const double uy=wrapped.y/dy();
    const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;
    const auto ix1=(ix0+1U)%config_.nx;const auto iy1=(iy0+1U)%config_.ny;const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);
    return nodal[flat_index_2d(ix0,iy0,config_.nx)]*(1.0-fx)*(1.0-fy)+nodal[flat_index_2d(ix1,iy0,config_.nx)]*fx*(1.0-fy)+nodal[flat_index_2d(ix0,iy1,config_.nx)]*(1.0-fx)*fy+nodal[flat_index_2d(ix1,iy1,config_.nx)]*fx*fy;
}

double ElectromagneticPic2D::ddx(std::span<const double> values,std::size_t ix,std::size_t iy) const {
    const std::size_t left=(ix+config_.nx-1U)%config_.nx,right=(ix+1U)%config_.nx;return (values[flat_index_2d(right,iy,config_.nx)]-values[flat_index_2d(left,iy,config_.nx)])/(2.0*dx());
}

double ElectromagneticPic2D::ddy(std::span<const double> values,std::size_t ix,std::size_t iy) const {
    const std::size_t down=(iy+config_.ny-1U)%config_.ny,up=(iy+1U)%config_.ny;return (values[flat_index_2d(ix,up,config_.nx)]-values[flat_index_2d(ix,down,config_.nx)])/(2.0*dy());
}

void ElectromagneticPic2D::deposit_sources(){
    charge_density_=deposit_cic_charge_density_2d(particles_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m);
    std::fill(current_x_.begin(),current_x_.end(),0.0);std::fill(current_y_.begin(),current_y_.end(),0.0);std::fill(current_z_.begin(),current_z_.end(),0.0);
    if(config_.solve_longitudinal_poisson){const auto field=periodic_electric_field_from_charge_density_2d(charge_density_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m,config_.relative_permittivity,config_.neutralize_mean_charge);electric_x_=field.electric_x_v_per_m;electric_y_=field.electric_y_v_per_m;}
}

ElectromagneticField ElectromagneticPic2D::gather_field(Vec2 position_m) const {
    return {{bilinear(electric_x_,position_m),bilinear(electric_y_,position_m),bilinear(electric_z_,position_m)},
            {bilinear(magnetic_x_,position_m),bilinear(magnetic_y_,position_m),bilinear(magnetic_z_,position_m)}};
}

void ElectromagneticPic2D::advance_magnetic(double dt_s){
    std::vector<double> next_x(magnetic_x_),next_y(magnetic_y_),next_z(magnetic_z_);
    for(std::size_t iy=0;iy<config_.ny;++iy){
        for(std::size_t ix=0;ix<config_.nx;++ix){
            const std::size_t idx=flat_index_2d(ix,iy,config_.nx);
            next_x[idx]-=dt_s*ddy(electric_z_,ix,iy);
            next_y[idx]+=dt_s*ddx(electric_z_,ix,iy);
            next_z[idx]-=dt_s*(ddx(electric_y_,ix,iy)-ddy(electric_x_,ix,iy));
        }
    }
    magnetic_x_.swap(next_x);magnetic_y_.swap(next_y);magnetic_z_.swap(next_z);
}

void ElectromagneticPic2D::deposit_transverse_current_z(std::span<const PicParticle2D> particles_old){
    std::fill(current_z_.begin(),current_z_.end(),0.0);const double cell_area=dx()*dy();
    for(std::size_t particle_index=0;particle_index<particles_.size();++particle_index){
        Vec2 old_position=particles_old[particle_index].position_m;Vec2 new_position=particles_[particle_index].position_m;
        double delta_x=new_position.x-old_position.x;if(delta_x>0.5*config_.length_x_m)old_position.x+=config_.length_x_m;else if(delta_x<-0.5*config_.length_x_m)new_position.x+=config_.length_x_m;
        double delta_y=new_position.y-old_position.y;if(delta_y>0.5*config_.length_y_m)old_position.y+=config_.length_y_m;else if(delta_y<-0.5*config_.length_y_m)new_position.y+=config_.length_y_m;
        const Vec2 mid=wrap_position({0.5*(old_position.x+new_position.x),0.5*(old_position.y+new_position.y)});
        const double ux=mid.x/dx();const double uy=mid.y/dy();const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;const auto ix1=(ix0+1U)%config_.nx;const auto iy1=(iy0+1U)%config_.ny;const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);
        const double vz=0.5*(particles_old[particle_index].velocity_m_per_s.z+particles_[particle_index].velocity_m_per_s.z);
        const double current=particles_[particle_index].charge_c*particles_[particle_index].weight*vz/cell_area;
        current_z_[flat_index_2d(ix0,iy0,config_.nx)]+=current*(1.0-fx)*(1.0-fy);
        current_z_[flat_index_2d(ix1,iy0,config_.nx)]+=current*fx*(1.0-fy);
        current_z_[flat_index_2d(ix0,iy1,config_.nx)]+=current*(1.0-fx)*fy;
        current_z_[flat_index_2d(ix1,iy1,config_.nx)]+=current*fx*fy;
    }
}

void ElectromagneticPic2D::step(std::size_t steps,const NeutralCollisionModel* collision_model){
    for(std::size_t step_index=0;step_index<steps;++step_index){
        const auto particles_old=particles_;advance_magnetic(0.5*config_.dt_s);
        for(auto& particle:particles_){
            ChargedParticle proxy;proxy.position_m={particle.position_m.x,particle.position_m.y,0.0};proxy.velocity_m_per_s=particle.velocity_m_per_s;proxy.charge_c=particle.charge_c;proxy.mass_kg=particle.mass_kg;proxy.weight=particle.weight;
            vay_push(proxy,gather_field(particle.position_m),config_.dt_s);
            particle.velocity_m_per_s=proxy.velocity_m_per_s;particle.position_m=wrap_position({proxy.position_m.x,proxy.position_m.y});
        }
        const auto cc=deposit_charge_conserving_current_2d(particles_old,particles_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m,config_.dt_s);charge_density_=cc.charge_density_new_c_per_m3;current_x_=cc.current_x_a_per_m2;current_y_=cc.current_y_a_per_m2;last_continuity_residual_=cc.continuity_linf_residual;deposit_transverse_current_z(particles_old);
        std::vector<double> next_ex(electric_x_),next_ey(electric_y_),next_ez(electric_z_);const double eps=epsilon();const double permeability=mu();
        for(std::size_t iy=0;iy<config_.ny;++iy){
            for(std::size_t ix=0;ix<config_.nx;++ix){
                const std::size_t idx=flat_index_2d(ix,iy,config_.nx);
                next_ex[idx]+=config_.dt_s*((ddy(magnetic_z_,ix,iy)/permeability)-current_x_[idx])/eps;
                next_ey[idx]+=config_.dt_s*((-ddx(magnetic_z_,ix,iy)/permeability)-current_y_[idx])/eps;
                next_ez[idx]+=config_.dt_s*(((ddx(magnetic_y_,ix,iy)-ddy(magnetic_x_,ix,iy))/permeability)-current_z_[idx])/eps;
            }
        }
        electric_x_.swap(next_ex);electric_y_.swap(next_ey);electric_z_.swap(next_ez);advance_magnetic(0.5*config_.dt_s);
        if(config_.solve_longitudinal_poisson){const auto field=periodic_electric_field_from_charge_density_2d(charge_density_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m,config_.relative_permittivity,config_.neutralize_mean_charge);electric_x_=field.electric_x_v_per_m;electric_y_=field.electric_y_v_per_m;}
        if(collision_model){
            std::vector<ChargedParticle> proxy;proxy.reserve(particles_.size());for(const auto& particle:particles_){ChargedParticle value;value.position_m={particle.position_m.x,particle.position_m.y,0.0};value.velocity_m_per_s=particle.velocity_m_per_s;value.charge_c=particle.charge_c;value.mass_kg=particle.mass_kg;value.weight=particle.weight;proxy.push_back(value);}const auto stats=apply_monte_carlo_collisions(proxy,config_.dt_s,*collision_model);accumulated_collisions_.elastic_events+=stats.elastic_events;accumulated_collisions_.ionization_events+=stats.ionization_events;accumulated_collisions_.energy_loss_j+=stats.energy_loss_j;particles_.clear();particles_.reserve(proxy.size());for(const auto& value:proxy){PicParticle2D particle;particle.position_m=wrap_position({value.position_m.x,value.position_m.y});particle.velocity_m_per_s=value.velocity_m_per_s;particle.charge_c=value.charge_c;particle.mass_kg=value.mass_kg;particle.weight=value.weight;particles_.push_back(particle);}
        }
        time_s_+=config_.dt_s;
    }
}

ElectromagneticPic2DDiagnostics ElectromagneticPic2D::diagnostics() const {
    ElectromagneticPic2DDiagnostics result;const double cell_area=dx()*dy();const double eps=epsilon();const double permeability=mu();
    for(std::size_t idx=0;idx<electric_x_.size();++idx){const double e2=electric_x_[idx]*electric_x_[idx]+electric_y_[idx]*electric_y_[idx]+electric_z_[idx]*electric_z_[idx];const double b2=magnetic_x_[idx]*magnetic_x_[idx]+magnetic_y_[idx]*magnetic_y_[idx]+magnetic_z_[idx]*magnetic_z_[idx];result.field_energy_j+=0.5*eps*e2*cell_area+0.5*b2*cell_area/permeability;}
    for(const auto& particle:particles_){ChargedParticle proxy;proxy.velocity_m_per_s=particle.velocity_m_per_s;proxy.mass_kg=particle.mass_kg;proxy.charge_c=particle.charge_c;proxy.weight=particle.weight;result.particle_kinetic_energy_j+=particle.weight*kinetic_energy_j(proxy);}result.total_energy_j=result.field_energy_j+result.particle_kinetic_energy_j;result.charge_continuity_linf_residual=last_continuity_residual_;result.collision_statistics=accumulated_collisions_;result.particle_count=particles_.size();return result;
}

ParticleBoundary2DResult apply_particle_box_boundary_2d(PicParticle2D& particle,const PicBoundaryBox2D& box,ParticleWallMode mode,const SecondaryEmissionModel* secondary){
    if(!(box.maximum_m.x>box.minimum_m.x&&box.maximum_m.y>box.minimum_m.y))throw std::invalid_argument("invalid 2-D particle boundary box");
    ParticleBoundary2DResult result;auto outside=[](double value,double lo,double hi){return value<lo||value>hi;};
    if(!outside(particle.position_m.x,box.minimum_m.x,box.maximum_m.x)&&!outside(particle.position_m.y,box.minimum_m.y,box.maximum_m.y))return result;
    result.impacted=true;ChargedParticle proxy;proxy.mass_kg=particle.mass_kg;proxy.charge_c=particle.charge_c;proxy.weight=particle.weight;proxy.velocity_m_per_s=particle.velocity_m_per_s;const double incident_ev=kinetic_energy_j(proxy)/elementary_charge;
    auto handle_axis=[&](double& position,double& velocity,double lo,double hi,Vec2 low_normal,Vec2 high_normal){if(position>=lo&&position<=hi)return;if(position<lo)result.impact_normal=low_normal;else result.impact_normal=high_normal;if(mode==ParticleWallMode::absorb){result.alive=false;return;}while(position<lo||position>hi){if(position<lo){position=2.0*lo-position;velocity=-velocity;}else{position=2.0*hi-position;velocity=-velocity;}}};
    handle_axis(particle.position_m.x,particle.velocity_m_per_s.x,box.minimum_m.x,box.maximum_m.x,{-1.0,0.0},{1.0,0.0});if(result.alive)handle_axis(particle.position_m.y,particle.velocity_m_per_s.y,box.minimum_m.y,box.maximum_m.y,{0.0,-1.0},{0.0,1.0});if(!result.alive&&secondary)result.secondary_macro_weight=particle.weight*secondary_electron_yield(incident_ev,*secondary);return result;
}

void apply_electrostatic_field_boundary_2d(std::vector<double>& electric_x_v_per_m,std::vector<double>& electric_y_v_per_m,
                                           std::size_t nx,std::size_t ny,const GridBoundary2DConfig& boundary){
    if(electric_x_v_per_m.size()!=nx*ny||electric_y_v_per_m.size()!=nx*ny||nx<2U||ny<2U)throw std::invalid_argument("invalid 2-D field boundary arrays");
    if(boundary.mode==GridBoundaryMode2D::periodic)return;
    if(boundary.mode==GridBoundaryMode2D::electric_wall){
        for(std::size_t iy=0;iy<ny;++iy){electric_y_v_per_m[flat_index_2d(0U,iy,nx)]=0.0;electric_y_v_per_m[flat_index_2d(nx-1U,iy,nx)]=0.0;}
        for(std::size_t ix=0;ix<nx;++ix){electric_x_v_per_m[flat_index_2d(ix,0U,nx)]=0.0;electric_x_v_per_m[flat_index_2d(ix,ny-1U,nx)]=0.0;}
        return;
    }
    if(boundary.mode==GridBoundaryMode2D::absorbing_sponge){
        const std::size_t layers=std::min(boundary.sponge_cells,std::min(nx,ny)/2U);if(layers==0U)return;if(!(boundary.sponge_strength>=0.0))throw std::invalid_argument("invalid absorbing sponge strength");
        for(std::size_t iy=0;iy<ny;++iy){
            for(std::size_t ix=0;ix<nx;++ix){
                const std::size_t edge=std::min(std::min(ix,nx-1U-ix),std::min(iy,ny-1U-iy));
                if(edge>=layers){continue;}
                const double normalized=1.0-static_cast<double>(edge)/static_cast<double>(layers);const double damping=std::exp(-boundary.sponge_strength*normalized*normalized);const std::size_t idx=flat_index_2d(ix,iy,nx);electric_x_v_per_m[idx]*=damping;electric_y_v_per_m[idx]*=damping;
            }
        }
    }
}


double StaggeredElectromagneticPic2D::epsilon() const noexcept {return epsilon0*config_.relative_permittivity;}
double StaggeredElectromagneticPic2D::mu() const noexcept {return mu0*config_.relative_permeability;}
double StaggeredElectromagneticPic2D::medium_light_speed() const noexcept {return 1.0/std::sqrt(epsilon()*mu());}
bool StaggeredElectromagneticPic2D::periodic_fields() const noexcept {return config_.field_boundary.mode==GridBoundaryMode2D::periodic;}

StaggeredElectromagneticPic2D::StaggeredElectromagneticPic2D(StaggeredElectromagneticPic2DConfig config):config_(config){
    if(config_.nx<2U||config_.ny<2U||!(config_.length_x_m>0.0)||!(config_.length_y_m>0.0)
       ||!(config_.relative_permittivity>0.0)||!(config_.relative_permeability>0.0)||!(config_.dt_s>0.0)||!std::isfinite(config_.dt_s))
        throw std::invalid_argument("invalid staggered 2-D EM PIC configuration");
    const double cfl=config_.dt_s*medium_light_speed()*std::sqrt(1.0/(dx()*dx())+1.0/(dy()*dy()));
    if(cfl>1.0)throw std::invalid_argument("staggered 2-D EM PIC CFL must be <= 1");
    const std::size_t n=config_.nx*config_.ny;
    charge_density_.assign(n,0.0);current_x_.assign(n,0.0);current_y_.assign(n,0.0);current_z_.assign(n,0.0);
    electric_x_.assign(n,0.0);electric_y_.assign(n,0.0);electric_z_.assign(n,0.0);
    magnetic_x_.assign(n,0.0);magnetic_y_.assign(n,0.0);magnetic_z_.assign(n,0.0);
}

Vec2 StaggeredElectromagneticPic2D::place_particle(Vec2 position_m) const {
    if(config_.periodic_particles){position_m.x=wrap_periodic_value(position_m.x,config_.length_x_m);position_m.y=wrap_periodic_value(position_m.y,config_.length_y_m);}return position_m;
}

void StaggeredElectromagneticPic2D::set_particles(std::vector<PicParticle2D> particles){
    for(auto& particle:particles){
        if(!(particle.mass_kg>0.0)||!std::isfinite(particle.mass_kg)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))
            throw std::invalid_argument("invalid staggered 2-D EM PIC particle");
        particle.position_m=place_particle(particle.position_m);
    }
    particles_=std::move(particles);
}

void StaggeredElectromagneticPic2D::set_fields(std::span<const double> electric_x_v_per_m,
                                                std::span<const double> electric_y_v_per_m,
                                                std::span<const double> electric_z_v_per_m,
                                                std::span<const double> magnetic_x_t,
                                                std::span<const double> magnetic_y_t,
                                                std::span<const double> magnetic_z_t){
    const std::size_t n=config_.nx*config_.ny;
    if(electric_x_v_per_m.size()!=n||electric_y_v_per_m.size()!=n||electric_z_v_per_m.size()!=n||magnetic_x_t.size()!=n||magnetic_y_t.size()!=n||magnetic_z_t.size()!=n)
        throw std::invalid_argument("staggered 2-D EM PIC field arrays must match grid size");
    electric_x_.assign(electric_x_v_per_m.begin(),electric_x_v_per_m.end());electric_y_.assign(electric_y_v_per_m.begin(),electric_y_v_per_m.end());electric_z_.assign(electric_z_v_per_m.begin(),electric_z_v_per_m.end());
    magnetic_x_.assign(magnetic_x_t.begin(),magnetic_x_t.end());magnetic_y_.assign(magnetic_y_t.begin(),magnetic_y_t.end());magnetic_z_.assign(magnetic_z_t.begin(),magnetic_z_t.end());
    apply_field_boundary();
}

void StaggeredElectromagneticPic2D::initialize_tm_z_mode(double electric_z_amplitude_v_per_m,std::size_t mode_x,std::size_t mode_y){
    if(mode_x==0U&&mode_y==0U)throw std::invalid_argument("staggered 2-D EM PIC mode index must be non-zero");
    const double kx=2.0*std::numbers::pi*static_cast<double>(mode_x)/config_.length_x_m;
    const double ky=2.0*std::numbers::pi*static_cast<double>(mode_y)/config_.length_y_m;
    const double k=std::hypot(kx,ky);const double omega=medium_light_speed()*k;
    for(std::size_t iy=0;iy<config_.ny;++iy){
        for(std::size_t ix=0;ix<config_.nx;++ix){
            const std::size_t idx=flat_index_2d(ix,iy,config_.nx);
            const double x_ez=static_cast<double>(ix)*dx();const double y_ez=static_cast<double>(iy)*dy();
            const double phase_ez=kx*x_ez+ky*y_ez;electric_x_[idx]=0.0;electric_y_[idx]=0.0;electric_z_[idx]=electric_z_amplitude_v_per_m*std::sin(phase_ez);
            const double phase_bx=kx*static_cast<double>(ix)*dx()+ky*(static_cast<double>(iy)+0.5)*dy();
            const double phase_by=kx*(static_cast<double>(ix)+0.5)*dx()+ky*static_cast<double>(iy)*dy();
            magnetic_x_[idx]=(ky/omega)*electric_z_amplitude_v_per_m*std::sin(phase_bx);
            magnetic_y_[idx]=-(kx/omega)*electric_z_amplitude_v_per_m*std::sin(phase_by);
            magnetic_z_[idx]=0.0;
        }
    }
    apply_field_boundary();
}

static std::size_t bounded_index_from_double(double value,std::size_t size){
    if(value<=0.0)return 0U;
    const double upper=static_cast<double>(size-1U);
    if(value>=upper)return size-1U;
    return static_cast<std::size_t>(std::floor(value));
}

double StaggeredElectromagneticPic2D::sample(std::span<const double> values,Vec2 position_m,double offset_x,double offset_y) const {
    if(values.size()!=config_.nx*config_.ny)throw std::invalid_argument("invalid staggered sample array");
    double ux{},uy{};
    if(periodic_fields()){
        position_m.x=wrap_periodic_value(position_m.x,config_.length_x_m);position_m.y=wrap_periodic_value(position_m.y,config_.length_y_m);
        ux=position_m.x/dx()-offset_x;uy=position_m.y/dy()-offset_y;
        ux-=std::floor(ux/static_cast<double>(config_.nx))*static_cast<double>(config_.nx);
        uy-=std::floor(uy/static_cast<double>(config_.ny))*static_cast<double>(config_.ny);
        const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;
        const auto ix1=(ix0+1U)%config_.nx;const auto iy1=(iy0+1U)%config_.ny;const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);
        return values[flat_index_2d(ix0,iy0,config_.nx)]*(1.0-fx)*(1.0-fy)+values[flat_index_2d(ix1,iy0,config_.nx)]*fx*(1.0-fy)+values[flat_index_2d(ix0,iy1,config_.nx)]*(1.0-fx)*fy+values[flat_index_2d(ix1,iy1,config_.nx)]*fx*fy;
    }
    position_m.x=std::clamp(position_m.x,0.0,config_.length_x_m);position_m.y=std::clamp(position_m.y,0.0,config_.length_y_m);
    ux=std::clamp(position_m.x/dx()-offset_x,0.0,static_cast<double>(config_.nx-1U));
    uy=std::clamp(position_m.y/dy()-offset_y,0.0,static_cast<double>(config_.ny-1U));
    const auto ix0=bounded_index_from_double(ux,config_.nx);const auto iy0=bounded_index_from_double(uy,config_.ny);
    const auto ix1=std::min(ix0+1U,config_.nx-1U);const auto iy1=std::min(iy0+1U,config_.ny-1U);const double fx=std::clamp(ux-static_cast<double>(ix0),0.0,1.0);const double fy=std::clamp(uy-static_cast<double>(iy0),0.0,1.0);
    return values[flat_index_2d(ix0,iy0,config_.nx)]*(1.0-fx)*(1.0-fy)+values[flat_index_2d(ix1,iy0,config_.nx)]*fx*(1.0-fy)+values[flat_index_2d(ix0,iy1,config_.nx)]*(1.0-fx)*fy+values[flat_index_2d(ix1,iy1,config_.nx)]*fx*fy;
}

std::size_t left_or_same(std::size_t index){return index==0U?0U:index-1U;}
std::size_t right_or_same(std::size_t index,std::size_t size){return index+1U<size?index+1U:index;}

double StaggeredElectromagneticPic2D::diff_x_forward(std::span<const double> values,std::size_t ix,std::size_t iy) const {const std::size_t r=periodic_fields()?((ix+1U)%config_.nx):right_or_same(ix,config_.nx);return (values[flat_index_2d(r,iy,config_.nx)]-values[flat_index_2d(ix,iy,config_.nx)])/dx();}
double StaggeredElectromagneticPic2D::diff_y_forward(std::span<const double> values,std::size_t ix,std::size_t iy) const {const std::size_t u=periodic_fields()?((iy+1U)%config_.ny):right_or_same(iy,config_.ny);return (values[flat_index_2d(ix,u,config_.nx)]-values[flat_index_2d(ix,iy,config_.nx)])/dy();}
double StaggeredElectromagneticPic2D::diff_x_backward(std::span<const double> values,std::size_t ix,std::size_t iy) const {const std::size_t l=periodic_fields()?((ix+config_.nx-1U)%config_.nx):left_or_same(ix);return (values[flat_index_2d(ix,iy,config_.nx)]-values[flat_index_2d(l,iy,config_.nx)])/dx();}
double StaggeredElectromagneticPic2D::diff_y_backward(std::span<const double> values,std::size_t ix,std::size_t iy) const {const std::size_t d=periodic_fields()?((iy+config_.ny-1U)%config_.ny):left_or_same(iy);return (values[flat_index_2d(ix,iy,config_.nx)]-values[flat_index_2d(ix,d,config_.nx)])/dy();}

ElectromagneticField StaggeredElectromagneticPic2D::gather_field(Vec2 position_m) const {
    return {{sample(electric_x_,position_m,0.5,0.0),sample(electric_y_,position_m,0.0,0.5),sample(electric_z_,position_m,0.0,0.0)},
            {sample(magnetic_x_,position_m,0.0,0.5),sample(magnetic_y_,position_m,0.5,0.0),sample(magnetic_z_,position_m,0.5,0.5)}};
}

void StaggeredElectromagneticPic2D::advance_magnetic(double dt_s){
    std::vector<double> next_x(magnetic_x_),next_y(magnetic_y_),next_z(magnetic_z_);
    for(std::size_t iy=0;iy<config_.ny;++iy){
        for(std::size_t ix=0;ix<config_.nx;++ix){
            const std::size_t idx=flat_index_2d(ix,iy,config_.nx);
            next_x[idx]-=dt_s*diff_y_forward(electric_z_,ix,iy);
            next_y[idx]+=dt_s*diff_x_forward(electric_z_,ix,iy);
            next_z[idx]-=dt_s*(diff_x_forward(electric_y_,ix,iy)-diff_y_forward(electric_x_,ix,iy));
        }
    }
    magnetic_x_.swap(next_x);magnetic_y_.swap(next_y);magnetic_z_.swap(next_z);apply_field_boundary();
}

void StaggeredElectromagneticPic2D::deposit_transverse_current_z(std::span<const PicParticle2D> particles_old){
    std::fill(current_z_.begin(),current_z_.end(),0.0);const double cell_area=dx()*dy();
    for(std::size_t particle_index=0;particle_index<particles_.size();++particle_index){
        Vec2 old_position=particles_old[particle_index].position_m;Vec2 new_position=particles_[particle_index].position_m;
        if(config_.periodic_particles){double delta_x=new_position.x-old_position.x;if(delta_x>0.5*config_.length_x_m)old_position.x+=config_.length_x_m;else if(delta_x<-0.5*config_.length_x_m)new_position.x+=config_.length_x_m;double delta_y=new_position.y-old_position.y;if(delta_y>0.5*config_.length_y_m)old_position.y+=config_.length_y_m;else if(delta_y<-0.5*config_.length_y_m)new_position.y+=config_.length_y_m;}
        const Vec2 mid=place_particle({0.5*(old_position.x+new_position.x),0.5*(old_position.y+new_position.y)});
        const double ux=std::clamp(mid.x/dx(),0.0,static_cast<double>(config_.nx-1U));const double uy=std::clamp(mid.y/dy(),0.0,static_cast<double>(config_.ny-1U));
        const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;const auto ix1=periodic_fields()?((ix0+1U)%config_.nx):std::min(ix0+1U,config_.nx-1U);const auto iy1=periodic_fields()?((iy0+1U)%config_.ny):std::min(iy0+1U,config_.ny-1U);const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);
        const double vz=0.5*(particles_old[particle_index].velocity_m_per_s.z+particles_[particle_index].velocity_m_per_s.z);const double current=particles_[particle_index].charge_c*particles_[particle_index].weight*vz/cell_area;
        current_z_[flat_index_2d(ix0,iy0,config_.nx)]+=current*(1.0-fx)*(1.0-fy);current_z_[flat_index_2d(ix1,iy0,config_.nx)]+=current*fx*(1.0-fy);current_z_[flat_index_2d(ix0,iy1,config_.nx)]+=current*(1.0-fx)*fy;current_z_[flat_index_2d(ix1,iy1,config_.nx)]+=current*fx*fy;
    }
}

void StaggeredElectromagneticPic2D::update_electric(double dt_s){
    std::vector<double> next_ex(electric_x_),next_ey(electric_y_),next_ez(electric_z_);const double eps=epsilon();const double permeability=mu();
    for(std::size_t iy=0;iy<config_.ny;++iy){
        for(std::size_t ix=0;ix<config_.nx;++ix){
            const std::size_t idx=flat_index_2d(ix,iy,config_.nx);
            next_ex[idx]+=dt_s*((diff_y_backward(magnetic_z_,ix,iy)/permeability)-current_x_[idx])/eps;
            next_ey[idx]+=dt_s*((-diff_x_backward(magnetic_z_,ix,iy)/permeability)-current_y_[idx])/eps;
            next_ez[idx]+=dt_s*(((diff_x_backward(magnetic_y_,ix,iy)-diff_y_backward(magnetic_x_,ix,iy))/permeability)-current_z_[idx])/eps;
        }
    }
    electric_x_.swap(next_ex);electric_y_.swap(next_ey);electric_z_.swap(next_ez);apply_field_boundary();
}

void StaggeredElectromagneticPic2D::apply_field_boundary(){
    if(config_.field_boundary.mode==GridBoundaryMode2D::periodic)return;
    if(config_.field_boundary.mode==GridBoundaryMode2D::electric_wall){
        for(std::size_t iy=0;iy<config_.ny;++iy){
            electric_y_[flat_index_2d(0U,iy,config_.nx)]=0.0;electric_z_[flat_index_2d(0U,iy,config_.nx)]=0.0;
            electric_y_[flat_index_2d(config_.nx-1U,iy,config_.nx)]=0.0;electric_z_[flat_index_2d(config_.nx-1U,iy,config_.nx)]=0.0;
            magnetic_x_[flat_index_2d(0U,iy,config_.nx)]=0.0;magnetic_x_[flat_index_2d(config_.nx-1U,iy,config_.nx)]=0.0;
        }
        for(std::size_t ix=0;ix<config_.nx;++ix){
            electric_x_[flat_index_2d(ix,0U,config_.nx)]=0.0;electric_z_[flat_index_2d(ix,0U,config_.nx)]=0.0;
            electric_x_[flat_index_2d(ix,config_.ny-1U,config_.nx)]=0.0;electric_z_[flat_index_2d(ix,config_.ny-1U,config_.nx)]=0.0;
            magnetic_y_[flat_index_2d(ix,0U,config_.nx)]=0.0;magnetic_y_[flat_index_2d(ix,config_.ny-1U,config_.nx)]=0.0;
        }
        return;
    }
    if(config_.field_boundary.mode==GridBoundaryMode2D::absorbing_sponge){
        const std::size_t layers=std::min(config_.field_boundary.sponge_cells,std::min(config_.nx,config_.ny)/2U);if(layers==0U)return;if(!(config_.field_boundary.sponge_strength>=0.0))throw std::invalid_argument("invalid staggered sponge strength");
        for(std::size_t iy=0;iy<config_.ny;++iy){
            for(std::size_t ix=0;ix<config_.nx;++ix){
                const std::size_t edge=std::min(std::min(ix,config_.nx-1U-ix),std::min(iy,config_.ny-1U-iy));if(edge>=layers)continue;const double normalized=1.0-static_cast<double>(edge)/static_cast<double>(layers);const double damping=std::exp(-config_.field_boundary.sponge_strength*normalized*normalized);const std::size_t idx=flat_index_2d(ix,iy,config_.nx);
                electric_x_[idx]*=damping;electric_y_[idx]*=damping;electric_z_[idx]*=damping;magnetic_x_[idx]*=damping;magnetic_y_[idx]*=damping;magnetic_z_[idx]*=damping;
            }
        }
    }
}

void StaggeredElectromagneticPic2D::apply_particle_boundary(){
    if(config_.periodic_particles){for(auto& particle:particles_)particle.position_m=place_particle(particle.position_m);return;}
    PicBoundaryBox2D box;box.minimum_m={0.0,0.0};box.maximum_m={config_.length_x_m,config_.length_y_m};
    std::vector<PicParticle2D> kept;kept.reserve(particles_.size());
    for(auto particle:particles_){
        const auto result=apply_particle_box_boundary_2d(particle,box,config_.particle_boundary,config_.enable_secondary_yield_report?&config_.secondary_model:nullptr);
        if(result.impacted)secondary_macro_weight_+=result.secondary_macro_weight;
        if(result.alive||!config_.remove_absorbed_particles)kept.push_back(particle);else ++absorbed_particles_;
    }
    particles_.swap(kept);
}

void StaggeredElectromagneticPic2D::step(std::size_t steps,const NeutralCollisionModel* collision_model){
    for(std::size_t step_index=0;step_index<steps;++step_index){
        const auto particles_old=particles_;advance_magnetic(0.5*config_.dt_s);
        for(auto& particle:particles_){ChargedParticle proxy;proxy.position_m={particle.position_m.x,particle.position_m.y,0.0};proxy.velocity_m_per_s=particle.velocity_m_per_s;proxy.charge_c=particle.charge_c;proxy.mass_kg=particle.mass_kg;proxy.weight=particle.weight;vay_push(proxy,gather_field(particle.position_m),config_.dt_s);particle.velocity_m_per_s=proxy.velocity_m_per_s;particle.position_m={proxy.position_m.x,proxy.position_m.y};}
        apply_particle_boundary();
        if(particles_.size()==particles_old.size()){const auto cc=deposit_charge_conserving_current_2d(particles_old,particles_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m,config_.dt_s);charge_density_=cc.charge_density_new_c_per_m3;current_x_=cc.current_x_a_per_m2;current_y_=cc.current_y_a_per_m2;last_continuity_residual_=cc.continuity_linf_residual;deposit_transverse_current_z(particles_old);}else{charge_density_=deposit_cic_charge_density_2d(particles_,config_.nx,config_.ny,config_.length_x_m,config_.length_y_m);std::fill(current_x_.begin(),current_x_.end(),0.0);std::fill(current_y_.begin(),current_y_.end(),0.0);std::fill(current_z_.begin(),current_z_.end(),0.0);last_continuity_residual_=0.0;}
        update_electric(config_.dt_s);advance_magnetic(0.5*config_.dt_s);
        if(collision_model){std::vector<ChargedParticle> proxy;proxy.reserve(particles_.size());for(const auto& particle:particles_){ChargedParticle value;value.position_m={particle.position_m.x,particle.position_m.y,0.0};value.velocity_m_per_s=particle.velocity_m_per_s;value.charge_c=particle.charge_c;value.mass_kg=particle.mass_kg;value.weight=particle.weight;proxy.push_back(value);}const auto stats=apply_monte_carlo_collisions(proxy,config_.dt_s,*collision_model);accumulated_collisions_.elastic_events+=stats.elastic_events;accumulated_collisions_.ionization_events+=stats.ionization_events;accumulated_collisions_.energy_loss_j+=stats.energy_loss_j;particles_.clear();particles_.reserve(proxy.size());for(const auto& value:proxy){PicParticle2D particle;particle.position_m=place_particle({value.position_m.x,value.position_m.y});particle.velocity_m_per_s=value.velocity_m_per_s;particle.charge_c=value.charge_c;particle.mass_kg=value.mass_kg;particle.weight=value.weight;particles_.push_back(particle);}apply_particle_boundary();}
        time_s_+=config_.dt_s;
    }
}

StaggeredElectromagneticPic2DDiagnostics StaggeredElectromagneticPic2D::diagnostics() const {
    StaggeredElectromagneticPic2DDiagnostics result;const double cell_area=dx()*dy();const double eps=epsilon();const double permeability=mu();
    for(std::size_t idx=0;idx<electric_x_.size();++idx){const double e2=electric_x_[idx]*electric_x_[idx]+electric_y_[idx]*electric_y_[idx]+electric_z_[idx]*electric_z_[idx];const double b2=magnetic_x_[idx]*magnetic_x_[idx]+magnetic_y_[idx]*magnetic_y_[idx]+magnetic_z_[idx]*magnetic_z_[idx];result.field_energy_j+=0.5*eps*e2*cell_area+0.5*b2*cell_area/permeability;}
    for(const auto& particle:particles_){ChargedParticle proxy;proxy.velocity_m_per_s=particle.velocity_m_per_s;proxy.mass_kg=particle.mass_kg;proxy.charge_c=particle.charge_c;proxy.weight=particle.weight;result.particle_kinetic_energy_j+=particle.weight*kinetic_energy_j(proxy);}result.total_energy_j=result.field_energy_j+result.particle_kinetic_energy_j;result.charge_continuity_linf_residual=last_continuity_residual_;result.particle_count=particles_.size();result.absorbed_particles=absorbed_particles_;result.secondary_macro_weight=secondary_macro_weight_;result.collision_statistics=accumulated_collisions_;return result;
}


namespace {
std::size_t flat_index_3d(std::size_t ix,std::size_t iy,std::size_t iz,std::size_t nx,std::size_t ny){return (iz*ny+iy)*nx+ix;}

std::complex<double> fourier_coefficient_3d(std::span<const double> values,std::size_t nx,std::size_t ny,std::size_t nz,
                                            std::size_t mode_x,std::size_t mode_y,std::size_t mode_z){
    using Complex=std::complex<double>;Complex coefficient{};const double two_pi=2.0*std::numbers::pi;
    for(std::size_t iz=0;iz<nz;++iz){
        for(std::size_t iy=0;iy<ny;++iy){
            for(std::size_t ix=0;ix<nx;++ix){
                const double angle=-two_pi*(static_cast<double>(mode_x*ix)/static_cast<double>(nx)
                                           +static_cast<double>(mode_y*iy)/static_cast<double>(ny)
                                           +static_cast<double>(mode_z*iz)/static_cast<double>(nz));
                coefficient+=values[flat_index_3d(ix,iy,iz,nx,ny)]*Complex{std::cos(angle),std::sin(angle)};
            }
        }
    }
    return coefficient/static_cast<double>(nx*ny*nz);
}

std::array<long,3> signed_modes_3d(std::size_t mode_x,std::size_t mode_y,std::size_t mode_z,
                                   std::size_t nx,std::size_t ny,std::size_t nz){
    const long sx=(mode_x<=nx/2U)?static_cast<long>(mode_x):static_cast<long>(mode_x)-static_cast<long>(nx);
    const long sy=(mode_y<=ny/2U)?static_cast<long>(mode_y):static_cast<long>(mode_y)-static_cast<long>(ny);
    const long sz=(mode_z<=nz/2U)?static_cast<long>(mode_z):static_cast<long>(mode_z)-static_cast<long>(nz);
    return {sx,sy,sz};
}
}

std::size_t particle_cell_index_3d(Vec3 position_m,std::size_t nx,std::size_t ny,std::size_t nz,
                                   double length_x_m,double length_y_m,double length_z_m){
    if(nx<1U||ny<1U||nz<1U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(length_z_m>0.0))throw std::invalid_argument("invalid 3-D particle cell grid");
    const double dx_value=length_x_m/static_cast<double>(nx);
    const double dy_value=length_y_m/static_cast<double>(ny);
    const double dz_value=length_z_m/static_cast<double>(nz);
    const double x=wrap_periodic_value(position_m.x,length_x_m);
    const double y=wrap_periodic_value(position_m.y,length_y_m);
    const double z=wrap_periodic_value(position_m.z,length_z_m);
    const auto ix=std::min<std::size_t>(nx-1U,static_cast<std::size_t>(std::floor(x/dx_value)));
    const auto iy=std::min<std::size_t>(ny-1U,static_cast<std::size_t>(std::floor(y/dy_value)));
    const auto iz=std::min<std::size_t>(nz-1U,static_cast<std::size_t>(std::floor(z/dz_value)));
    return flat_index_3d(ix,iy,iz,nx,ny);
}

ParticleCellSort3DResult sort_particles_by_cell_3d(std::span<const PicParticle3D> particles,
                                                    std::size_t nx,std::size_t ny,std::size_t nz,
                                                    double length_x_m,double length_y_m,double length_z_m){
    if(nx<1U||ny<1U||nz<1U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(length_z_m>0.0))throw std::invalid_argument("invalid 3-D particle sorting grid");
    const std::size_t cells=nx*ny*nz;
    ParticleCellSort3DResult result;result.cell_counts.assign(cells,0U);result.cell_offsets.assign(cells+1U,0U);
    std::vector<std::size_t> particle_cells(particles.size());
    for(std::size_t index=0;index<particles.size();++index){
        if(!(particles[index].mass_kg>0.0)||!std::isfinite(particles[index].charge_c)||!std::isfinite(particles[index].weight))throw std::invalid_argument("invalid particle for 3-D sorting");
        const auto cell=particle_cell_index_3d(particles[index].position_m,nx,ny,nz,length_x_m,length_y_m,length_z_m);
        particle_cells[index]=cell;++result.cell_counts[cell];
    }
    for(std::size_t cell=0;cell<cells;++cell){result.cell_offsets[cell+1U]=result.cell_offsets[cell]+result.cell_counts[cell];if(result.cell_counts[cell]>0U)++result.occupied_cells;}
    result.sorted_particles.resize(particles.size());result.original_indices.resize(particles.size());
    auto write=result.cell_offsets;
    for(std::size_t index=0;index<particles.size();++index){
        const std::size_t destination=write[particle_cells[index]]++;
        result.sorted_particles[destination]=particles[index];
        result.original_indices[destination]=index;
    }
    return result;
}

ParticleGuardHalo3DResult classify_particle_guard_halos_3d(std::span<const PicParticle3D> particles,
                                                            const ParticleGuardHalo3DConfig& config){
    if(config.nx<1U||config.ny<1U||config.nz<1U||config.guard_cells<1U||!(config.length_x_m>0.0)||!(config.length_y_m>0.0)||!(config.length_z_m>0.0))throw std::invalid_argument("invalid 3-D particle guard halo config");
    ParticleGuardHalo3DResult result;
    const double tx=std::min(config.length_x_m,static_cast<double>(config.guard_cells)*config.length_x_m/static_cast<double>(config.nx));
    const double ty=std::min(config.length_y_m,static_cast<double>(config.guard_cells)*config.length_y_m/static_cast<double>(config.ny));
    const double tz=std::min(config.length_z_m,static_cast<double>(config.guard_cells)*config.length_z_m/static_cast<double>(config.nz));
    for(std::size_t index=0;index<particles.size();++index){
        Vec3 p=particles[index].position_m;
        const bool outside=(p.x<0.0||p.x>=config.length_x_m||p.y<0.0||p.y>=config.length_y_m||p.z<0.0||p.z>=config.length_z_m);
        if(outside&&!config.periodic){result.outside_domain.push_back(index);continue;}
        if(config.periodic){p.x=wrap_periodic_value(p.x,config.length_x_m);p.y=wrap_periodic_value(p.y,config.length_y_m);p.z=wrap_periodic_value(p.z,config.length_z_m);}
        bool halo=false;
        if(p.x<tx){result.x_minus.push_back(index);halo=true;}
        if(p.x>=config.length_x_m-tx){result.x_plus.push_back(index);halo=true;}
        if(p.y<ty){result.y_minus.push_back(index);halo=true;}
        if(p.y>=config.length_y_m-ty){result.y_plus.push_back(index);halo=true;}
        if(p.z<tz){result.z_minus.push_back(index);halo=true;}
        if(p.z>=config.length_z_m-tz){result.z_plus.push_back(index);halo=true;}
        if(!halo)result.interior.push_back(index);
    }
    return result;
}

namespace {
void validate_domain_grid(const PicDomainGrid3DConfig& config){
    if(config.global_nx<1U||config.global_ny<1U||config.global_nz<1U||
       config.domains_x<1U||config.domains_y<1U||config.domains_z<1U||
       config.domains_x>config.global_nx||config.domains_y>config.global_ny||config.domains_z>config.global_nz||
       !(config.length_x_m>0.0)||!(config.length_y_m>0.0)||!(config.length_z_m>0.0))
        throw std::invalid_argument("invalid 3-D PIC domain grid");
}
std::size_t block_start(std::size_t block,std::size_t blocks,std::size_t cells){return (block*cells)/blocks;}
std::size_t block_count(std::size_t block,std::size_t blocks,std::size_t cells){return block_start(block+1U,blocks,cells)-block_start(block,blocks,cells);}
std::size_t wrap_cell_index(long cell,std::size_t cells){
    const long span=static_cast<long>(cells);long value=cell%span;if(value<0)value+=span;return static_cast<std::size_t>(value);
}
std::size_t locate_domain_from_cell(std::size_t gx,std::size_t gy,std::size_t gz,std::span<const PicDomain3D> domains){
    for(const auto& domain:domains){
        if(gx>=domain.first_x&&gx<domain.first_x+domain.cells_x&&
           gy>=domain.first_y&&gy<domain.first_y+domain.cells_y&&
           gz>=domain.first_z&&gz<domain.first_z+domain.cells_z)return domain.domain_id;
    }
    return invalid_pic_domain_id_3d;
}
}

std::vector<PicDomain3D> make_pic_domain_grid_3d(const PicDomainGrid3DConfig& config){
    validate_domain_grid(config);
    std::vector<PicDomain3D> domains;domains.reserve(config.domains_x*config.domains_y*config.domains_z);
    const double dx_value=config.length_x_m/static_cast<double>(config.global_nx);
    const double dy_value=config.length_y_m/static_cast<double>(config.global_ny);
    const double dz_value=config.length_z_m/static_cast<double>(config.global_nz);
    for(std::size_t dz_i=0;dz_i<config.domains_z;++dz_i){
        for(std::size_t dy_i=0;dy_i<config.domains_y;++dy_i){
            for(std::size_t dx_i=0;dx_i<config.domains_x;++dx_i){
                PicDomain3D domain;domain.ix=dx_i;domain.iy=dy_i;domain.iz=dz_i;
                domain.domain_id=flat_index_3d(dx_i,dy_i,dz_i,config.domains_x,config.domains_y);
                domain.first_x=block_start(dx_i,config.domains_x,config.global_nx);
                domain.first_y=block_start(dy_i,config.domains_y,config.global_ny);
                domain.first_z=block_start(dz_i,config.domains_z,config.global_nz);
                domain.cells_x=block_count(dx_i,config.domains_x,config.global_nx);
                domain.cells_y=block_count(dy_i,config.domains_y,config.global_ny);
                domain.cells_z=block_count(dz_i,config.domains_z,config.global_nz);
                domain.minimum_m={static_cast<double>(domain.first_x)*dx_value,static_cast<double>(domain.first_y)*dy_value,static_cast<double>(domain.first_z)*dz_value};
                domain.maximum_m={static_cast<double>(domain.first_x+domain.cells_x)*dx_value,static_cast<double>(domain.first_y+domain.cells_y)*dy_value,static_cast<double>(domain.first_z+domain.cells_z)*dz_value};
                domains.push_back(domain);
            }
        }
    }
    return domains;
}

std::size_t locate_pic_domain_3d(Vec3 position_m,const PicDomainGrid3DConfig& config){
    validate_domain_grid(config);
    if(config.periodic){position_m.x=wrap_periodic_value(position_m.x,config.length_x_m);position_m.y=wrap_periodic_value(position_m.y,config.length_y_m);position_m.z=wrap_periodic_value(position_m.z,config.length_z_m);}else{
        if(position_m.x<0.0||position_m.x>=config.length_x_m||position_m.y<0.0||position_m.y>=config.length_y_m||position_m.z<0.0||position_m.z>=config.length_z_m)return invalid_pic_domain_id_3d;
    }
    const auto gx=std::min<std::size_t>(config.global_nx-1U,static_cast<std::size_t>(std::floor(position_m.x/(config.length_x_m/static_cast<double>(config.global_nx)))));
    const auto gy=std::min<std::size_t>(config.global_ny-1U,static_cast<std::size_t>(std::floor(position_m.y/(config.length_y_m/static_cast<double>(config.global_ny)))));
    const auto gz=std::min<std::size_t>(config.global_nz-1U,static_cast<std::size_t>(std::floor(position_m.z/(config.length_z_m/static_cast<double>(config.global_nz)))));
    const auto domains=make_pic_domain_grid_3d(config);
    return locate_domain_from_cell(gx,gy,gz,domains);
}

ParticleDomainMigration3D plan_particle_domain_migration_3d(std::span<const PicParticle3D> particles,const PicDomainGrid3DConfig& config){
    const auto domains=make_pic_domain_grid_3d(config);
    ParticleDomainMigration3D result;result.particles_by_domain.resize(domains.size());result.source_indices_by_domain.resize(domains.size());result.destination_domain_by_particle.assign(particles.size(),invalid_pic_domain_id_3d);
    for(std::size_t index=0;index<particles.size();++index){
        const auto domain_id=locate_pic_domain_3d(particles[index].position_m,config);
        result.destination_domain_by_particle[index]=domain_id;
        if(domain_id==invalid_pic_domain_id_3d){result.outside_indices.push_back(index);continue;}
        PicParticle3D particle=particles[index];
        if(config.periodic){particle.position_m.x=wrap_periodic_value(particle.position_m.x,config.length_x_m);particle.position_m.y=wrap_periodic_value(particle.position_m.y,config.length_y_m);particle.position_m.z=wrap_periodic_value(particle.position_m.z,config.length_z_m);}
        result.particles_by_domain[domain_id].push_back(particle);result.source_indices_by_domain[domain_id].push_back(index);
    }
    return result;
}

std::vector<ScalarGuardedBlock3D> exchange_scalar_guard_cells_3d(const std::vector<std::vector<double>>& cell_values_by_domain,
                                                                 std::span<const PicDomain3D> domains,
                                                                 const PicDomainGrid3DConfig& config,
                                                                 std::size_t guard_cells,double exterior_value){
    validate_domain_grid(config);
    if(guard_cells<1U||guard_cells>config.global_nx||guard_cells>config.global_ny||guard_cells>config.global_nz)throw std::invalid_argument("invalid scalar guard-cell width");
    if(cell_values_by_domain.size()!=domains.size())throw std::invalid_argument("domain value count must match domains");
    std::vector<ScalarGuardedBlock3D> result(domains.size());
    for(const auto& domain:domains){
        if(domain.domain_id>=domains.size())throw std::invalid_argument("domain id out of range");
        const auto expected=domain.cells_x*domain.cells_y*domain.cells_z;
        if(cell_values_by_domain[domain.domain_id].size()!=expected)throw std::invalid_argument("scalar domain block size mismatch");
    }
    const long nx_global=static_cast<long>(config.global_nx),ny_global=static_cast<long>(config.global_ny),nz_global=static_cast<long>(config.global_nz);
    for(const auto& domain:domains){
        ScalarGuardedBlock3D block;block.domain_id=domain.domain_id;block.interior_nx=domain.cells_x;block.interior_ny=domain.cells_y;block.interior_nz=domain.cells_z;block.guard_cells=guard_cells;
        const std::size_t px=block.padded_nx(),py=block.padded_ny(),pz=block.padded_nz();block.values.assign(px*py*pz,exterior_value);
        for(std::size_t iz=0;iz<pz;++iz){for(std::size_t iy=0;iy<py;++iy){for(std::size_t ix=0;ix<px;++ix){
            long gx=static_cast<long>(domain.first_x)+static_cast<long>(ix)-static_cast<long>(guard_cells);
            long gy=static_cast<long>(domain.first_y)+static_cast<long>(iy)-static_cast<long>(guard_cells);
            long gz=static_cast<long>(domain.first_z)+static_cast<long>(iz)-static_cast<long>(guard_cells);
            if(config.periodic){gx=static_cast<long>(wrap_cell_index(gx,config.global_nx));gy=static_cast<long>(wrap_cell_index(gy,config.global_ny));gz=static_cast<long>(wrap_cell_index(gz,config.global_nz));}
            if(gx<0||gx>=nx_global||gy<0||gy>=ny_global||gz<0||gz>=nz_global){continue;}
            const auto owner=locate_domain_from_cell(static_cast<std::size_t>(gx),static_cast<std::size_t>(gy),static_cast<std::size_t>(gz),domains);
            if(owner==invalid_pic_domain_id_3d)continue;
            const auto& owner_domain=domains[owner];
            const std::size_t lx=static_cast<std::size_t>(gx)-owner_domain.first_x;
            const std::size_t ly=static_cast<std::size_t>(gy)-owner_domain.first_y;
            const std::size_t lz=static_cast<std::size_t>(gz)-owner_domain.first_z;
            const double value=cell_values_by_domain[owner][flat_index_3d(lx,ly,lz,owner_domain.cells_x,owner_domain.cells_y)];
            block.values[flat_index_3d(ix,iy,iz,px,py)]=value;
        }}}
        result[domain.domain_id]=std::move(block);
    }
    return result;
}



namespace {
ParticleMigrationMessage3D& find_or_add_particle_message(std::vector<ParticleMigrationMessage3D>& messages,
                                                          std::size_t source_domain_id,
                                                          std::size_t destination_domain_id){
    for(auto& message:messages){
        if(message.source_domain_id==source_domain_id&&message.destination_domain_id==destination_domain_id)return message;
    }
    ParticleMigrationMessage3D message;message.source_domain_id=source_domain_id;message.destination_domain_id=destination_domain_id;
    messages.push_back(std::move(message));return messages.back();
}

ScalarGuardCellMessage3D& find_or_add_scalar_guard_message(std::vector<ScalarGuardCellMessage3D>& messages,
                                                            std::size_t source_domain_id,
                                                            std::size_t destination_domain_id,
                                                            int offset_x,int offset_y,int offset_z){
    for(auto& message:messages){
        if(message.source_domain_id==source_domain_id&&message.destination_domain_id==destination_domain_id&&
           message.offset_x==offset_x&&message.offset_y==offset_y&&message.offset_z==offset_z)return message;
    }
    ScalarGuardCellMessage3D message;message.source_domain_id=source_domain_id;message.destination_domain_id=destination_domain_id;
    message.offset_x=offset_x;message.offset_y=offset_y;message.offset_z=offset_z;messages.push_back(std::move(message));return messages.back();
}

void validate_domain_span(std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config){
    validate_domain_grid(config);
    if(domains.size()!=config.domains_x*config.domains_y*config.domains_z)throw std::invalid_argument("PIC domain span size does not match grid");
    for(const auto& domain:domains){if(domain.domain_id>=domains.size())throw std::invalid_argument("PIC domain id out of range");}
}

int offset_sign_for_padded_index(std::size_t padded_index,std::size_t guard_cells,std::size_t interior_cells){
    if(padded_index<guard_cells)return -1;
    if(padded_index>=guard_cells+interior_cells)return 1;
    return 0;
}
}

PackedParticleMigration3D pack_particle_migration_messages_3d(const std::vector<std::vector<PicParticle3D>>& particles_by_source_domain,
                                                               std::span<const PicDomain3D> domains,
                                                               const PicDomainGrid3DConfig& config){
    validate_domain_span(domains,config);
    if(particles_by_source_domain.size()!=domains.size())throw std::invalid_argument("particle source bucket count must match PIC domains");
    PackedParticleMigration3D result;result.retained_by_domain.resize(domains.size());result.retained_source_indices_by_domain.resize(domains.size());
    for(const auto& domain:domains){
        const auto source_domain_id=domain.domain_id;
        const auto& source_particles=particles_by_source_domain[source_domain_id];
        for(std::size_t source_index=0;source_index<source_particles.size();++source_index){
            PicParticle3D particle=source_particles[source_index];
            const auto destination_domain_id=locate_pic_domain_3d(particle.position_m,config);
            if(destination_domain_id==invalid_pic_domain_id_3d){result.outside_particles.push_back({source_domain_id,source_index});continue;}
            if(config.periodic){particle.position_m.x=wrap_periodic_value(particle.position_m.x,config.length_x_m);particle.position_m.y=wrap_periodic_value(particle.position_m.y,config.length_y_m);particle.position_m.z=wrap_periodic_value(particle.position_m.z,config.length_z_m);}
            if(destination_domain_id==source_domain_id){
                result.retained_by_domain[source_domain_id].push_back(particle);
                result.retained_source_indices_by_domain[source_domain_id].push_back(source_index);
            }else{
                auto& message=find_or_add_particle_message(result.messages,source_domain_id,destination_domain_id);
                message.particles.push_back(particle);message.source_indices.push_back(source_index);
            }
        }
    }
    return result;
}

std::vector<std::vector<PicParticle3D>> apply_particle_migration_messages_3d(
    std::vector<std::vector<PicParticle3D>> retained_by_domain,std::span<const ParticleMigrationMessage3D> messages){
    for(const auto& message:messages){
        if(message.destination_domain_id>=retained_by_domain.size())throw std::invalid_argument("particle migration message destination out of range");
        if(message.particles.size()!=message.source_indices.size())throw std::invalid_argument("particle migration message source-index size mismatch");
        auto& destination=retained_by_domain[message.destination_domain_id];
        destination.insert(destination.end(),message.particles.begin(),message.particles.end());
    }
    return retained_by_domain;
}

std::vector<ScalarGuardCellMessage3D> pack_scalar_guard_cell_messages_3d(
    const std::vector<std::vector<double>>& cell_values_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    std::size_t guard_cells){
    validate_domain_span(domains,config);
    if(guard_cells<1U||guard_cells>config.global_nx||guard_cells>config.global_ny||guard_cells>config.global_nz)throw std::invalid_argument("invalid scalar guard-cell width");
    if(cell_values_by_domain.size()!=domains.size())throw std::invalid_argument("domain value count must match domains");
    for(const auto& domain:domains){
        const auto expected=domain.cells_x*domain.cells_y*domain.cells_z;
        if(cell_values_by_domain[domain.domain_id].size()!=expected)throw std::invalid_argument("scalar domain block size mismatch");
    }
    std::vector<ScalarGuardCellMessage3D> messages;
    const long nx_global=static_cast<long>(config.global_nx),ny_global=static_cast<long>(config.global_ny),nz_global=static_cast<long>(config.global_nz);
    for(const auto& destination_domain:domains){
        const std::size_t px=destination_domain.cells_x+2U*guard_cells;
        const std::size_t py=destination_domain.cells_y+2U*guard_cells;
        const std::size_t pz=destination_domain.cells_z+2U*guard_cells;
        for(std::size_t iz=0;iz<pz;++iz){for(std::size_t iy=0;iy<py;++iy){for(std::size_t ix=0;ix<px;++ix){
            const bool interior=(ix>=guard_cells&&ix<guard_cells+destination_domain.cells_x&&
                                 iy>=guard_cells&&iy<guard_cells+destination_domain.cells_y&&
                                 iz>=guard_cells&&iz<guard_cells+destination_domain.cells_z);
            if(interior)continue;
            long gx=static_cast<long>(destination_domain.first_x)+static_cast<long>(ix)-static_cast<long>(guard_cells);
            long gy=static_cast<long>(destination_domain.first_y)+static_cast<long>(iy)-static_cast<long>(guard_cells);
            long gz=static_cast<long>(destination_domain.first_z)+static_cast<long>(iz)-static_cast<long>(guard_cells);
            if(config.periodic){gx=static_cast<long>(wrap_cell_index(gx,config.global_nx));gy=static_cast<long>(wrap_cell_index(gy,config.global_ny));gz=static_cast<long>(wrap_cell_index(gz,config.global_nz));}
            if(gx<0||gx>=nx_global||gy<0||gy>=ny_global||gz<0||gz>=nz_global)continue;
            const auto source_domain_id=locate_domain_from_cell(static_cast<std::size_t>(gx),static_cast<std::size_t>(gy),static_cast<std::size_t>(gz),domains);
            if(source_domain_id==invalid_pic_domain_id_3d)continue;
            const auto& source_domain=domains[source_domain_id];
            const std::size_t sx=static_cast<std::size_t>(gx)-source_domain.first_x;
            const std::size_t sy=static_cast<std::size_t>(gy)-source_domain.first_y;
            const std::size_t sz=static_cast<std::size_t>(gz)-source_domain.first_z;
            const double value=cell_values_by_domain[source_domain_id][flat_index_3d(sx,sy,sz,source_domain.cells_x,source_domain.cells_y)];
            const int ox=offset_sign_for_padded_index(ix,guard_cells,destination_domain.cells_x);
            const int oy=offset_sign_for_padded_index(iy,guard_cells,destination_domain.cells_y);
            const int oz=offset_sign_for_padded_index(iz,guard_cells,destination_domain.cells_z);
            auto& message=find_or_add_scalar_guard_message(messages,source_domain_id,destination_domain.domain_id,ox,oy,oz);
            message.values.push_back({ix,iy,iz,value});
        }}}
    }
    return messages;
}


ElectromagneticGuardedFieldBlocks3D exchange_electromagnetic_field_guard_cells_3d(
    const ElectromagneticFieldBlocks3D& fields_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    std::size_t guard_cells,double exterior_value){
    const auto domain_count=domains.size();
    auto check_component=[&](const std::vector<std::vector<double>>& component,const char* name){
        if(component.size()!=domain_count)throw std::invalid_argument(std::string(name)+" field domain count mismatch");
        for(const auto& domain:domains){
            const auto expected=domain.cells_x*domain.cells_y*domain.cells_z;
            if(component[domain.domain_id].size()!=expected)throw std::invalid_argument(std::string(name)+" field block size mismatch");
        }
    };
    check_component(fields_by_domain.electric_x_by_domain,"Ex");
    check_component(fields_by_domain.electric_y_by_domain,"Ey");
    check_component(fields_by_domain.electric_z_by_domain,"Ez");
    check_component(fields_by_domain.magnetic_x_by_domain,"Bx");
    check_component(fields_by_domain.magnetic_y_by_domain,"By");
    check_component(fields_by_domain.magnetic_z_by_domain,"Bz");
    ElectromagneticGuardedFieldBlocks3D result;
    result.electric_x=exchange_scalar_guard_cells_3d(fields_by_domain.electric_x_by_domain,domains,config,guard_cells,exterior_value);
    result.electric_y=exchange_scalar_guard_cells_3d(fields_by_domain.electric_y_by_domain,domains,config,guard_cells,exterior_value);
    result.electric_z=exchange_scalar_guard_cells_3d(fields_by_domain.electric_z_by_domain,domains,config,guard_cells,exterior_value);
    result.magnetic_x=exchange_scalar_guard_cells_3d(fields_by_domain.magnetic_x_by_domain,domains,config,guard_cells,exterior_value);
    result.magnetic_y=exchange_scalar_guard_cells_3d(fields_by_domain.magnetic_y_by_domain,domains,config,guard_cells,exterior_value);
    result.magnetic_z=exchange_scalar_guard_cells_3d(fields_by_domain.magnetic_z_by_domain,domains,config,guard_cells,exterior_value);
    return result;
}

std::vector<ScalarGuardedBlock3D> apply_scalar_guard_cell_messages_3d(
    const std::vector<std::vector<double>>& cell_values_by_domain,
    std::span<const PicDomain3D> domains,std::size_t guard_cells,double exterior_value,
    std::span<const ScalarGuardCellMessage3D> messages){
    if(guard_cells<1U)throw std::invalid_argument("invalid scalar guard-cell width");
    if(cell_values_by_domain.size()!=domains.size())throw std::invalid_argument("domain value count must match domains");
    std::vector<ScalarGuardedBlock3D> blocks(domains.size());
    for(const auto& domain:domains){
        if(domain.domain_id>=domains.size())throw std::invalid_argument("domain id out of range");
        const auto expected=domain.cells_x*domain.cells_y*domain.cells_z;
        if(cell_values_by_domain[domain.domain_id].size()!=expected)throw std::invalid_argument("scalar domain block size mismatch");
        ScalarGuardedBlock3D block;block.domain_id=domain.domain_id;block.interior_nx=domain.cells_x;block.interior_ny=domain.cells_y;block.interior_nz=domain.cells_z;block.guard_cells=guard_cells;
        const std::size_t px=block.padded_nx(),py=block.padded_ny(),pz=block.padded_nz();block.values.assign(px*py*pz,exterior_value);
        for(std::size_t iz=0;iz<domain.cells_z;++iz){for(std::size_t iy=0;iy<domain.cells_y;++iy){for(std::size_t ix=0;ix<domain.cells_x;++ix){
            block.values[flat_index_3d(ix+guard_cells,iy+guard_cells,iz+guard_cells,px,py)]=cell_values_by_domain[domain.domain_id][flat_index_3d(ix,iy,iz,domain.cells_x,domain.cells_y)];
        }}}
        blocks[domain.domain_id]=std::move(block);
    }
    for(const auto& message:messages){
        if(message.destination_domain_id>=blocks.size())throw std::invalid_argument("scalar guard message destination out of range");
        auto& block=blocks[message.destination_domain_id];
        const std::size_t px=block.padded_nx(),py=block.padded_ny(),pz=block.padded_nz();
        for(const auto& entry:message.values){
            if(entry.padded_x>=px||entry.padded_y>=py||entry.padded_z>=pz)throw std::invalid_argument("scalar guard message entry out of range");
            block.values[flat_index_3d(entry.padded_x,entry.padded_y,entry.padded_z,px,py)]=entry.value;
        }
    }
    return blocks;
}


PicRankTopology3D make_pic_rank_topology_3d(std::span<const PicDomain3D> domains,
                                            std::size_t ranks_x,std::size_t ranks_y,std::size_t ranks_z){
    if(domains.empty()||ranks_x==0U||ranks_y==0U||ranks_z==0U)throw std::invalid_argument("invalid PIC rank topology dimensions");
    std::size_t domain_count_x=0U,domain_count_y=0U,domain_count_z=0U;
    for(const auto& domain:domains){
        if(domain.domain_id>=domains.size())throw std::invalid_argument("PIC domain id out of range for rank topology");
        domain_count_x=std::max(domain_count_x,domain.ix+1U);
        domain_count_y=std::max(domain_count_y,domain.iy+1U);
        domain_count_z=std::max(domain_count_z,domain.iz+1U);
    }
    if(ranks_x>domain_count_x||ranks_y>domain_count_y||ranks_z>domain_count_z)throw std::invalid_argument("PIC rank topology cannot exceed domain grid dimensions");
    PicRankTopology3D topology;topology.rank_count=ranks_x*ranks_y*ranks_z;topology.domain_to_rank.assign(domains.size(),0U);
    for(const auto& domain:domains){
        const std::size_t rx=std::min(ranks_x-1U,domain.ix*ranks_x/domain_count_x);
        const std::size_t ry=std::min(ranks_y-1U,domain.iy*ranks_y/domain_count_y);
        const std::size_t rz=std::min(ranks_z-1U,domain.iz*ranks_z/domain_count_z);
        topology.domain_to_rank[domain.domain_id]=(rz*ranks_y+ry)*ranks_x+rx;
    }
    return topology;
}

std::size_t rank_for_pic_domain_3d(const PicRankTopology3D& topology,std::size_t domain_id){
    if(domain_id>=topology.domain_to_rank.size())throw std::invalid_argument("PIC domain id out of range for rank lookup");
    const std::size_t rank=topology.domain_to_rank[domain_id];
    if(rank>=topology.rank_count)throw std::invalid_argument("PIC topology contains invalid rank id");
    return rank;
}

InMemoryPicTransport3D::InMemoryPicTransport3D(std::size_t rank_count):rank_count_(rank_count),inbox_by_rank_(rank_count){
    if(rank_count==0U)throw std::invalid_argument("PIC in-memory transport requires at least one rank");
}

void InMemoryPicTransport3D::post(std::span<const PicTransportEnvelope3D> envelopes){
    for(const auto& envelope:envelopes){
        if(envelope.destination_rank>=rank_count_)throw std::invalid_argument("PIC transport destination rank out of range");
        if(envelope.source_rank>=rank_count_)throw std::invalid_argument("PIC transport source rank out of range");
        inbox_by_rank_[envelope.destination_rank].push_back(envelope);
    }
}

std::vector<PicTransportEnvelope3D> InMemoryPicTransport3D::receive_rank(std::size_t rank){
    if(rank>=rank_count_)throw std::invalid_argument("PIC transport receive rank out of range");
    std::vector<PicTransportEnvelope3D> result;result.swap(inbox_by_rank_[rank]);return result;
}

std::size_t InMemoryPicTransport3D::pending_count() const noexcept{
    std::size_t total=0U;for(const auto& inbox:inbox_by_rank_)total+=inbox.size();return total;
}

void InMemoryPicTransport3D::clear() noexcept{for(auto& inbox:inbox_by_rank_)inbox.clear();}

std::vector<PicTransportEnvelope3D> make_pic_transport_envelopes_3d(
    std::span<const ParticleMigrationMessage3D> particle_messages,
    std::span<const ScalarGuardCellMessage3D> scalar_guard_messages,const PicRankTopology3D& topology){
    if(topology.rank_count==0U)throw std::invalid_argument("PIC transport topology has no ranks");
    std::vector<PicTransportEnvelope3D> envelopes;envelopes.reserve(particle_messages.size()+scalar_guard_messages.size());
    for(const auto& message:particle_messages){
        PicTransportEnvelope3D envelope;envelope.kind=PicTransportPayloadKind3D::particle_migration;
        envelope.source_domain_id=message.source_domain_id;envelope.destination_domain_id=message.destination_domain_id;
        envelope.source_rank=rank_for_pic_domain_3d(topology,message.source_domain_id);
        envelope.destination_rank=rank_for_pic_domain_3d(topology,message.destination_domain_id);
        envelope.particle_migration=message;envelopes.push_back(std::move(envelope));
    }
    for(const auto& message:scalar_guard_messages){
        PicTransportEnvelope3D envelope;envelope.kind=PicTransportPayloadKind3D::scalar_guard_cells;
        envelope.source_domain_id=message.source_domain_id;envelope.destination_domain_id=message.destination_domain_id;
        envelope.source_rank=rank_for_pic_domain_3d(topology,message.source_domain_id);
        envelope.destination_rank=rank_for_pic_domain_3d(topology,message.destination_domain_id);
        envelope.scalar_guard=message;envelopes.push_back(std::move(envelope));
    }
    return envelopes;
}

PicTransportDiagnostics3D summarize_pic_transport_3d(std::span<const PicTransportEnvelope3D> envelopes){
    PicTransportDiagnostics3D diagnostics;diagnostics.envelopes=envelopes.size();
    for(const auto& envelope:envelopes){
        if(envelope.source_rank==envelope.destination_rank)++diagnostics.same_rank_messages;else ++diagnostics.remote_rank_messages;
        if(envelope.kind==PicTransportPayloadKind3D::particle_migration){
            ++diagnostics.particle_messages;diagnostics.particle_payload_count+=envelope.particle_migration.particles.size();
        }else{
            ++diagnostics.scalar_guard_messages;diagnostics.scalar_guard_value_count+=envelope.scalar_guard.values.size();
        }
    }
    return diagnostics;
}

PicTransportedExchange3D exchange_pic_messages_in_memory_3d(
    std::vector<std::vector<PicParticle3D>> retained_by_domain,
    std::span<const ParticleMigrationMessage3D> particle_messages,
    const std::vector<std::vector<double>>& cell_values_by_domain,
    std::span<const ScalarGuardCellMessage3D> scalar_guard_messages,
    std::span<const PicDomain3D> domains,std::size_t guard_cells,double exterior_value,
    const PicRankTopology3D& topology){
    const auto envelopes=make_pic_transport_envelopes_3d(particle_messages,scalar_guard_messages,topology);
    InMemoryPicTransport3D transport(topology.rank_count);transport.post(envelopes);
    std::vector<ParticleMigrationMessage3D> delivered_particle_messages;
    std::vector<ScalarGuardCellMessage3D> delivered_guard_messages;
    for(std::size_t rank=0;rank<topology.rank_count;++rank){
        const auto inbox=transport.receive_rank(rank);
        for(const auto& envelope:inbox){
            if(envelope.kind==PicTransportPayloadKind3D::particle_migration)delivered_particle_messages.push_back(envelope.particle_migration);
            else delivered_guard_messages.push_back(envelope.scalar_guard);
        }
    }
    PicTransportedExchange3D result;
    result.particles_by_domain=apply_particle_migration_messages_3d(std::move(retained_by_domain),delivered_particle_messages);
    result.guarded_blocks=apply_scalar_guard_cell_messages_3d(cell_values_by_domain,domains,guard_cells,exterior_value,delivered_guard_messages);
    result.diagnostics=summarize_pic_transport_3d(envelopes);
    return result;
}

namespace {
constexpr std::uint32_t pic_transport_magic_3d = 0x50494333U; // PIC3
constexpr std::uint32_t pic_transport_version_3d = 1U;

class ByteWriter3D {
public:
    template<class T> void write(T value) {
        static_assert(std::is_trivially_copyable_v<T>, "PIC transport writes only trivially copyable values");
        const auto* raw = reinterpret_cast<const std::uint8_t*>(&value);
        bytes_.insert(bytes_.end(), raw, raw + sizeof(T));
    }
    void write_vec3(Vec3 value){write(value.x);write(value.y);write(value.z);} 
    void write_particle(const PicParticle3D& particle){write_vec3(particle.position_m);write_vec3(particle.velocity_m_per_s);write(particle.charge_c);write(particle.mass_kg);write(particle.weight);} 
    void write_bytes(std::span<const std::uint8_t> payload){bytes_.insert(bytes_.end(),payload.begin(),payload.end());}
    [[nodiscard]] PicSerializedEnvelope3D finish() && { return {std::move(bytes_)}; }
private:
    std::vector<std::uint8_t> bytes_;
};

class ByteReader3D {
public:
    explicit ByteReader3D(std::span<const std::uint8_t> bytes):bytes_(bytes){}
    template<class T> [[nodiscard]] T read(){
        static_assert(std::is_trivially_copyable_v<T>, "PIC transport reads only trivially copyable values");
        if(offset_ + sizeof(T) > bytes_.size())throw std::invalid_argument("truncated PIC transport envelope");
        T value{};std::memcpy(&value,bytes_.data()+offset_,sizeof(T));offset_+=sizeof(T);return value;
    }
    [[nodiscard]] Vec3 read_vec3(){return {read<double>(),read<double>(),read<double>()};}
    [[nodiscard]] PicParticle3D read_particle(){PicParticle3D p;p.position_m=read_vec3();p.velocity_m_per_s=read_vec3();p.charge_c=read<double>();p.mass_kg=read<double>();p.weight=read<double>();return p;}
    void require_finished() const { if(offset_!=bytes_.size())throw std::invalid_argument("PIC transport envelope has trailing bytes"); }
private:
    std::span<const std::uint8_t> bytes_;
    std::size_t offset_{};
};

void validate_serialized_rank(std::size_t rank,std::size_t rank_count,const char* label){
    if(rank>=rank_count)throw std::invalid_argument(std::string("PIC serialized exchange ")+label+" rank out of range");
}
}

PicSerializedEnvelope3D serialize_pic_transport_envelope_3d(const PicTransportEnvelope3D& envelope){
    ByteWriter3D writer;writer.write(pic_transport_magic_3d);writer.write(pic_transport_version_3d);
    writer.write(static_cast<std::uint32_t>(envelope.kind));writer.write(envelope.source_domain_id);writer.write(envelope.destination_domain_id);writer.write(envelope.source_rank);writer.write(envelope.destination_rank);
    if(envelope.kind==PicTransportPayloadKind3D::particle_migration){
        const auto& message=envelope.particle_migration;
        if(message.particles.size()!=message.source_indices.size())throw std::invalid_argument("PIC particle message size mismatch before serialization");
        if(message.source_domain_id!=envelope.source_domain_id||message.destination_domain_id!=envelope.destination_domain_id)throw std::invalid_argument("PIC particle envelope domain mismatch");
        writer.write(message.source_domain_id);writer.write(message.destination_domain_id);writer.write(message.particles.size());
        for(std::size_t i=0;i<message.particles.size();++i){writer.write_particle(message.particles[i]);writer.write(message.source_indices[i]);}
    }else{
        const auto& message=envelope.scalar_guard;
        if(message.source_domain_id!=envelope.source_domain_id||message.destination_domain_id!=envelope.destination_domain_id)throw std::invalid_argument("PIC scalar guard envelope domain mismatch");
        writer.write(message.source_domain_id);writer.write(message.destination_domain_id);writer.write(message.offset_x);writer.write(message.offset_y);writer.write(message.offset_z);writer.write(message.values.size());
        for(const auto& value:message.values){writer.write(value.padded_x);writer.write(value.padded_y);writer.write(value.padded_z);writer.write(value.value);} 
    }
    return std::move(writer).finish();
}

PicTransportEnvelope3D deserialize_pic_transport_envelope_3d(const PicSerializedEnvelope3D& serialized){
    if(serialized.bytes.empty())throw std::invalid_argument("empty PIC transport envelope");
    ByteReader3D reader(serialized.bytes);const auto magic=reader.read<std::uint32_t>();const auto version=reader.read<std::uint32_t>();
    if(magic!=pic_transport_magic_3d||version!=pic_transport_version_3d)throw std::invalid_argument("unsupported PIC transport envelope version");
    PicTransportEnvelope3D envelope;const auto kind_value=reader.read<std::uint32_t>();
    if(kind_value>static_cast<std::uint32_t>(PicTransportPayloadKind3D::scalar_guard_cells))throw std::invalid_argument("invalid PIC transport payload kind");
    envelope.kind=static_cast<PicTransportPayloadKind3D>(kind_value);envelope.source_domain_id=reader.read<std::size_t>();envelope.destination_domain_id=reader.read<std::size_t>();envelope.source_rank=reader.read<std::size_t>();envelope.destination_rank=reader.read<std::size_t>();
    if(envelope.kind==PicTransportPayloadKind3D::particle_migration){
        ParticleMigrationMessage3D message;message.source_domain_id=reader.read<std::size_t>();message.destination_domain_id=reader.read<std::size_t>();const auto count=reader.read<std::size_t>();
        if(count>static_cast<std::size_t>(1U<<28U))throw std::invalid_argument("unreasonable PIC particle message payload");
        message.particles.reserve(count);message.source_indices.reserve(count);
        for(std::size_t i=0;i<count;++i){message.particles.push_back(reader.read_particle());message.source_indices.push_back(reader.read<std::size_t>());} 
        if(message.source_domain_id!=envelope.source_domain_id||message.destination_domain_id!=envelope.destination_domain_id)throw std::invalid_argument("PIC particle message/envelope domain mismatch after deserialize");
        envelope.particle_migration=std::move(message);
    }else{
        ScalarGuardCellMessage3D message;message.source_domain_id=reader.read<std::size_t>();message.destination_domain_id=reader.read<std::size_t>();message.offset_x=reader.read<int>();message.offset_y=reader.read<int>();message.offset_z=reader.read<int>();const auto count=reader.read<std::size_t>();
        if(count>static_cast<std::size_t>(1U<<28U))throw std::invalid_argument("unreasonable PIC scalar guard message payload");
        message.values.reserve(count);for(std::size_t i=0;i<count;++i){ScalarGuardCellValue3D v;v.padded_x=reader.read<std::size_t>();v.padded_y=reader.read<std::size_t>();v.padded_z=reader.read<std::size_t>();v.value=reader.read<double>();message.values.push_back(v);} 
        if(message.source_domain_id!=envelope.source_domain_id||message.destination_domain_id!=envelope.destination_domain_id)throw std::invalid_argument("PIC scalar guard message/envelope domain mismatch after deserialize");
        envelope.scalar_guard=std::move(message);
    }
    reader.require_finished();return envelope;
}

std::vector<PicSerializedEnvelope3D> serialize_pic_transport_envelopes_3d(std::span<const PicTransportEnvelope3D> envelopes){
    std::vector<PicSerializedEnvelope3D> result;result.reserve(envelopes.size());for(const auto& envelope:envelopes)result.push_back(serialize_pic_transport_envelope_3d(envelope));return result;
}

std::vector<PicTransportEnvelope3D> deserialize_pic_transport_envelopes_3d(std::span<const PicSerializedEnvelope3D> serialized){
    std::vector<PicTransportEnvelope3D> result;result.reserve(serialized.size());for(const auto& item:serialized)result.push_back(deserialize_pic_transport_envelope_3d(item));return result;
}

PicSerializedExchangePlan3D plan_serialized_pic_exchange_3d(std::span<const PicTransportEnvelope3D> envelopes,std::size_t rank_count){
    if(rank_count==0U)throw std::invalid_argument("PIC serialized exchange plan requires ranks");
    PicSerializedExchangePlan3D plan;plan.rank_count=rank_count;plan.envelopes_by_destination_rank.resize(rank_count);plan.message_count_by_destination_rank.assign(rank_count,0U);plan.byte_count_by_destination_rank.assign(rank_count,0U);
    for(const auto& envelope:envelopes){validate_serialized_rank(envelope.source_rank,rank_count,"source");validate_serialized_rank(envelope.destination_rank,rank_count,"destination");auto serialized=serialize_pic_transport_envelope_3d(envelope);plan.byte_count_by_destination_rank[envelope.destination_rank]+=serialized.bytes.size();++plan.message_count_by_destination_rank[envelope.destination_rank];plan.envelopes_by_destination_rank[envelope.destination_rank].push_back(std::move(serialized));}
    return plan;
}

PicTransportDiagnostics3D summarize_pic_serialized_transport_3d(std::span<const PicSerializedEnvelope3D> serialized){
    return summarize_pic_transport_3d(deserialize_pic_transport_envelopes_3d(serialized));
}


DistributedPicExchangeRound3D run_serialized_distributed_pic_exchange_round_3d(
    const std::vector<std::vector<PicParticle3D>>& particles_by_domain,
    const std::vector<std::vector<double>>& scalar_cell_values_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    std::size_t guard_cells,double exterior_value,const PicRankTopology3D& topology){
    validate_domain_span(domains,config);
    if(particles_by_domain.size()!=domains.size())throw std::invalid_argument("distributed PIC exchange particle domain count mismatch");
    if(scalar_cell_values_by_domain.size()!=domains.size())throw std::invalid_argument("distributed PIC exchange scalar domain count mismatch");
    if(topology.domain_to_rank.size()!=domains.size())throw std::invalid_argument("distributed PIC exchange topology/domain mismatch");
    const auto packed_particles=pack_particle_migration_messages_3d(particles_by_domain,domains,config);
    const auto guard_messages=pack_scalar_guard_cell_messages_3d(scalar_cell_values_by_domain,domains,config,guard_cells);
    const auto envelopes=make_pic_transport_envelopes_3d(packed_particles.messages,guard_messages,topology);
    auto plan=plan_serialized_pic_exchange_3d(envelopes,topology.rank_count);
    std::vector<PicTransportEnvelope3D> delivered_envelopes;
    for(const auto& bucket:plan.envelopes_by_destination_rank){
        const auto decoded=deserialize_pic_transport_envelopes_3d(bucket);
        delivered_envelopes.insert(delivered_envelopes.end(),decoded.begin(),decoded.end());
    }
    std::vector<ParticleMigrationMessage3D> delivered_particles;
    std::vector<ScalarGuardCellMessage3D> delivered_guards;
    for(const auto& envelope:delivered_envelopes){
        if(envelope.kind==PicTransportPayloadKind3D::particle_migration)delivered_particles.push_back(envelope.particle_migration);
        else delivered_guards.push_back(envelope.scalar_guard);
    }
    DistributedPicExchangeRound3D result;
    result.particles_by_domain=apply_particle_migration_messages_3d(packed_particles.retained_by_domain,delivered_particles);
    result.scalar_guarded_blocks=apply_scalar_guard_cell_messages_3d(scalar_cell_values_by_domain,domains,guard_cells,exterior_value,delivered_guards);
    result.transport_diagnostics=summarize_pic_transport_3d(envelopes);
    result.serialized_plan=std::move(plan);
    for(std::size_t rank=0;rank<result.serialized_plan.rank_count;++rank){
        result.serialized_message_count+=result.serialized_plan.message_count_by_destination_rank[rank];
        result.serialized_byte_count+=result.serialized_plan.byte_count_by_destination_rank[rank];
    }
    result.outside_particles=packed_particles.outside_particles;
    return result;
}


namespace {
std::size_t count_pic_particles_by_domain_3d(const std::vector<std::vector<PicParticle3D>>& buckets){
    std::size_t total=0U;for(const auto& bucket:buckets)total+=bucket.size();return total;
}
void validate_em_domain_fields_3d(const ElectromagneticFieldBlocks3D& fields,std::span<const PicDomain3D> domains){
    const std::size_t n=domains.size();
    if(fields.electric_x_by_domain.size()!=n||fields.electric_y_by_domain.size()!=n||fields.electric_z_by_domain.size()!=n||
       fields.magnetic_x_by_domain.size()!=n||fields.magnetic_y_by_domain.size()!=n||fields.magnetic_z_by_domain.size()!=n)
        throw std::invalid_argument("distributed PIC step field domain count mismatch");
    for(const auto& domain:domains){
        const std::size_t cells=domain.cells_x*domain.cells_y*domain.cells_z;
        if(fields.electric_x_by_domain[domain.domain_id].size()!=cells||fields.electric_y_by_domain[domain.domain_id].size()!=cells||
           fields.electric_z_by_domain[domain.domain_id].size()!=cells||fields.magnetic_x_by_domain[domain.domain_id].size()!=cells||
           fields.magnetic_y_by_domain[domain.domain_id].size()!=cells||fields.magnetic_z_by_domain[domain.domain_id].size()!=cells)
            throw std::invalid_argument("distributed PIC step field array size mismatch");
    }
}
double local_sample_field_3d(std::span<const double> values,const PicDomain3D& domain,Vec3 position){
    if(values.empty())return 0.0;
    const double dx=(domain.maximum_m.x-domain.minimum_m.x)/static_cast<double>(domain.cells_x);
    const double dy=(domain.maximum_m.y-domain.minimum_m.y)/static_cast<double>(domain.cells_y);
    const double dz=(domain.maximum_m.z-domain.minimum_m.z)/static_cast<double>(domain.cells_z);
    const auto clamp_cell=[](double u,std::size_t n)->std::pair<std::size_t,double>{
        if(u<=0.0)return {0U,0.0};
        const double upper=static_cast<double>(n)-1.0e-12;
        if(u>=upper){return {n-1U,0.0};}
        const double floor_u=std::floor(u);
        return {static_cast<std::size_t>(floor_u),u-floor_u};
    };
    const auto [ix0,fx]=clamp_cell((position.x-domain.minimum_m.x)/dx,domain.cells_x);
    const auto [iy0,fy]=clamp_cell((position.y-domain.minimum_m.y)/dy,domain.cells_y);
    const auto [iz0,fz]=clamp_cell((position.z-domain.minimum_m.z)/dz,domain.cells_z);
    const std::size_t ix1=std::min(ix0+1U,domain.cells_x-1U);
    const std::size_t iy1=std::min(iy0+1U,domain.cells_y-1U);
    const std::size_t iz1=std::min(iz0+1U,domain.cells_z-1U);
    double value=0.0;
    for(std::size_t dz_i=0;dz_i<2U;++dz_i){const std::size_t iz=dz_i?iz1:iz0;const double wz=dz_i?fz:(1.0-fz);
        for(std::size_t dy_i=0;dy_i<2U;++dy_i){const std::size_t iy=dy_i?iy1:iy0;const double wy=dy_i?fy:(1.0-fy);
            for(std::size_t dx_i=0;dx_i<2U;++dx_i){const std::size_t ix=dx_i?ix1:ix0;const double wx=dx_i?fx:(1.0-fx);
                value+=values[flat_index_3d(ix,iy,iz,domain.cells_x,domain.cells_y)]*wx*wy*wz;
            }
        }
    }
    return value;
}
ElectromagneticField local_domain_field_3d(const ElectromagneticFieldBlocks3D& fields,const PicDomain3D& domain,Vec3 position){
    return {{local_sample_field_3d(fields.electric_x_by_domain[domain.domain_id],domain,position),
             local_sample_field_3d(fields.electric_y_by_domain[domain.domain_id],domain,position),
             local_sample_field_3d(fields.electric_z_by_domain[domain.domain_id],domain,position)},
            {local_sample_field_3d(fields.magnetic_x_by_domain[domain.domain_id],domain,position),
             local_sample_field_3d(fields.magnetic_y_by_domain[domain.domain_id],domain,position),
             local_sample_field_3d(fields.magnetic_z_by_domain[domain.domain_id],domain,position)}};
}
}

DistributedStaggeredPicStep3D run_serialized_distributed_staggered_pic_step_3d(
    const std::vector<std::vector<PicParticle3D>>& particles_by_domain,
    const ElectromagneticFieldBlocks3D& fields_by_domain,
    std::span<const PicDomain3D> domains,const PicDomainGrid3DConfig& config,
    const DistributedStaggeredPicStep3DConfig& step_config,
    const PicRankTopology3D& topology){
    validate_domain_span(domains,config);
    if(!(step_config.dt_s>0.0)||!std::isfinite(step_config.dt_s))throw std::invalid_argument("distributed PIC step requires positive dt");
    if(particles_by_domain.size()!=domains.size())throw std::invalid_argument("distributed PIC step particle domain count mismatch");
    if(topology.domain_to_rank.size()!=domains.size())throw std::invalid_argument("distributed PIC step topology/domain mismatch");
    validate_em_domain_fields_3d(fields_by_domain,domains);
    auto advanced=particles_by_domain;
    DistributedStaggeredPicStep3D result;result.particle_count_before=count_pic_particles_by_domain_3d(particles_by_domain);
    for(const auto& domain:domains){
        auto& bucket=advanced[domain.domain_id];
        for(auto& particle:bucket){
            const Vec3 old_position=particle.position_m;
            ChargedParticle charged;charged.position_m=particle.position_m;charged.velocity_m_per_s=particle.velocity_m_per_s;charged.charge_c=particle.charge_c;charged.mass_kg=particle.mass_kg;charged.weight=particle.weight;
            vay_push(charged,local_domain_field_3d(fields_by_domain,domain,particle.position_m),step_config.dt_s);
            particle.position_m=charged.position_m;particle.velocity_m_per_s=charged.velocity_m_per_s;
            const Vec3 delta{particle.position_m.x-old_position.x,particle.position_m.y-old_position.y,particle.position_m.z-old_position.z};
            result.max_particle_displacement_m=std::max(result.max_particle_displacement_m,std::sqrt(norm2(delta)));
        }
    }
    result.particle_count_after_push=count_pic_particles_by_domain_3d(advanced);
    result.guarded_fields=step_config.exchange_em_field_guards
        ?exchange_electromagnetic_field_guard_cells_3d(fields_by_domain,domains,config,step_config.guard_cells,step_config.exterior_value)
        :ElectromagneticGuardedFieldBlocks3D{};
    result.exchange_round=run_serialized_distributed_pic_exchange_round_3d(
        advanced,fields_by_domain.electric_x_by_domain,domains,config,step_config.guard_cells,step_config.exterior_value,topology);
    result.particles_by_domain=result.exchange_round.particles_by_domain;
    result.particle_count_after_migration=count_pic_particles_by_domain_3d(result.particles_by_domain);
    return result;
}

#if defined(CFD_HAS_MPI)
namespace {
void check_pic_mpi_3d(int code,const char* operation){
    if(code==MPI_SUCCESS)return;char message[MPI_MAX_ERROR_STRING]{};int length=0;MPI_Error_string(code,message,&length);throw std::runtime_error(std::string(operation)+": "+std::string(message,static_cast<std::size_t>(length)));
}
int mpi_byte_count_3d(std::size_t n){if(n>static_cast<std::size_t>(std::numeric_limits<int>::max()))throw std::overflow_error("PIC MPI payload exceeds INT_MAX");return static_cast<int>(n);} 
}
std::vector<PicSerializedEnvelope3D> exchange_pic_serialized_envelopes_mpi_3d(std::span<const PicSerializedEnvelope3D> outbound,MPI_Comm communicator){
    int rank=0,size=1;check_pic_mpi_3d(MPI_Comm_rank(communicator,&rank),"MPI_Comm_rank");check_pic_mpi_3d(MPI_Comm_size(communicator,&size),"MPI_Comm_size");
    ByteWriter3D writer;writer.write(static_cast<std::size_t>(outbound.size()));for(const auto& envelope:outbound){writer.write(envelope.bytes.size());writer.write_bytes(envelope.bytes);}auto local=std::move(writer).finish().bytes;
    std::vector<int> counts(static_cast<std::size_t>(size));const int local_count=mpi_byte_count_3d(local.size());check_pic_mpi_3d(MPI_Allgather(&local_count,1,MPI_INT,counts.data(),1,MPI_INT,communicator),"MPI_Allgather");
    std::vector<int> displacements(static_cast<std::size_t>(size));int total=0;for(int i=0;i<size;++i){displacements[static_cast<std::size_t>(i)]=total;total+=counts[static_cast<std::size_t>(i)];}
    std::vector<std::uint8_t> global(static_cast<std::size_t>(total));check_pic_mpi_3d(MPI_Allgatherv(local.data(),local_count,MPI_UNSIGNED_CHAR,global.data(),counts.data(),displacements.data(),MPI_UNSIGNED_CHAR,communicator),"MPI_Allgatherv");
    std::vector<PicSerializedEnvelope3D> received;
    for(int source=0;source<size;++source){const auto begin=global.data()+displacements[static_cast<std::size_t>(source)];const auto count=static_cast<std::size_t>(counts[static_cast<std::size_t>(source)]);ByteReader3D reader(std::span<const std::uint8_t>(begin,count));const auto records=reader.read<std::size_t>();for(std::size_t record=0;record<records;++record){const auto nbytes=reader.read<std::size_t>();if(nbytes>count)throw std::invalid_argument("PIC MPI serialized record length out of range");PicSerializedEnvelope3D payload;payload.bytes.resize(nbytes);for(std::size_t i=0;i<nbytes;++i)payload.bytes[i]=reader.read<std::uint8_t>();const auto decoded=deserialize_pic_transport_envelope_3d(payload);if(decoded.destination_rank==static_cast<std::size_t>(rank))received.push_back(std::move(payload));}reader.require_finished();}
    return received;
}
#endif

std::vector<double> deposit_cic_charge_density_3d(std::span<const PicParticle3D> particles,std::size_t nx,std::size_t ny,std::size_t nz,
                                                   double length_x_m,double length_y_m,double length_z_m){
    if(nx<2U||ny<2U||nz<2U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(length_z_m>0.0))throw std::invalid_argument("invalid 3-D CIC grid");
    std::vector<double> rho(nx*ny*nz,0.0);const double dx_value=length_x_m/static_cast<double>(nx);const double dy_value=length_y_m/static_cast<double>(ny);const double dz_value=length_z_m/static_cast<double>(nz);
    for(const auto& particle:particles){
        if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid 3-D PIC particle");
        const double ux=wrap_periodic_value(particle.position_m.x,length_x_m)/dx_value;const double uy=wrap_periodic_value(particle.position_m.y,length_y_m)/dy_value;const double uz=wrap_periodic_value(particle.position_m.z,length_z_m)/dz_value;
        const auto ix0=static_cast<std::size_t>(std::floor(ux))%nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%ny;const auto iz0=static_cast<std::size_t>(std::floor(uz))%nz;
        const auto ix1=(ix0+1U)%nx;const auto iy1=(iy0+1U)%ny;const auto iz1=(iz0+1U)%nz;const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);const double fz=uz-std::floor(uz);
        const double density=particle.charge_c*particle.weight/(dx_value*dy_value*dz_value);
        for(std::size_t dz_i=0;dz_i<2U;++dz_i){
            const std::size_t iz=dz_i?iz1:iz0;const double wz=dz_i?fz:(1.0-fz);
            for(std::size_t dy_i=0;dy_i<2U;++dy_i){
                const std::size_t iy=dy_i?iy1:iy0;const double wy=dy_i?fy:(1.0-fy);
                for(std::size_t dx_i=0;dx_i<2U;++dx_i){
                    const std::size_t ix=dx_i?ix1:ix0;const double wx=dx_i?fx:(1.0-fx);
                    rho[flat_index_3d(ix,iy,iz,nx,ny)]+=density*wx*wy*wz;
                }
            }
        }
    }
    return rho;
}

ElectrostaticField3D periodic_electric_field_from_charge_density_3d(std::span<const double> rho,std::size_t nx,std::size_t ny,std::size_t nz,
                                                                     double length_x_m,double length_y_m,double length_z_m,double relative_permittivity,
                                                                     bool remove_mean){
    if(rho.size()!=nx*ny*nz||nx<2U||ny<2U||nz<2U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(length_z_m>0.0)||!(relative_permittivity>0.0))throw std::invalid_argument("invalid 3-D periodic Poisson grid");
    double mean_value=0.0,max_abs_value=0.0;for(double value:rho){if(!std::isfinite(value))throw std::invalid_argument("non-finite 3-D charge density");mean_value+=value;max_abs_value=std::max(max_abs_value,std::abs(value));}mean_value/=static_cast<double>(rho.size());
    if(!remove_mean&&std::abs(mean_value)>1.0e-12*std::max(1.0,max_abs_value))throw std::invalid_argument("3-D periodic Poisson requires mean-zero charge density");
    std::vector<double> centered(rho.begin(),rho.end());if(remove_mean)for(double& value:centered)value-=mean_value;
    using Complex=std::complex<double>;std::vector<Complex> ex_hat(nx*ny*nz),ey_hat(nx*ny*nz),ez_hat(nx*ny*nz);const double epsilon=epsilon0*relative_permittivity;const double two_pi=2.0*std::numbers::pi;
    for(std::size_t mz=0;mz<nz;++mz){
        for(std::size_t my=0;my<ny;++my){
            for(std::size_t mx=0;mx<nx;++mx){
                if(mx==0U&&my==0U&&mz==0U)continue;
                const auto modes=signed_modes_3d(mx,my,mz,nx,ny,nz);const double kx=two_pi*static_cast<double>(modes[0])/length_x_m;const double ky=two_pi*static_cast<double>(modes[1])/length_y_m;const double kz=two_pi*static_cast<double>(modes[2])/length_z_m;const double k2=kx*kx+ky*ky+kz*kz;if(k2==0.0)continue;
                const Complex rho_hat=fourier_coefficient_3d(centered,nx,ny,nz,mx,my,mz);const Complex factor=Complex{0.0,-1.0}*rho_hat/(epsilon*k2);const std::size_t k_index=flat_index_3d(mx,my,mz,nx,ny);
                ex_hat[k_index]=factor*kx;ey_hat[k_index]=factor*ky;ez_hat[k_index]=factor*kz;
            }
        }
    }
    ElectrostaticField3D field;field.electric_x_v_per_m.assign(nx*ny*nz,0.0);field.electric_y_v_per_m.assign(nx*ny*nz,0.0);field.electric_z_v_per_m.assign(nx*ny*nz,0.0);
    for(std::size_t iz=0;iz<nz;++iz){
        for(std::size_t iy=0;iy<ny;++iy){
            for(std::size_t ix=0;ix<nx;++ix){
                Complex ex_value{},ey_value{},ez_value{};
                for(std::size_t mz=0;mz<nz;++mz){
                    for(std::size_t my=0;my<ny;++my){
                        for(std::size_t mx=0;mx<nx;++mx){
                            const double angle=two_pi*(static_cast<double>(mx*ix)/static_cast<double>(nx)+static_cast<double>(my*iy)/static_cast<double>(ny)+static_cast<double>(mz*iz)/static_cast<double>(nz));const Complex phase{std::cos(angle),std::sin(angle)};const std::size_t k_index=flat_index_3d(mx,my,mz,nx,ny);
                            ex_value+=ex_hat[k_index]*phase;ey_value+=ey_hat[k_index]*phase;ez_value+=ez_hat[k_index]*phase;
                        }
                    }
                }
                const std::size_t idx=flat_index_3d(ix,iy,iz,nx,ny);field.electric_x_v_per_m[idx]=ex_value.real();field.electric_y_v_per_m[idx]=ey_value.real();field.electric_z_v_per_m[idx]=ez_value.real();
            }
        }
    }
    return field;
}

ChargeConservingCurrent3D deposit_charge_conserving_current_3d(std::span<const PicParticle3D> particles_old,
                                                                std::span<const PicParticle3D> particles_new,
                                                                std::size_t nx,std::size_t ny,std::size_t nz,
                                                                double length_x_m,double length_y_m,double length_z_m,double dt_s){
    if(particles_old.size()!=particles_new.size())throw std::invalid_argument("old/new 3-D particle arrays must have equal size");
    if(nx<2U||ny<2U||nz<2U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(length_z_m>0.0)||!(dt_s>0.0)||!std::isfinite(dt_s))throw std::invalid_argument("invalid 3-D charge-conserving deposition parameters");
    ChargeConservingCurrent3D result;result.charge_density_old_c_per_m3=deposit_cic_charge_density_3d(particles_old,nx,ny,nz,length_x_m,length_y_m,length_z_m);result.charge_density_new_c_per_m3=deposit_cic_charge_density_3d(particles_new,nx,ny,nz,length_x_m,length_y_m,length_z_m);result.current_x_a_per_m2.assign(nx*ny*nz,0.0);result.current_y_a_per_m2.assign(nx*ny*nz,0.0);result.current_z_a_per_m2.assign(nx*ny*nz,0.0);
    std::vector<double> delta(nx*ny*nz);for(std::size_t idx=0;idx<delta.size();++idx)delta[idx]=result.charge_density_new_c_per_m3[idx]-result.charge_density_old_c_per_m3[idx];
    using Complex=std::complex<double>;std::vector<Complex> jx_hat(nx*ny*nz),jy_hat(nx*ny*nz),jz_hat(nx*ny*nz),div_hat(nx*ny*nz);const double two_pi=2.0*std::numbers::pi;
    for(std::size_t mz=0;mz<nz;++mz){
        for(std::size_t my=0;my<ny;++my){
            for(std::size_t mx=0;mx<nx;++mx){
                if(mx==0U&&my==0U&&mz==0U)continue;
                const auto modes=signed_modes_3d(mx,my,mz,nx,ny,nz);const double kx=two_pi*static_cast<double>(modes[0])/length_x_m;const double ky=two_pi*static_cast<double>(modes[1])/length_y_m;const double kz=two_pi*static_cast<double>(modes[2])/length_z_m;const double k2=kx*kx+ky*ky+kz*kz;if(k2==0.0)continue;
                const Complex delta_hat=fourier_coefficient_3d(delta,nx,ny,nz,mx,my,mz);const Complex factor=Complex{0.0,1.0}*delta_hat/(dt_s*k2);const std::size_t k_index=flat_index_3d(mx,my,mz,nx,ny);
                jx_hat[k_index]=factor*kx;jy_hat[k_index]=factor*ky;jz_hat[k_index]=factor*kz;div_hat[k_index]=Complex{0.0,1.0}*(kx*jx_hat[k_index]+ky*jy_hat[k_index]+kz*jz_hat[k_index]);
            }
        }
    }
    std::vector<double> divergence(nx*ny*nz,0.0);
    for(std::size_t iz=0;iz<nz;++iz){
        for(std::size_t iy=0;iy<ny;++iy){
            for(std::size_t ix=0;ix<nx;++ix){
                Complex jx_value{},jy_value{},jz_value{},div_value{};
                for(std::size_t mz=0;mz<nz;++mz){
                    for(std::size_t my=0;my<ny;++my){
                        for(std::size_t mx=0;mx<nx;++mx){
                            const double angle=two_pi*(static_cast<double>(mx*ix)/static_cast<double>(nx)+static_cast<double>(my*iy)/static_cast<double>(ny)+static_cast<double>(mz*iz)/static_cast<double>(nz));const Complex phase{std::cos(angle),std::sin(angle)};const std::size_t k_index=flat_index_3d(mx,my,mz,nx,ny);
                            jx_value+=jx_hat[k_index]*phase;jy_value+=jy_hat[k_index]*phase;jz_value+=jz_hat[k_index]*phase;div_value+=div_hat[k_index]*phase;
                        }
                    }
                }
                const std::size_t idx=flat_index_3d(ix,iy,iz,nx,ny);result.current_x_a_per_m2[idx]=jx_value.real();result.current_y_a_per_m2[idx]=jy_value.real();result.current_z_a_per_m2[idx]=jz_value.real();divergence[idx]=div_value.real();
            }
        }
    }
    double residual=0.0;for(std::size_t idx=0;idx<delta.size();++idx){const double continuity=delta[idx]/dt_s+divergence[idx];residual=std::max(residual,std::abs(continuity));}result.continuity_linf_residual=residual;return result;
}


ChargeConservingCurrent3D deposit_charge_conserving_current_3d_local(std::span<const PicParticle3D> particles_old,
                                                                      std::span<const PicParticle3D> particles_new,
                                                                      std::size_t nx,std::size_t ny,std::size_t nz,
                                                                      double length_x_m,double length_y_m,double length_z_m,double dt_s){
    if(particles_old.size()!=particles_new.size())throw std::invalid_argument("old/new 3-D particle arrays must have equal size");
    if(nx<2U||ny<2U||nz<2U||!(length_x_m>0.0)||!(length_y_m>0.0)||!(length_z_m>0.0)||!(dt_s>0.0)||!std::isfinite(dt_s))throw std::invalid_argument("invalid local 3-D current deposition parameters");
    ChargeConservingCurrent3D result;result.charge_density_old_c_per_m3=deposit_cic_charge_density_3d(particles_old,nx,ny,nz,length_x_m,length_y_m,length_z_m);result.charge_density_new_c_per_m3=deposit_cic_charge_density_3d(particles_new,nx,ny,nz,length_x_m,length_y_m,length_z_m);
    const std::size_t n=nx*ny*nz;result.current_x_a_per_m2.assign(n,0.0);result.current_y_a_per_m2.assign(n,0.0);result.current_z_a_per_m2.assign(n,0.0);
    std::vector<double> delta(n);for(std::size_t idx=0;idx<n;++idx)delta[idx]=result.charge_density_new_c_per_m3[idx]-result.charge_density_old_c_per_m3[idx];
    const double dx_value=length_x_m/static_cast<double>(nx),dy_value=length_y_m/static_cast<double>(ny),dz_value=length_z_m/static_cast<double>(nz);
    std::vector<double> mean_x(ny*nz,0.0),mean_xy(nz,0.0);
    for(std::size_t iz=0;iz<nz;++iz){
        for(std::size_t iy=0;iy<ny;++iy){
            double sum=0.0;for(std::size_t ix=0;ix<nx;++ix)sum+=delta[flat_index_3d(ix,iy,iz,nx,ny)];
            mean_x[iz*ny+iy]=sum/static_cast<double>(nx);
        }
    }
    for(std::size_t iz=0;iz<nz;++iz){
        double sum=0.0;for(std::size_t iy=0;iy<ny;++iy)sum+=mean_x[iz*ny+iy];
        mean_xy[iz]=sum/static_cast<double>(ny);
    }
    double global_mean=0.0;for(double value:mean_xy)global_mean+=value;global_mean/=static_cast<double>(nz);

    for(std::size_t iz=0;iz<nz;++iz){
        for(std::size_t iy=0;iy<ny;++iy){
            double cumulative=0.0;
            const double row_mean=mean_x[iz*ny+iy];
            for(std::size_t ix=0;ix<nx;++ix){
                const std::size_t idx=flat_index_3d(ix,iy,iz,nx,ny);
                cumulative+=-dx_value*(delta[idx]-row_mean)/dt_s;
                result.current_x_a_per_m2[idx]=cumulative;
            }
        }
    }
    for(std::size_t iz=0;iz<nz;++iz){
        double cumulative=0.0;
        for(std::size_t iy=0;iy<ny;++iy){
            cumulative+=-dy_value*(mean_x[iz*ny+iy]-mean_xy[iz])/dt_s;
            for(std::size_t ix=0;ix<nx;++ix)result.current_y_a_per_m2[flat_index_3d(ix,iy,iz,nx,ny)]=cumulative;
        }
    }
    double cumulative_z=0.0;
    for(std::size_t iz=0;iz<nz;++iz){
        cumulative_z+=-dz_value*(mean_xy[iz]-global_mean)/dt_s;
        for(std::size_t iy=0;iy<ny;++iy)for(std::size_t ix=0;ix<nx;++ix)result.current_z_a_per_m2[flat_index_3d(ix,iy,iz,nx,ny)]=cumulative_z;
    }

    double residual=0.0;
    for(std::size_t iz=0;iz<nz;++iz){
        const std::size_t iz_back=(iz+nz-1U)%nz;
        for(std::size_t iy=0;iy<ny;++iy){
            const std::size_t iy_back=(iy+ny-1U)%ny;
            for(std::size_t ix=0;ix<nx;++ix){
                const std::size_t ix_back=(ix+nx-1U)%nx;
                const std::size_t idx=flat_index_3d(ix,iy,iz,nx,ny);
                const double divergence=(result.current_x_a_per_m2[idx]-result.current_x_a_per_m2[flat_index_3d(ix_back,iy,iz,nx,ny)])/dx_value
                    +(result.current_y_a_per_m2[idx]-result.current_y_a_per_m2[flat_index_3d(ix,iy_back,iz,nx,ny)])/dy_value
                    +(result.current_z_a_per_m2[idx]-result.current_z_a_per_m2[flat_index_3d(ix,iy,iz_back,nx,ny)])/dz_value;
                const double continuity=delta[idx]/dt_s+divergence;
                residual=std::max(residual,std::abs(continuity));
            }
        }
    }
    result.continuity_linf_residual=residual;return result;
}

ElectrostaticPic3D::ElectrostaticPic3D(ElectrostaticPic3DConfig config):config_(config){
    if(config_.nx<2U||config_.ny<2U||config_.nz<2U||!(config_.length_x_m>0.0)||!(config_.length_y_m>0.0)||!(config_.length_z_m>0.0)||!(config_.relative_permittivity>0.0)||!(config_.dt_s>0.0)||!std::isfinite(config_.dt_s))throw std::invalid_argument("invalid electrostatic PIC 3-D configuration");
    const std::size_t n=config_.nx*config_.ny*config_.nz;charge_density_.assign(n,0.0);electric_x_.assign(n,0.0);electric_y_.assign(n,0.0);electric_z_.assign(n,0.0);
}

Vec3 ElectrostaticPic3D::wrap_position(Vec3 position_m) const {position_m.x=wrap_periodic_value(position_m.x,config_.length_x_m);position_m.y=wrap_periodic_value(position_m.y,config_.length_y_m);position_m.z=wrap_periodic_value(position_m.z,config_.length_z_m);return position_m;}

void ElectrostaticPic3D::set_particles(std::vector<PicParticle3D> particles){for(auto& particle:particles){if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid electrostatic PIC 3-D particle");particle.position_m=wrap_position(particle.position_m);}particles_=std::move(particles);}

void ElectrostaticPic3D::deposit_and_solve(){charge_density_=deposit_cic_charge_density_3d(particles_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m);const auto field=periodic_electric_field_from_charge_density_3d(charge_density_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m,config_.relative_permittivity,config_.neutralize_mean_charge);electric_x_=field.electric_x_v_per_m;electric_y_=field.electric_y_v_per_m;electric_z_=field.electric_z_v_per_m;}

double ElectrostaticPic3D::trilinear(std::span<const double> nodal,Vec3 position_m) const {
    const Vec3 wrapped=wrap_position(position_m);const double ux=wrapped.x/dx();const double uy=wrapped.y/dy();const double uz=wrapped.z/dz();const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;const auto iz0=static_cast<std::size_t>(std::floor(uz))%config_.nz;const auto ix1=(ix0+1U)%config_.nx;const auto iy1=(iy0+1U)%config_.ny;const auto iz1=(iz0+1U)%config_.nz;const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);const double fz=uz-std::floor(uz);
    double value=0.0;for(std::size_t dz_i=0;dz_i<2U;++dz_i){const std::size_t iz=dz_i?iz1:iz0;const double wz=dz_i?fz:(1.0-fz);for(std::size_t dy_i=0;dy_i<2U;++dy_i){const std::size_t iy=dy_i?iy1:iy0;const double wy=dy_i?fy:(1.0-fy);for(std::size_t dx_i=0;dx_i<2U;++dx_i){const std::size_t ix=dx_i?ix1:ix0;const double wx=dx_i?fx:(1.0-fx);value+=nodal[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)]*wx*wy*wz;}}}return value;
}

Vec3 ElectrostaticPic3D::gather_electric_field(Vec3 position_m) const {return {trilinear(electric_x_,position_m),trilinear(electric_y_,position_m),trilinear(electric_z_,position_m)};}

void ElectrostaticPic3D::step(std::size_t steps){for(std::size_t step_index=0;step_index<steps;++step_index){deposit_and_solve();for(auto& particle:particles_){const Vec3 field=gather_electric_field(particle.position_m);particle.velocity_m_per_s.x+=(particle.charge_c/particle.mass_kg)*field.x*config_.dt_s;particle.velocity_m_per_s.y+=(particle.charge_c/particle.mass_kg)*field.y*config_.dt_s;particle.velocity_m_per_s.z+=(particle.charge_c/particle.mass_kg)*field.z*config_.dt_s;particle.position_m.x+=particle.velocity_m_per_s.x*config_.dt_s;particle.position_m.y+=particle.velocity_m_per_s.y*config_.dt_s;particle.position_m.z+=particle.velocity_m_per_s.z*config_.dt_s;particle.position_m=wrap_position(particle.position_m);}time_s_+=config_.dt_s;}deposit_and_solve();}


ElectromagneticPic3D::ElectromagneticPic3D(ElectromagneticPic3DConfig config):config_(config){
    if(config_.nx<2U||config_.ny<2U||config_.nz<2U||!(config_.length_x_m>0.0)||!(config_.length_y_m>0.0)||!(config_.length_z_m>0.0)||!(config_.relative_permittivity>0.0)||!(config_.relative_permeability>0.0)||!(config_.dt_s>0.0)||!std::isfinite(config_.dt_s))throw std::invalid_argument("invalid electromagnetic PIC 3-D configuration");
    const double c=medium_light_speed();
    const double inv_cfl=std::sqrt(1.0/(dx()*dx())+1.0/(dy()*dy())+1.0/(dz()*dz()));
    if(config_.dt_s*c*inv_cfl>=0.55)throw std::invalid_argument("3-D EM-PIC time step violates compact centered-CFL guard");
    const std::size_t n=config_.nx*config_.ny*config_.nz;
    charge_density_.assign(n,0.0);current_x_.assign(n,0.0);current_y_.assign(n,0.0);current_z_.assign(n,0.0);
    electric_x_.assign(n,0.0);electric_y_.assign(n,0.0);electric_z_.assign(n,0.0);
    magnetic_x_.assign(n,0.0);magnetic_y_.assign(n,0.0);magnetic_z_.assign(n,0.0);
}

double ElectromagneticPic3D::epsilon() const noexcept {return epsilon0*config_.relative_permittivity;}
double ElectromagneticPic3D::mu() const noexcept {return mu0*config_.relative_permeability;}
double ElectromagneticPic3D::medium_light_speed() const noexcept {return 1.0/std::sqrt(epsilon()*mu());}

Vec3 ElectromagneticPic3D::wrap_position(Vec3 position_m) const {
    position_m.x=wrap_periodic_value(position_m.x,config_.length_x_m);
    position_m.y=wrap_periodic_value(position_m.y,config_.length_y_m);
    position_m.z=wrap_periodic_value(position_m.z,config_.length_z_m);
    return position_m;
}

void ElectromagneticPic3D::set_particles(std::vector<PicParticle3D> particles){
    for(auto& particle:particles){
        if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid electromagnetic PIC 3-D particle");
        particle.position_m=wrap_position(particle.position_m);
    }
    particles_=std::move(particles);
}

void ElectromagneticPic3D::set_fields(std::span<const double> electric_x_v_per_m,
                                      std::span<const double> electric_y_v_per_m,
                                      std::span<const double> electric_z_v_per_m,
                                      std::span<const double> magnetic_x_t,
                                      std::span<const double> magnetic_y_t,
                                      std::span<const double> magnetic_z_t){
    const std::size_t n=config_.nx*config_.ny*config_.nz;
    if(electric_x_v_per_m.size()!=n||electric_y_v_per_m.size()!=n||electric_z_v_per_m.size()!=n||magnetic_x_t.size()!=n||magnetic_y_t.size()!=n||magnetic_z_t.size()!=n)throw std::invalid_argument("3-D EM-PIC field arrays must match grid size");
    electric_x_.assign(electric_x_v_per_m.begin(),electric_x_v_per_m.end());
    electric_y_.assign(electric_y_v_per_m.begin(),electric_y_v_per_m.end());
    electric_z_.assign(electric_z_v_per_m.begin(),electric_z_v_per_m.end());
    magnetic_x_.assign(magnetic_x_t.begin(),magnetic_x_t.end());
    magnetic_y_.assign(magnetic_y_t.begin(),magnetic_y_t.end());
    magnetic_z_.assign(magnetic_z_t.begin(),magnetic_z_t.end());
}

void ElectromagneticPic3D::initialize_z_polarized_mode(double electric_z_amplitude_v_per_m,std::size_t mode_x,std::size_t mode_y,std::size_t mode_z){
    if(mode_x==0U&&mode_y==0U&&mode_z==0U)throw std::invalid_argument("3-D EM-PIC mode cannot be DC");
    std::fill(electric_x_.begin(),electric_x_.end(),0.0);std::fill(electric_y_.begin(),electric_y_.end(),0.0);std::fill(electric_z_.begin(),electric_z_.end(),0.0);
    std::fill(magnetic_x_.begin(),magnetic_x_.end(),0.0);std::fill(magnetic_y_.begin(),magnetic_y_.end(),0.0);std::fill(magnetic_z_.begin(),magnetic_z_.end(),0.0);
    const double kx=2.0*std::numbers::pi*static_cast<double>(mode_x)/config_.length_x_m;
    const double ky=2.0*std::numbers::pi*static_cast<double>(mode_y)/config_.length_y_m;
    const double kz=2.0*std::numbers::pi*static_cast<double>(mode_z)/config_.length_z_m;
    const double k=std::sqrt(kx*kx+ky*ky+kz*kz);
    for(std::size_t iz=0;iz<config_.nz;++iz){
        const double z=config_.length_z_m*static_cast<double>(iz)/static_cast<double>(config_.nz);
        for(std::size_t iy=0;iy<config_.ny;++iy){
            const double y=config_.length_y_m*static_cast<double>(iy)/static_cast<double>(config_.ny);
            for(std::size_t ix=0;ix<config_.nx;++ix){
                const double x=config_.length_x_m*static_cast<double>(ix)/static_cast<double>(config_.nx);
                const double phase=kx*x+ky*y+kz*z;
                const double ez=electric_z_amplitude_v_per_m*std::sin(phase);
                const double cphase=std::sin(phase);
                const std::size_t idx=flat_index_3d(ix,iy,iz,config_.nx,config_.ny);
                electric_z_[idx]=ez;
                // B = (k-hat x E)/c = (mu/Z)*(k-hat x E). For z-polarized E this creates transverse B.
                magnetic_x_[idx]= (ky/k)*(electric_z_amplitude_v_per_m/medium_light_speed())*cphase;
                magnetic_y_[idx]=-(kx/k)*(electric_z_amplitude_v_per_m/medium_light_speed())*cphase;
            }
        }
    }
}

void ElectromagneticPic3D::deposit_sources(){
    charge_density_=deposit_cic_charge_density_3d(particles_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m);
    if(config_.solve_longitudinal_poisson){
        const auto field=periodic_electric_field_from_charge_density_3d(charge_density_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m,config_.relative_permittivity,config_.neutralize_mean_charge);
        electric_x_=field.electric_x_v_per_m;electric_y_=field.electric_y_v_per_m;electric_z_=field.electric_z_v_per_m;
    }
}

double ElectromagneticPic3D::trilinear(std::span<const double> nodal,Vec3 position_m) const {
    const Vec3 wrapped=wrap_position(position_m);const double ux=wrapped.x/dx();const double uy=wrapped.y/dy();const double uz=wrapped.z/dz();
    const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;const auto iz0=static_cast<std::size_t>(std::floor(uz))%config_.nz;
    const auto ix1=(ix0+1U)%config_.nx;const auto iy1=(iy0+1U)%config_.ny;const auto iz1=(iz0+1U)%config_.nz;const double fx=ux-std::floor(ux);const double fy=uy-std::floor(uy);const double fz=uz-std::floor(uz);
    double value=0.0;
    for(std::size_t dz_i=0;dz_i<2U;++dz_i){const std::size_t iz=dz_i?iz1:iz0;const double wz=dz_i?fz:(1.0-fz);
        for(std::size_t dy_i=0;dy_i<2U;++dy_i){const std::size_t iy=dy_i?iy1:iy0;const double wy=dy_i?fy:(1.0-fy);
            for(std::size_t dx_i=0;dx_i<2U;++dx_i){const std::size_t ix=dx_i?ix1:ix0;const double wx=dx_i?fx:(1.0-fx);value+=nodal[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)]*wx*wy*wz;}}}
    return value;
}

ElectromagneticField ElectromagneticPic3D::gather_field(Vec3 position_m) const {
    return {{trilinear(electric_x_,position_m),trilinear(electric_y_,position_m),trilinear(electric_z_,position_m)},
            {trilinear(magnetic_x_,position_m),trilinear(magnetic_y_,position_m),trilinear(magnetic_z_,position_m)}};
}

double ElectromagneticPic3D::ddx(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t r=(ix+1U)%config_.nx,l=(ix+config_.nx-1U)%config_.nx;return (values[flat_index_3d(r,iy,iz,config_.nx,config_.ny)]-values[flat_index_3d(l,iy,iz,config_.nx,config_.ny)])/(2.0*dx());}
double ElectromagneticPic3D::ddy(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t u=(iy+1U)%config_.ny,d=(iy+config_.ny-1U)%config_.ny;return (values[flat_index_3d(ix,u,iz,config_.nx,config_.ny)]-values[flat_index_3d(ix,d,iz,config_.nx,config_.ny)])/(2.0*dy());}
double ElectromagneticPic3D::ddz(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t f=(iz+1U)%config_.nz,b=(iz+config_.nz-1U)%config_.nz;return (values[flat_index_3d(ix,iy,f,config_.nx,config_.ny)]-values[flat_index_3d(ix,iy,b,config_.nx,config_.ny)])/(2.0*dz());}

void ElectromagneticPic3D::advance_magnetic(double dt_s){
    auto bx=magnetic_x_,by=magnetic_y_,bz=magnetic_z_;
    for(std::size_t iz=0;iz<config_.nz;++iz){
        for(std::size_t iy=0;iy<config_.ny;++iy){
            for(std::size_t ix=0;ix<config_.nx;++ix){
                const std::size_t idx=flat_index_3d(ix,iy,iz,config_.nx,config_.ny);
                const double curl_e_x=ddy(electric_z_,ix,iy,iz)-ddz(electric_y_,ix,iy,iz);
                const double curl_e_y=ddz(electric_x_,ix,iy,iz)-ddx(electric_z_,ix,iy,iz);
                const double curl_e_z=ddx(electric_y_,ix,iy,iz)-ddy(electric_x_,ix,iy,iz);
                bx[idx]-=dt_s*curl_e_x;by[idx]-=dt_s*curl_e_y;bz[idx]-=dt_s*curl_e_z;
            }
        }
    }
    magnetic_x_.swap(bx);magnetic_y_.swap(by);magnetic_z_.swap(bz);
}

void ElectromagneticPic3D::update_electric(double dt_s){
    auto ex=electric_x_,ey=electric_y_,ez=electric_z_;const double eps=epsilon();const double permeability=mu();
    for(std::size_t iz=0;iz<config_.nz;++iz){
        for(std::size_t iy=0;iy<config_.ny;++iy){
            for(std::size_t ix=0;ix<config_.nx;++ix){
                const std::size_t idx=flat_index_3d(ix,iy,iz,config_.nx,config_.ny);
                const double curl_h_x=(ddy(magnetic_z_,ix,iy,iz)-ddz(magnetic_y_,ix,iy,iz))/permeability;
                const double curl_h_y=(ddz(magnetic_x_,ix,iy,iz)-ddx(magnetic_z_,ix,iy,iz))/permeability;
                const double curl_h_z=(ddx(magnetic_y_,ix,iy,iz)-ddy(magnetic_x_,ix,iy,iz))/permeability;
                ex[idx]+=dt_s*(curl_h_x-current_x_[idx])/eps;ey[idx]+=dt_s*(curl_h_y-current_y_[idx])/eps;ez[idx]+=dt_s*(curl_h_z-current_z_[idx])/eps;
            }
        }
    }
    electric_x_.swap(ex);electric_y_.swap(ey);electric_z_.swap(ez);
}

void ElectromagneticPic3D::step(std::size_t steps,const NeutralCollisionModel* collision_model){
    for(std::size_t step_index=0;step_index<steps;++step_index){
        const auto particles_old=particles_;
        advance_magnetic(0.5*config_.dt_s);
        for(auto& particle:particles_){
            ChargedParticle charged;charged.position_m=particle.position_m;charged.velocity_m_per_s=particle.velocity_m_per_s;charged.charge_c=particle.charge_c;charged.mass_kg=particle.mass_kg;charged.weight=particle.weight;
            vay_push(charged,gather_field(particle.position_m),config_.dt_s);
            particle.position_m=wrap_position(charged.position_m);particle.velocity_m_per_s=charged.velocity_m_per_s;
        }
        const auto current=deposit_charge_conserving_current_3d(particles_old,particles_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m,config_.dt_s);
        charge_density_=current.charge_density_new_c_per_m3;current_x_=current.current_x_a_per_m2;current_y_=current.current_y_a_per_m2;current_z_=current.current_z_a_per_m2;last_continuity_residual_=current.continuity_linf_residual;
        update_electric(config_.dt_s);
        advance_magnetic(0.5*config_.dt_s);
        if(config_.solve_longitudinal_poisson){const auto field=periodic_electric_field_from_charge_density_3d(charge_density_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m,config_.relative_permittivity,config_.neutralize_mean_charge);electric_x_=field.electric_x_v_per_m;electric_y_=field.electric_y_v_per_m;electric_z_=field.electric_z_v_per_m;}
        if(collision_model){std::vector<ChargedParticle> charged_particles;charged_particles.reserve(particles_.size());for(const auto& particle:particles_){charged_particles.push_back({particle.position_m,particle.velocity_m_per_s,particle.charge_c,particle.mass_kg,particle.weight});}const auto stats=apply_monte_carlo_collisions(charged_particles,config_.dt_s,*collision_model);accumulated_collisions_.elastic_events+=stats.elastic_events;accumulated_collisions_.ionization_events+=stats.ionization_events;accumulated_collisions_.energy_loss_j+=stats.energy_loss_j;particles_.clear();particles_.reserve(charged_particles.size());for(const auto& charged:charged_particles){PicParticle3D particle;particle.position_m=wrap_position(charged.position_m);particle.velocity_m_per_s=charged.velocity_m_per_s;particle.charge_c=charged.charge_c;particle.mass_kg=charged.mass_kg;particle.weight=charged.weight;particles_.push_back(particle);}}
        time_s_+=config_.dt_s;
    }
}

ElectromagneticPic3DDiagnostics ElectromagneticPic3D::diagnostics() const {
    ElectromagneticPic3DDiagnostics result;const double volume=dx()*dy()*dz();const double eps=epsilon();const double permeability=mu();
    for(std::size_t idx=0;idx<electric_x_.size();++idx){const double e2=electric_x_[idx]*electric_x_[idx]+electric_y_[idx]*electric_y_[idx]+electric_z_[idx]*electric_z_[idx];const double b2=magnetic_x_[idx]*magnetic_x_[idx]+magnetic_y_[idx]*magnetic_y_[idx]+magnetic_z_[idx]*magnetic_z_[idx];result.field_energy_j+=0.5*eps*e2*volume+0.5*b2*volume/permeability;}
    for(const auto& particle:particles_){ChargedParticle charged;charged.position_m=particle.position_m;charged.velocity_m_per_s=particle.velocity_m_per_s;charged.charge_c=particle.charge_c;charged.mass_kg=particle.mass_kg;charged.weight=particle.weight;result.particle_kinetic_energy_j+=particle.weight*kinetic_energy_j(charged);}    
    result.total_energy_j=result.field_energy_j+result.particle_kinetic_energy_j;result.charge_continuity_linf_residual=last_continuity_residual_;result.collision_statistics=accumulated_collisions_;result.particle_count=particles_.size();return result;
}

void apply_electromagnetic_field_boundary_3d(std::vector<double>& electric_x_v_per_m,
                                             std::vector<double>& electric_y_v_per_m,
                                             std::vector<double>& electric_z_v_per_m,
                                             std::vector<double>& magnetic_x_t,
                                             std::vector<double>& magnetic_y_t,
                                             std::vector<double>& magnetic_z_t,
                                             std::size_t nx,std::size_t ny,std::size_t nz,
                                             const GridBoundary3DConfig& boundary){
    if(electric_x_v_per_m.size()!=nx*ny*nz||electric_y_v_per_m.size()!=nx*ny*nz||electric_z_v_per_m.size()!=nx*ny*nz||
       magnetic_x_t.size()!=nx*ny*nz||magnetic_y_t.size()!=nx*ny*nz||magnetic_z_t.size()!=nx*ny*nz)throw std::invalid_argument("invalid 3-D field boundary array size");
    if(boundary.mode==GridBoundaryMode3D::periodic)return;
    if(boundary.mode==GridBoundaryMode3D::electric_wall){
        for(std::size_t iz=0;iz<nz;++iz){
            for(std::size_t iy=0;iy<ny;++iy){
                for(std::size_t ix=0;ix<nx;++ix){
                    if(ix==0U||iy==0U||iz==0U||ix+1U==nx||iy+1U==ny||iz+1U==nz){
                        const std::size_t idx=flat_index_3d(ix,iy,iz,nx,ny);
                        electric_x_v_per_m[idx]=0.0;electric_y_v_per_m[idx]=0.0;electric_z_v_per_m[idx]=0.0;
                    }
                }
            }
        }
        return;
    }
    const std::size_t cells=std::max<std::size_t>(1U,boundary.sponge_cells);
    const double strength=std::max(0.0,boundary.sponge_strength);
    const double polynomial_order=std::max(0.25,boundary.sponge_polynomial_order);
    for(std::size_t iz=0;iz<nz;++iz){
        const std::size_t dz_edge=std::min(iz,nz-1U-iz);
        for(std::size_t iy=0;iy<ny;++iy){
            const std::size_t dy_edge=std::min(iy,ny-1U-iy);
            for(std::size_t ix=0;ix<nx;++ix){
                const std::size_t dx_edge=std::min(ix,nx-1U-ix);
                const std::size_t edge=std::min(dx_edge,std::min(dy_edge,dz_edge));
                if(edge<cells){
                    const double s=static_cast<double>(cells-edge)/static_cast<double>(cells);
                    const double factor=std::exp(-strength*std::pow(s,polynomial_order));
                    const std::size_t idx=flat_index_3d(ix,iy,iz,nx,ny);
                    electric_x_v_per_m[idx]*=factor;electric_y_v_per_m[idx]*=factor;electric_z_v_per_m[idx]*=factor;
                    magnetic_x_t[idx]*=factor;magnetic_y_t[idx]*=factor;magnetic_z_t[idx]*=factor;
                }
            }
        }
    }
}

StaggeredElectromagneticPic3D::StaggeredElectromagneticPic3D(StaggeredElectromagneticPic3DConfig config):config_(config){
    if(config_.nx<2U||config_.ny<2U||config_.nz<2U||!(config_.length_x_m>0.0)||!(config_.length_y_m>0.0)||!(config_.length_z_m>0.0)||!(config_.relative_permittivity>0.0)||!(config_.relative_permeability>0.0)||!(config_.dt_s>0.0)||!std::isfinite(config_.dt_s))throw std::invalid_argument("invalid staggered electromagnetic PIC 3-D configuration");
    const double c=medium_light_speed();
    const double inv_cfl=std::sqrt(1.0/(dx()*dx())+1.0/(dy()*dy())+1.0/(dz()*dz()));
    if(config_.dt_s*c*inv_cfl>=0.98)throw std::invalid_argument("staggered 3-D EM-PIC time step violates Yee CFL guard");
    const std::size_t n=config_.nx*config_.ny*config_.nz;
    charge_density_.assign(n,0.0);current_x_.assign(n,0.0);current_y_.assign(n,0.0);current_z_.assign(n,0.0);
    electric_x_.assign(n,0.0);electric_y_.assign(n,0.0);electric_z_.assign(n,0.0);
    magnetic_x_.assign(n,0.0);magnetic_y_.assign(n,0.0);magnetic_z_.assign(n,0.0);
}

double StaggeredElectromagneticPic3D::epsilon() const noexcept {return epsilon0*config_.relative_permittivity;}
double StaggeredElectromagneticPic3D::mu() const noexcept {return mu0*config_.relative_permeability;}
double StaggeredElectromagneticPic3D::medium_light_speed() const noexcept {return 1.0/std::sqrt(epsilon()*mu());}
bool StaggeredElectromagneticPic3D::periodic_fields() const noexcept {return config_.field_boundary.mode==GridBoundaryMode3D::periodic;}

Vec3 StaggeredElectromagneticPic3D::place_particle(Vec3 position_m) const {
    if(config_.periodic_particles||periodic_fields()){
        position_m.x=wrap_periodic_value(position_m.x,config_.length_x_m);
        position_m.y=wrap_periodic_value(position_m.y,config_.length_y_m);
        position_m.z=wrap_periodic_value(position_m.z,config_.length_z_m);
    }else{
        position_m.x=std::clamp(position_m.x,0.0,config_.length_x_m);
        position_m.y=std::clamp(position_m.y,0.0,config_.length_y_m);
        position_m.z=std::clamp(position_m.z,0.0,config_.length_z_m);
    }
    return position_m;
}

void StaggeredElectromagneticPic3D::set_particles(std::vector<PicParticle3D> particles){
    for(auto& particle:particles){if(!(particle.mass_kg>0.0)||!std::isfinite(particle.charge_c)||!std::isfinite(particle.weight))throw std::invalid_argument("invalid staggered 3-D PIC particle");particle.position_m=place_particle(particle.position_m);}particles_=std::move(particles);
}

void StaggeredElectromagneticPic3D::set_fields(std::span<const double> electric_x_v_per_m,
                                               std::span<const double> electric_y_v_per_m,
                                               std::span<const double> electric_z_v_per_m,
                                               std::span<const double> magnetic_x_t,
                                               std::span<const double> magnetic_y_t,
                                               std::span<const double> magnetic_z_t){
    const std::size_t n=config_.nx*config_.ny*config_.nz;
    if(electric_x_v_per_m.size()!=n||electric_y_v_per_m.size()!=n||electric_z_v_per_m.size()!=n||magnetic_x_t.size()!=n||magnetic_y_t.size()!=n||magnetic_z_t.size()!=n)throw std::invalid_argument("staggered 3-D field arrays must match grid size");
    electric_x_.assign(electric_x_v_per_m.begin(),electric_x_v_per_m.end());electric_y_.assign(electric_y_v_per_m.begin(),electric_y_v_per_m.end());electric_z_.assign(electric_z_v_per_m.begin(),electric_z_v_per_m.end());
    magnetic_x_.assign(magnetic_x_t.begin(),magnetic_x_t.end());magnetic_y_.assign(magnetic_y_t.begin(),magnetic_y_t.end());magnetic_z_.assign(magnetic_z_t.begin(),magnetic_z_t.end());apply_field_boundary();
}

void StaggeredElectromagneticPic3D::initialize_z_polarized_mode(double electric_z_amplitude_v_per_m,std::size_t mode_x,std::size_t mode_y,std::size_t mode_z){
    std::fill(electric_x_.begin(),electric_x_.end(),0.0);std::fill(electric_y_.begin(),electric_y_.end(),0.0);std::fill(electric_z_.begin(),electric_z_.end(),0.0);
    std::fill(magnetic_x_.begin(),magnetic_x_.end(),0.0);std::fill(magnetic_y_.begin(),magnetic_y_.end(),0.0);std::fill(magnetic_z_.begin(),magnetic_z_.end(),0.0);
    if(mode_x==0U&&mode_y==0U&&mode_z==0U)throw std::invalid_argument("staggered 3-D EM-PIC mode cannot be DC");
    const double kx=2.0*std::numbers::pi*static_cast<double>(mode_x)/config_.length_x_m;
    const double ky=2.0*std::numbers::pi*static_cast<double>(mode_y)/config_.length_y_m;
    const double kz=2.0*std::numbers::pi*static_cast<double>(mode_z)/config_.length_z_m;
    const double k=std::sqrt(kx*kx+ky*ky+kz*kz);
    for(std::size_t iz=0;iz<config_.nz;++iz){
        const double z_ez=(static_cast<double>(iz)+0.5)*dz();
        const double z_bx=(static_cast<double>(iz)+0.5)*dz();
        const double z_by=(static_cast<double>(iz)+0.5)*dz();
        for(std::size_t iy=0;iy<config_.ny;++iy){
            const double y_ez=static_cast<double>(iy)*dy();
            const double y_bx=(static_cast<double>(iy)+0.5)*dy();
            const double y_by=static_cast<double>(iy)*dy();
            for(std::size_t ix=0;ix<config_.nx;++ix){
                const double x_ez=static_cast<double>(ix)*dx();
                const double x_bx=static_cast<double>(ix)*dx();
                const double x_by=(static_cast<double>(ix)+0.5)*dx();
                const std::size_t idx=flat_index_3d(ix,iy,iz,config_.nx,config_.ny);
                const double phase_e=kx*x_ez+ky*y_ez+kz*z_ez;
                electric_z_[idx]=electric_z_amplitude_v_per_m*std::sin(phase_e);
                magnetic_x_[idx]=(ky/k)*(electric_z_amplitude_v_per_m/medium_light_speed())*std::sin(kx*x_bx+ky*y_bx+kz*z_bx);
                magnetic_y_[idx]=-(kx/k)*(electric_z_amplitude_v_per_m/medium_light_speed())*std::sin(kx*x_by+ky*y_by+kz*z_by);
            }
        }
    }
    apply_field_boundary();
}

double StaggeredElectromagneticPic3D::sample(std::span<const double> values,Vec3 position_m,double offset_x,double offset_y,double offset_z) const {
    const auto wrap=[&](double value,double length){return periodic_fields()?wrap_periodic_value(value,length):std::clamp(value,0.0,std::nextafter(length,0.0));};
    const double x=wrap(position_m.x-offset_x*dx(),config_.length_x_m);
    const double y=wrap(position_m.y-offset_y*dy(),config_.length_y_m);
    const double z=wrap(position_m.z-offset_z*dz(),config_.length_z_m);
    const double ux=x/dx(),uy=y/dy(),uz=z/dz();
    const auto ix0=static_cast<std::size_t>(std::floor(ux))%config_.nx;const auto iy0=static_cast<std::size_t>(std::floor(uy))%config_.ny;const auto iz0=static_cast<std::size_t>(std::floor(uz))%config_.nz;
    const auto ix1=periodic_fields()?(ix0+1U)%config_.nx:std::min(ix0+1U,config_.nx-1U);const auto iy1=periodic_fields()?(iy0+1U)%config_.ny:std::min(iy0+1U,config_.ny-1U);const auto iz1=periodic_fields()?(iz0+1U)%config_.nz:std::min(iz0+1U,config_.nz-1U);
    const double fx=ux-std::floor(ux),fy=uy-std::floor(uy),fz=uz-std::floor(uz);
    double value=0.0;for(std::size_t dz_i=0;dz_i<2U;++dz_i){const std::size_t iz=dz_i?iz1:iz0;const double wz=dz_i?fz:(1.0-fz);for(std::size_t dy_i=0;dy_i<2U;++dy_i){const std::size_t iy=dy_i?iy1:iy0;const double wy=dy_i?fy:(1.0-fy);for(std::size_t dx_i=0;dx_i<2U;++dx_i){const std::size_t ix=dx_i?ix1:ix0;const double wx=dx_i?fx:(1.0-fx);value+=values[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)]*wx*wy*wz;}}}return value;
}

ElectromagneticField StaggeredElectromagneticPic3D::gather_field(Vec3 position_m) const {
    return {{sample(electric_x_,position_m,0.5,0.0,0.0),sample(electric_y_,position_m,0.0,0.5,0.0),sample(electric_z_,position_m,0.0,0.0,0.5)},
            {sample(magnetic_x_,position_m,0.0,0.5,0.5),sample(magnetic_y_,position_m,0.5,0.0,0.5),sample(magnetic_z_,position_m,0.5,0.5,0.0)}};
}

double StaggeredElectromagneticPic3D::diff_x_forward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t r=(ix+1U<config_.nx)?ix+1U:(periodic_fields()?0U:ix);return (values[flat_index_3d(r,iy,iz,config_.nx,config_.ny)]-values[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)])/dx();}
double StaggeredElectromagneticPic3D::diff_y_forward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t u=(iy+1U<config_.ny)?iy+1U:(periodic_fields()?0U:iy);return (values[flat_index_3d(ix,u,iz,config_.nx,config_.ny)]-values[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)])/dy();}
double StaggeredElectromagneticPic3D::diff_z_forward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t f=(iz+1U<config_.nz)?iz+1U:(periodic_fields()?0U:iz);return (values[flat_index_3d(ix,iy,f,config_.nx,config_.ny)]-values[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)])/dz();}
double StaggeredElectromagneticPic3D::diff_x_backward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t l=(ix>0U)?ix-1U:(periodic_fields()?config_.nx-1U:ix);return (values[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)]-values[flat_index_3d(l,iy,iz,config_.nx,config_.ny)])/dx();}
double StaggeredElectromagneticPic3D::diff_y_backward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t d=(iy>0U)?iy-1U:(periodic_fields()?config_.ny-1U:iy);return (values[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)]-values[flat_index_3d(ix,d,iz,config_.nx,config_.ny)])/dy();}
double StaggeredElectromagneticPic3D::diff_z_backward(std::span<const double> values,std::size_t ix,std::size_t iy,std::size_t iz) const {const std::size_t b=(iz>0U)?iz-1U:(periodic_fields()?config_.nz-1U:iz);return (values[flat_index_3d(ix,iy,iz,config_.nx,config_.ny)]-values[flat_index_3d(ix,iy,b,config_.nx,config_.ny)])/dz();}

void StaggeredElectromagneticPic3D::advance_magnetic(double dt_s){
    auto bx=magnetic_x_,by=magnetic_y_,bz=magnetic_z_;
    for(std::size_t iz=0;iz<config_.nz;++iz){for(std::size_t iy=0;iy<config_.ny;++iy){for(std::size_t ix=0;ix<config_.nx;++ix){
        const std::size_t idx=flat_index_3d(ix,iy,iz,config_.nx,config_.ny);
        bx[idx]-=dt_s*(diff_y_forward(electric_z_,ix,iy,iz)-diff_z_forward(electric_y_,ix,iy,iz));
        by[idx]-=dt_s*(diff_z_forward(electric_x_,ix,iy,iz)-diff_x_forward(electric_z_,ix,iy,iz));
        bz[idx]-=dt_s*(diff_x_forward(electric_y_,ix,iy,iz)-diff_y_forward(electric_x_,ix,iy,iz));
    }}}
    magnetic_x_.swap(bx);magnetic_y_.swap(by);magnetic_z_.swap(bz);apply_field_boundary();
}

void StaggeredElectromagneticPic3D::update_electric(double dt_s){
    auto ex=electric_x_,ey=electric_y_,ez=electric_z_;const double eps=epsilon();const double permeability=mu();
    for(std::size_t iz=0;iz<config_.nz;++iz){for(std::size_t iy=0;iy<config_.ny;++iy){for(std::size_t ix=0;ix<config_.nx;++ix){
        const std::size_t idx=flat_index_3d(ix,iy,iz,config_.nx,config_.ny);
        const double curl_h_x=(diff_y_backward(magnetic_z_,ix,iy,iz)-diff_z_backward(magnetic_y_,ix,iy,iz))/permeability;
        const double curl_h_y=(diff_z_backward(magnetic_x_,ix,iy,iz)-diff_x_backward(magnetic_z_,ix,iy,iz))/permeability;
        const double curl_h_z=(diff_x_backward(magnetic_y_,ix,iy,iz)-diff_y_backward(magnetic_x_,ix,iy,iz))/permeability;
        ex[idx]+=dt_s*(curl_h_x-current_x_[idx])/eps;ey[idx]+=dt_s*(curl_h_y-current_y_[idx])/eps;ez[idx]+=dt_s*(curl_h_z-current_z_[idx])/eps;
    }}}
    electric_x_.swap(ex);electric_y_.swap(ey);electric_z_.swap(ez);apply_field_boundary();
}

void StaggeredElectromagneticPic3D::apply_field_boundary(){
    apply_electromagnetic_field_boundary_3d(electric_x_,electric_y_,electric_z_,magnetic_x_,magnetic_y_,magnetic_z_,config_.nx,config_.ny,config_.nz,config_.field_boundary);
}

void StaggeredElectromagneticPic3D::apply_particle_boundary(){
    if(config_.periodic_particles||periodic_fields()){for(auto& particle:particles_)particle.position_m=place_particle(particle.position_m);return;}
    AxisAlignedParticleBox box;box.minimum_m={0.0,0.0,0.0};box.maximum_m={config_.length_x_m,config_.length_y_m,config_.length_z_m};
    std::vector<PicParticle3D> survivors;survivors.reserve(particles_.size());
    for(auto& particle:particles_){
        ChargedParticle charged;charged.position_m=particle.position_m;charged.velocity_m_per_s=particle.velocity_m_per_s;charged.charge_c=particle.charge_c;charged.mass_kg=particle.mass_kg;charged.weight=particle.weight;
        const auto interaction=apply_particle_box_boundary(charged,box,config_.particle_boundary,config_.enable_secondary_yield_report?&config_.secondary_model:nullptr);
        if(interaction.impacted){++absorbed_particles_;secondary_macro_weight_+=interaction.secondary_macro_weight;}
        particle.position_m=charged.position_m;particle.velocity_m_per_s=charged.velocity_m_per_s;
        if(interaction.alive||!config_.remove_absorbed_particles)survivors.push_back(particle);
    }
    particles_.swap(survivors);
}

void StaggeredElectromagneticPic3D::maybe_sort_particles(){
    ++step_counter_;
    if(!config_.sort_particles_by_cell||config_.particle_sort_interval==0U||step_counter_%config_.particle_sort_interval!=0U)return;
    auto sorted=sort_particles_by_cell_3d(particles_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m);
    occupied_particle_cells_=sorted.occupied_cells;
    particles_=std::move(sorted.sorted_particles);
    ++particle_sort_passes_;
}

void StaggeredElectromagneticPic3D::step(std::size_t steps,const NeutralCollisionModel* collision_model){
    for(std::size_t step_index=0;step_index<steps;++step_index){
        const auto particles_old=particles_;
        advance_magnetic(0.5*config_.dt_s);
        for(auto& particle:particles_){
            ChargedParticle charged;charged.position_m=particle.position_m;charged.velocity_m_per_s=particle.velocity_m_per_s;charged.charge_c=particle.charge_c;charged.mass_kg=particle.mass_kg;charged.weight=particle.weight;
            vay_push(charged,gather_field(particle.position_m),config_.dt_s);
            particle.position_m=charged.position_m;particle.velocity_m_per_s=charged.velocity_m_per_s;
        }
        apply_particle_boundary();
        if(particles_old.size()==particles_.size()){
            const auto current=(config_.current_deposition==CurrentDeposition3DMode::local_finite_volume)
                ?deposit_charge_conserving_current_3d_local(particles_old,particles_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m,config_.dt_s)
                :deposit_charge_conserving_current_3d(particles_old,particles_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m,config_.dt_s);
            charge_density_=current.charge_density_new_c_per_m3;current_x_=current.current_x_a_per_m2;current_y_=current.current_y_a_per_m2;current_z_=current.current_z_a_per_m2;last_continuity_residual_=current.continuity_linf_residual;
        }else{
            charge_density_=deposit_cic_charge_density_3d(particles_,config_.nx,config_.ny,config_.nz,config_.length_x_m,config_.length_y_m,config_.length_z_m);
            std::fill(current_x_.begin(),current_x_.end(),0.0);std::fill(current_y_.begin(),current_y_.end(),0.0);std::fill(current_z_.begin(),current_z_.end(),0.0);last_continuity_residual_=0.0;
        }
        update_electric(config_.dt_s);
        advance_magnetic(0.5*config_.dt_s);
        if(collision_model){
            std::vector<ChargedParticle> charged_particles;charged_particles.reserve(particles_.size());
            for(const auto& particle:particles_)charged_particles.push_back({particle.position_m,particle.velocity_m_per_s,particle.charge_c,particle.mass_kg,particle.weight});
            const auto stats=apply_monte_carlo_collisions(charged_particles,config_.dt_s,*collision_model);accumulated_collisions_.elastic_events+=stats.elastic_events;accumulated_collisions_.ionization_events+=stats.ionization_events;accumulated_collisions_.energy_loss_j+=stats.energy_loss_j;
            particles_.clear();particles_.reserve(charged_particles.size());for(const auto& charged:charged_particles){PicParticle3D particle;particle.position_m=place_particle(charged.position_m);particle.velocity_m_per_s=charged.velocity_m_per_s;particle.charge_c=charged.charge_c;particle.mass_kg=charged.mass_kg;particle.weight=charged.weight;particles_.push_back(particle);}apply_particle_boundary();
        }
        time_s_+=config_.dt_s;
        maybe_sort_particles();
    }
}

StaggeredElectromagneticPic3DDiagnostics StaggeredElectromagneticPic3D::diagnostics() const {
    StaggeredElectromagneticPic3DDiagnostics result;const double volume=dx()*dy()*dz();const double eps=epsilon();const double permeability=mu();
    for(std::size_t idx=0;idx<electric_x_.size();++idx){const double e2=electric_x_[idx]*electric_x_[idx]+electric_y_[idx]*electric_y_[idx]+electric_z_[idx]*electric_z_[idx];const double b2=magnetic_x_[idx]*magnetic_x_[idx]+magnetic_y_[idx]*magnetic_y_[idx]+magnetic_z_[idx]*magnetic_z_[idx];result.field_energy_j+=0.5*eps*e2*volume+0.5*b2*volume/permeability;}
    for(const auto& particle:particles_){ChargedParticle charged;charged.position_m=particle.position_m;charged.velocity_m_per_s=particle.velocity_m_per_s;charged.charge_c=particle.charge_c;charged.mass_kg=particle.mass_kg;charged.weight=particle.weight;result.particle_kinetic_energy_j+=particle.weight*kinetic_energy_j(charged);}
    result.total_energy_j=result.field_energy_j+result.particle_kinetic_energy_j;result.charge_continuity_linf_residual=last_continuity_residual_;result.particle_count=particles_.size();result.absorbed_particles=absorbed_particles_;result.secondary_macro_weight=secondary_macro_weight_;result.collision_statistics=accumulated_collisions_;result.particle_sort_passes=particle_sort_passes_;result.occupied_particle_cells=occupied_particle_cells_;return result;
}

} // namespace cfd::particle
