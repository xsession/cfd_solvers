#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace cfd::particle {
namespace {
constexpr double epsilon0=8.8541878128e-12;
Vec3 add(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Vec3 scale(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
Vec3 cross(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double norm2(Vec3 a){return a.x*a.x+a.y*a.y+a.z*a.z;}
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

} // namespace cfd::particle
