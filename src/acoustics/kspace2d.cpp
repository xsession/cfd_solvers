#include "cfd/acoustics/kspace2d.hpp"
#include "cfd/spectral/fft.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace cfd::acoustics {
namespace {
[[nodiscard]] double sinc(double x){return std::abs(x)<1.0e-12?1.0:std::sin(x)/x;}
[[nodiscard]] std::size_t index(std::size_t i,std::size_t j,std::size_t nx){return j*nx+i;}
}

KSpaceAcoustic2D::KSpaceAcoustic2D(KSpace2DConfig config):config_(config){
    if(!cfd::spectral::is_power_of_two(config_.nx)||!cfd::spectral::is_power_of_two(config_.ny)||
       !(config_.dx>0.0)||!(config_.dy>0.0)||!(config_.cfl>0.0)||
       !(config_.pml_target_reflection>0.0&&config_.pml_target_reflection<1.0)||!(config_.pml_order>0.0)||
       config_.pml_cells*2U>=std::min(config_.nx,config_.ny))
        throw std::invalid_argument("invalid 2-D k-space acoustic configuration");
    const auto n=cells();
    c0_.assign(n,1500.0);rho0_.assign(n,1000.0);bon_a_.assign(n,0.0);
    ux_.assign(n,0.0);uy_.assign(n,0.0);rho_x_.assign(n,0.0);rho_y_.assign(n,0.0);rho_.assign(n,0.0);pressure_.assign(n,0.0);
    damp_x_.resize(config_.nx);damp_y_.resize(config_.ny);
    kx_=cfd::spectral::wavenumbers(config_.nx,config_.dx);
    ky_=cfd::spectral::wavenumbers(config_.ny,config_.dy);
    kappa_.resize(n);
    rebuild_time_step_and_kappa();rebuild_pml();
}

void KSpaceAcoustic2D::set_uniform_medium(double c,double rho){
    if(!(c>0.0)||!(rho>0.0)||!std::isfinite(c)||!std::isfinite(rho)) throw std::invalid_argument("invalid acoustic medium");
    std::fill(c0_.begin(),c0_.end(),c);std::fill(rho0_.begin(),rho0_.end(),rho);rebuild_time_step_and_kappa();rebuild_pml();
}
void KSpaceAcoustic2D::set_medium(std::span<const double> c,std::span<const double> rho){
    if(c.size()!=cells()||rho.size()!=cells()) throw std::invalid_argument("acoustic medium size mismatch");
    for(std::size_t q=0;q<cells();++q) if(!(c[q]>0.0)||!(rho[q]>0.0)||!std::isfinite(c[q])||!std::isfinite(rho[q])) throw std::invalid_argument("invalid heterogeneous acoustic medium");
    c0_.assign(c.begin(),c.end());rho0_.assign(rho.begin(),rho.end());rebuild_time_step_and_kappa();rebuild_pml();
}
void KSpaceAcoustic2D::set_nonlinearity(std::span<const double> b){if(b.size()!=cells())throw std::invalid_argument("BonA size mismatch");for(double value:b)if(!std::isfinite(value))throw std::invalid_argument("invalid BonA field");bon_a_.assign(b.begin(),b.end());}
void KSpaceAcoustic2D::set_uniform_nonlinearity(double b){if(!std::isfinite(b))throw std::invalid_argument("invalid BonA");std::fill(bon_a_.begin(),bon_a_.end(),b);}
void KSpaceAcoustic2D::set_power_law_absorption(PowerLawAbsorption2D a){if(!(a.alpha0_np_per_m>=0.0)||!(a.exponent>0.0)||!(a.reference_frequency_hz>0.0))throw std::invalid_argument("invalid power-law absorption");absorption_=a;absorption_enabled_=a.alpha0_np_per_m>0.0;}
void KSpaceAcoustic2D::clear_absorption() noexcept {absorption_enabled_=false;}

void KSpaceAcoustic2D::initialize_pressure(std::span<const double> p){
    if(p.size()!=cells())throw std::invalid_argument("initial pressure size mismatch");
    pressure_.assign(p.begin(),p.end());
    for(std::size_t q=0;q<cells();++q){rho_[q]=pressure_[q]/(c0_[q]*c0_[q]);rho_x_[q]=0.5*rho_[q];rho_y_[q]=0.5*rho_[q];}
    std::fill(ux_.begin(),ux_.end(),0.0);std::fill(uy_.begin(),uy_.end(),0.0);time_=0.0;record_sensors();
}
void KSpaceAcoustic2D::initialize_gaussian_pressure(double cx,double cy,double sigma,double amplitude){
    if(!(sigma>0.0)||!std::isfinite(amplitude))throw std::invalid_argument("invalid Gaussian pressure");
    std::vector<double> p(cells());const double x0=cx*static_cast<double>(config_.nx)*config_.dx,y0=cy*static_cast<double>(config_.ny)*config_.dy;
    for(std::size_t j=0;j<config_.ny;++j)for(std::size_t i=0;i<config_.nx;++i){const double x=(static_cast<double>(i)+0.5)*config_.dx,y=(static_cast<double>(j)+0.5)*config_.dy;const double r2=(x-x0)*(x-x0)+(y-y0)*(y-y0);p[index(i,j,config_.nx)]=amplitude*std::exp(-0.5*r2/(sigma*sigma));}
    initialize_pressure(p);
}
void KSpaceAcoustic2D::add_pressure_source(std::size_t cell,std::function<double(double)> value){if(cell>=cells()||!value)throw std::invalid_argument("invalid acoustic source");sources_.push_back({cell,std::move(value)});}
void KSpaceAcoustic2D::add_plane_source_x(std::size_t i,std::function<double(double)> value){if(i>=config_.nx||!value)throw std::invalid_argument("invalid acoustic plane source");for(std::size_t j=0;j<config_.ny;++j)sources_.push_back({index(i,j,config_.nx),value});}
void KSpaceAcoustic2D::add_delayed_source_array(std::span<const std::size_t> cells_in,std::span<const double> delays,std::function<double(double)> waveform){if(cells_in.empty()||cells_in.size()!=delays.size()||!waveform)throw std::invalid_argument("invalid delayed acoustic source array");for(std::size_t n=0;n<cells_in.size();++n){if(cells_in[n]>=cells()||delays[n]<0.0||!std::isfinite(delays[n]))throw std::invalid_argument("invalid delayed acoustic source entry");const double delay=delays[n];sources_.push_back({cells_in[n],[waveform,delay](double t){return t>=delay?waveform(t-delay):0.0;}});}}
std::size_t KSpaceAcoustic2D::add_sensor(std::size_t cell){if(cell>=cells())throw std::out_of_range("sensor cell out of range");sensors_.push_back({cell,{}});return sensors_.size()-1U;}
std::vector<std::size_t> KSpaceAcoustic2D::add_sensor_array(std::span<const std::size_t> cells_in){std::vector<std::size_t> ids;ids.reserve(cells_in.size());for(const auto cell:cells_in)ids.push_back(add_sensor(cell));return ids;}
void KSpaceAcoustic2D::clear_sensors(){sensors_.clear();}

std::vector<double> KSpaceAcoustic2D::time_reversal_reconstruct(
    std::span<const std::size_t> sensor_cells,
    std::span<const std::vector<double>> sensor_pressure,
    bool positivity) {
    if(sensor_cells.empty()||sensor_cells.size()!=sensor_pressure.size())
        throw std::invalid_argument("invalid time-reversal sensor data");
    const std::size_t samples=sensor_pressure.front().size();
    if(samples<2U)throw std::invalid_argument("time-reversal traces are too short");
    for(std::size_t sidx=0;sidx<sensor_cells.size();++sidx){
        if(sensor_cells[sidx]>=cells()||sensor_pressure[sidx].size()!=samples)
            throw std::invalid_argument("time-reversal trace size/cell mismatch");
        for(double v:sensor_pressure[sidx])if(!std::isfinite(v))throw std::invalid_argument("non-finite time-reversal pressure");
    }
    std::fill(ux_.begin(),ux_.end(),0.0);std::fill(uy_.begin(),uy_.end(),0.0);
    std::fill(rho_x_.begin(),rho_x_.end(),0.0);std::fill(rho_y_.begin(),rho_y_.end(),0.0);
    std::fill(rho_.begin(),rho_.end(),0.0);std::fill(pressure_.begin(),pressure_.end(),0.0);time_=0.0;
    std::vector<double> gx,gy,dux_dx,duy_dy,unused;
    for(std::size_t n=0;n<samples;++n){
        spectral_gradient(pressure_,gx,gy);
        for(std::size_t j=0;j<config_.ny;++j)for(std::size_t i=0;i<config_.nx;++i){const auto q=index(i,j,config_.nx);ux_[q]=(ux_[q]-dt_*gx[q]/rho0_[q])*damp_x_[i];uy_[q]=(uy_[q]-dt_*gy[q]/rho0_[q])*damp_y_[j];}
        spectral_gradient(ux_,dux_dx,unused);spectral_gradient(uy_,unused,duy_dy);
        for(std::size_t j=0;j<config_.ny;++j)for(std::size_t i=0;i<config_.nx;++i){const auto q=index(i,j,config_.nx);rho_x_[q]=(rho_x_[q]-dt_*rho0_[q]*dux_dx[q])*damp_x_[i];rho_y_[q]=(rho_y_[q]-dt_*rho0_[q]*duy_dy[q])*damp_y_[j];}
        update_pressure_from_density();
        // Time reversal is run without forward sources or attenuation. Recorded
        // data are imposed as a Dirichlet boundary in reversed time order.
        const std::size_t src=samples-1U-n;
        for(std::size_t k=0;k<sensor_cells.size();++k){const auto q=sensor_cells[k];const double p=sensor_pressure[k][src];pressure_[q]=p;const double dr=p/(c0_[q]*c0_[q]);rho_[q]=dr;rho_x_[q]=0.5*dr;rho_y_[q]=0.5*dr;}
        time_+=dt_;
    }
    if(positivity)for(double& p:pressure_)p=std::max(0.0,p);
    return pressure_;
}

void KSpaceAcoustic2D::rebuild_time_step_and_kappa(){
    const double cmax=*std::max_element(c0_.begin(),c0_.end());
    const double cref=config_.reference_sound_speed>0.0?config_.reference_sound_speed:cmax;
    dt_=config_.dt>0.0?config_.dt:config_.cfl*std::min(config_.dx,config_.dy)/cmax;
    if(!(dt_>0.0)||dt_*cmax/std::min(config_.dx,config_.dy)>1.0)throw std::invalid_argument("unstable/invalid acoustic time step");
    for(std::size_t j=0;j<config_.ny;++j)for(std::size_t i=0;i<config_.nx;++i){const double kmag=std::hypot(kx_[i],ky_[j]);kappa_[index(i,j,config_.nx)]=sinc(0.5*cref*kmag*dt_);}
}
void KSpaceAcoustic2D::rebuild_pml(){
    std::fill(damp_x_.begin(),damp_x_.end(),1.0);std::fill(damp_y_.begin(),damp_y_.end(),1.0);
    if(config_.pml_cells==0U)return;
    const double cmax=*std::max_element(c0_.begin(),c0_.end());
    const double thickness=static_cast<double>(config_.pml_cells)*std::min(config_.dx,config_.dy);
    const double sigma_max=-(config_.pml_order+1.0)*cmax*std::log(config_.pml_target_reflection)/(2.0*thickness);
    auto fill=[&](std::vector<double>& d,std::size_t n){for(std::size_t i=0;i<n;++i){const std::size_t depth=std::min(i,n-1U-i);if(depth<config_.pml_cells){const double s=static_cast<double>(config_.pml_cells-depth)/static_cast<double>(config_.pml_cells);d[i]=std::exp(-sigma_max*std::pow(s,config_.pml_order)*dt_);}}};
    fill(damp_x_,config_.nx);fill(damp_y_,config_.ny);
}

void KSpaceAcoustic2D::spectral_gradient(std::span<const double> input,std::vector<double>& gx,std::vector<double>& gy) const{
    std::vector<cfd::spectral::Complex> f(cells());for(std::size_t q=0;q<cells();++q)f[q]=input[q];
    cfd::spectral::fft2_inplace(f,config_.nx,config_.ny,false);
    std::vector<cfd::spectral::Complex> fx=f,fy=f;const cfd::spectral::Complex imag{0.0,1.0};
    for(std::size_t j=0;j<config_.ny;++j)for(std::size_t i=0;i<config_.nx;++i){const auto q=index(i,j,config_.nx);const double kdx=(i==config_.nx/2U)?0.0:kx_[i];const double kdy=(j==config_.ny/2U)?0.0:ky_[j];fx[q]*=imag*kdx*kappa_[q];fy[q]*=imag*kdy*kappa_[q];}
    cfd::spectral::fft2_inplace(fx,config_.nx,config_.ny,true);cfd::spectral::fft2_inplace(fy,config_.nx,config_.ny,true);
    gx.resize(cells());gy.resize(cells());for(std::size_t q=0;q<cells();++q){gx[q]=fx[q].real();gy[q]=fy[q].real();}
}

void KSpaceAcoustic2D::update_pressure_from_density(){
    for(std::size_t q=0;q<cells();++q){rho_[q]=rho_x_[q]+rho_y_[q];const double nonlinear=1.0+0.5*bon_a_[q]*rho_[q]/rho0_[q];pressure_[q]=c0_[q]*c0_[q]*rho_[q]*nonlinear;}
}

void KSpaceAcoustic2D::apply_absorption(){
    if(!absorption_enabled_)return;
    const double cref=config_.reference_sound_speed>0.0?config_.reference_sound_speed:*std::max_element(c0_.begin(),c0_.end());
    pressure_=power_law_attenuate_periodic_2d(pressure_,config_.nx,config_.ny,config_.dx,config_.dy,cref,dt_,absorption_);
    for(std::size_t q=0;q<cells();++q){rho_[q]=pressure_[q]/(c0_[q]*c0_[q]);rho_x_[q]=0.5*rho_[q];rho_y_[q]=0.5*rho_[q];}
}

void KSpaceAcoustic2D::step(std::size_t count){
    std::vector<double> gx,gy,dux_dx,duy_dy;
    for(std::size_t n=0;n<count;++n){
        spectral_gradient(pressure_,gx,gy);
        for(std::size_t j=0;j<config_.ny;++j)for(std::size_t i=0;i<config_.nx;++i){const auto q=index(i,j,config_.nx);ux_[q]=(ux_[q]-dt_*gx[q]/rho0_[q])*damp_x_[i];uy_[q]=(uy_[q]-dt_*gy[q]/rho0_[q])*damp_y_[j];}
        std::vector<double> unused;
        spectral_gradient(ux_,dux_dx,unused);
        spectral_gradient(uy_,unused,duy_dy);
        for(std::size_t j=0;j<config_.ny;++j)for(std::size_t i=0;i<config_.nx;++i){const auto q=index(i,j,config_.nx);rho_x_[q]=(rho_x_[q]-dt_*rho0_[q]*dux_dx[q])*damp_x_[i];rho_y_[q]=(rho_y_[q]-dt_*rho0_[q]*duy_dy[q])*damp_y_[j];}
        update_pressure_from_density();apply_absorption();
        const double next_time=time_+dt_;
        for(const auto& source:sources_){
            const double dp=source.value(next_time);
            pressure_[source.cell]+=dp;
            const double dr=dp/(c0_[source.cell]*c0_[source.cell]);
            rho_[source.cell]+=dr;rho_x_[source.cell]+=0.5*dr;rho_y_[source.cell]+=0.5*dr;
        }
        time_=next_time;record_sensors();
    }
}

void KSpaceAcoustic2D::record_sensors(){for(auto& s:sensors_)s.pressure.push_back(pressure_[s.cell]);}

double KSpaceAcoustic2D::total_energy() const{
    double e=0.0;const double dv=config_.dx*config_.dy;
    for(std::size_t q=0;q<cells();++q)e+=0.5*(rho0_[q]*(ux_[q]*ux_[q]+uy_[q]*uy_[q])+pressure_[q]*pressure_[q]/(rho0_[q]*c0_[q]*c0_[q]))*dv;
    return e;
}

std::vector<double> fractional_laplacian_periodic_2d(std::span<const double> field,std::size_t nx,std::size_t ny,double dx,double dy,double order){
    if(field.size()!=nx*ny||!cfd::spectral::is_power_of_two(nx)||!cfd::spectral::is_power_of_two(ny)||!(dx>0.0)||!(dy>0.0)||!(order>=0.0))throw std::invalid_argument("invalid 2-D fractional Laplacian request");
    const auto kx=cfd::spectral::wavenumbers(nx,dx),ky=cfd::spectral::wavenumbers(ny,dy);
    std::vector<cfd::spectral::Complex> f(field.size());for(std::size_t q=0;q<field.size();++q)f[q]=field[q];cfd::spectral::fft2_inplace(f,nx,ny,false);
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i){const double k2=kx[i]*kx[i]+ky[j]*ky[j];const double factor=k2==0.0?(order==0.0?1.0:0.0):std::pow(k2,order);f[index(i,j,nx)]*=factor;}
    cfd::spectral::fft2_inplace(f,nx,ny,true);std::vector<double> out(field.size());for(std::size_t q=0;q<field.size();++q)out[q]=f[q].real();return out;
}

std::vector<double> power_law_attenuate_periodic_2d(std::span<const double> field,std::size_t nx,std::size_t ny,double dx,double dy,double c_ref,double dt,PowerLawAbsorption2D a){
    if(field.size()!=nx*ny||!cfd::spectral::is_power_of_two(nx)||!cfd::spectral::is_power_of_two(ny)||!(c_ref>0.0)||!(dt>=0.0)||!(a.alpha0_np_per_m>=0.0)||!(a.exponent>0.0)||!(a.reference_frequency_hz>0.0))throw std::invalid_argument("invalid power-law attenuation request");
    const auto kx=cfd::spectral::wavenumbers(nx,dx),ky=cfd::spectral::wavenumbers(ny,dy);
    std::vector<cfd::spectral::Complex> f(field.size());for(std::size_t q=0;q<field.size();++q)f[q]=field[q];cfd::spectral::fft2_inplace(f,nx,ny,false);
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i){const double kmag=std::hypot(kx[i],ky[j]);const double freq=c_ref*kmag/(2.0*std::numbers::pi);const double alpha=a.alpha0_np_per_m*std::pow(freq/a.reference_frequency_hz,a.exponent);f[index(i,j,nx)]*=std::exp(-alpha*c_ref*dt);}
    cfd::spectral::fft2_inplace(f,nx,ny,true);std::vector<double> out(field.size());for(std::size_t q=0;q<field.size();++q)out[q]=f[q].real();return out;
}

} // namespace cfd::acoustics
