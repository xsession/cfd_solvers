#include "cfd/solvers/fvm/compressible1d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::fvm {
namespace {

using Vector3=std::array<double,3>;
using Matrix3=std::array<Vector3,3>;

double minmod(double a,double b){return a*b<=0.0?0.0:(std::abs(a)<std::abs(b)?a:b);}

EulerState1D conservative(EulerPrimitive1D q,const IdealGasEquationOfState& eos){
    if(!(q.density>0.0)||!(q.pressure>0.0)||!std::isfinite(q.density)||!std::isfinite(q.velocity)||!std::isfinite(q.pressure))
        throw std::invalid_argument("non-positive or non-finite compressible primitive state");
    return {q.density,q.density*q.velocity,q.pressure/(eos.gamma-1.0)+0.5*q.density*q.velocity*q.velocity};
}

double raw_pressure(const EulerState1D& u,const IdealGasEquationOfState& eos) noexcept {
    if(!(u.density>0.0)||!std::isfinite(u.density)||!std::isfinite(u.momentum)||!std::isfinite(u.energy))
        return -std::numeric_limits<double>::infinity();
    return (eos.gamma-1.0)*(u.energy-0.5*u.momentum*u.momentum/u.density);
}

EulerState1D flux(const EulerState1D& u,const IdealGasEquationOfState& eos){
    const double p=eos.pressure(u.density,u.momentum,u.energy);
    const double v=u.momentum/u.density;
    return {u.momentum,u.momentum*v+p,(u.energy+p)*v};
}

EulerState1D add(EulerState1D a,EulerState1D b,double scale=1.0){return {a.density+scale*b.density,a.momentum+scale*b.momentum,a.energy+scale*b.energy};}
EulerState1D mul(EulerState1D a,double s){return {a.density*s,a.momentum*s,a.energy*s};}

Vector3 vector_of(const EulerState1D& u){return {u.density,u.momentum,u.energy};}
EulerState1D state_of(const Vector3& v){return {v[0],v[1],v[2]};}

Vector3 multiply(const Matrix3& a,const Vector3& x){
    Vector3 y{};
    for(std::size_t r=0;r<3U;++r)for(std::size_t c=0;c<3U;++c)y[r]+=a[r][c]*x[c];
    return y;
}

Matrix3 inverse(const Matrix3& a){
    const double det=
        a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
       -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
       +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
    if(!(std::abs(det)>1.0e-14)||!std::isfinite(det))throw std::runtime_error("singular Euler characteristic basis");
    const double s=1.0/det;
    return {{{
        (a[1][1]*a[2][2]-a[1][2]*a[2][1])*s,
        (a[0][2]*a[2][1]-a[0][1]*a[2][2])*s,
        (a[0][1]*a[1][2]-a[0][2]*a[1][1])*s},
        {
        (a[1][2]*a[2][0]-a[1][0]*a[2][2])*s,
        (a[0][0]*a[2][2]-a[0][2]*a[2][0])*s,
        (a[0][2]*a[1][0]-a[0][0]*a[1][2])*s},
        {
        (a[1][0]*a[2][1]-a[1][1]*a[2][0])*s,
        (a[0][1]*a[2][0]-a[0][0]*a[2][1])*s,
        (a[0][0]*a[1][1]-a[0][1]*a[1][0])*s
    }}};
}

struct CharacteristicBasis {
    Matrix3 right{};
    Matrix3 left{};
};

CharacteristicBasis roe_basis(const EulerState1D& ul,const EulerState1D& ur,const IdealGasEquationOfState& eos){
    const double pl=eos.pressure(ul.density,ul.momentum,ul.energy);
    const double pr=eos.pressure(ur.density,ur.momentum,ur.energy);
    const double vl=ul.momentum/ul.density;
    const double vr=ur.momentum/ur.density;
    const double hl=(ul.energy+pl)/ul.density;
    const double hr=(ur.energy+pr)/ur.density;
    const double rl=std::sqrt(ul.density);
    const double rr=std::sqrt(ur.density);
    const double denom=rl+rr;
    const double v=(rl*vl+rr*vr)/denom;
    const double h=(rl*hl+rr*hr)/denom;
    const double c2=(eos.gamma-1.0)*(h-0.5*v*v);
    if(!(c2>0.0)||!std::isfinite(c2))throw std::runtime_error("invalid Roe-averaged sound speed");
    const double c=std::sqrt(c2);
    Matrix3 r{{
        {{1.0,1.0,1.0}},
        {{v-c,v,v+c}},
        {{h-v*c,0.5*v*v,h+v*c}}
    }};
    return {r,inverse(r)};
}

double weno5_left(const std::array<double,5>& q,double epsilon){
    const double p0=(2.0*q[0]-7.0*q[1]+11.0*q[2])/6.0;
    const double p1=(-q[1]+5.0*q[2]+2.0*q[3])/6.0;
    const double p2=(2.0*q[2]+5.0*q[3]-q[4])/6.0;
    const double b0=(13.0/12.0)*std::pow(q[0]-2.0*q[1]+q[2],2.0)+0.25*std::pow(q[0]-4.0*q[1]+3.0*q[2],2.0);
    const double b1=(13.0/12.0)*std::pow(q[1]-2.0*q[2]+q[3],2.0)+0.25*std::pow(q[1]-q[3],2.0);
    const double b2=(13.0/12.0)*std::pow(q[2]-2.0*q[3]+q[4],2.0)+0.25*std::pow(3.0*q[2]-4.0*q[3]+q[4],2.0);
    const double a0=0.1/std::pow(epsilon+b0,2.0);
    const double a1=0.6/std::pow(epsilon+b1,2.0);
    const double a2=0.3/std::pow(epsilon+b2,2.0);
    const double inv=1.0/(a0+a1+a2);
    return (a0*p0+a1*p1+a2*p2)*inv;
}

EulerState1D rusanov_flux(const EulerState1D& ul,const EulerState1D& ur,const IdealGasEquationOfState& eos){
    const double pl=eos.pressure(ul.density,ul.momentum,ul.energy);
    const double pr=eos.pressure(ur.density,ur.momentum,ur.energy);
    const double vl=ul.momentum/ul.density;
    const double vr=ur.momentum/ur.density;
    const double a=std::max(std::abs(vl)+eos.sound_speed(ul.density,pl),std::abs(vr)+eos.sound_speed(ur.density,pr));
    const auto fl=flux(ul,eos);
    const auto fr=flux(ur,eos);
    return add(mul(add(fl,fr),0.5),add(ur,ul,-1.0),-0.5*a);
}

std::size_t wrapped_index(long long index,std::size_t n,bool periodic){
    if(periodic){
        const auto nn=static_cast<long long>(n);
        long long value=index%nn;
        if(value<0)value+=nn;
        return static_cast<std::size_t>(value);
    }
    if(index<0)return 0U;
    const auto u=static_cast<std::size_t>(index);
    return std::min(u,n-1U);
}

} // namespace

double IdealGasEquationOfState::pressure(double rho,double mom,double e)const{
    if(!(gamma>1.0)||!(gas_constant>0.0)||!(rho>0.0))throw std::invalid_argument("invalid ideal-gas state");
    const double p=(gamma-1.0)*(e-0.5*mom*mom/rho);
    if(!(p>0.0)||!std::isfinite(p))throw std::runtime_error("non-positive ideal-gas pressure");
    return p;
}

double IdealGasEquationOfState::sound_speed(double rho,double p)const{
    if(!(rho>0.0)||!(p>0.0))throw std::invalid_argument("invalid sound-speed state");
    return std::sqrt(gamma*p/rho);
}

double IdealGasEquationOfState::temperature(double rho,double p)const{
    if(!(rho>0.0)||!(p>0.0))throw std::invalid_argument("invalid temperature state");
    return p/(rho*gas_constant);
}

CompressibleEuler1D::CompressibleEuler1D(Compressible1DConfig c):config_(c),state_(c.cells),next_(c.cells){
    if(c.cells<8U||!(c.length>0.0)||!(c.cfl>0.0&&c.cfl<1.0)||!(c.eos.gamma>1.0)||
       !(c.eos.gas_constant>0.0)||!(c.weno_epsilon>0.0)||!(c.density_floor>0.0)||!(c.pressure_floor>0.0))
        throw std::invalid_argument("invalid compressible 1-D controls");
    if(c.reconstruction==CompressibleReconstruction::characteristic_weno5&&c.cells<8U)
        throw std::invalid_argument("WENO5 requires at least eight cells");
    initialize_uniform({1.0,0.0,1.0});
}

void CompressibleEuler1D::initialize(EulerPrimitive1D left,EulerPrimitive1D right,double fraction){
    if(!(fraction>0.0&&fraction<1.0))throw std::invalid_argument("invalid Riemann interface");
    const auto ul=conservative(left,config_.eos),ur=conservative(right,config_.eos);
    for(std::size_t i=0;i<state_.size();++i)
        state_[i]=(static_cast<double>(i)+0.5)/static_cast<double>(state_.size())<fraction?ul:ur;
}

void CompressibleEuler1D::initialize_uniform(EulerPrimitive1D q){
    const auto u=conservative(q,config_.eos);
    std::fill(state_.begin(),state_.end(),u);
}

void CompressibleEuler1D::initialize_profile(std::span<const EulerPrimitive1D> states){
    if(states.size()!=state_.size())throw std::invalid_argument("compressible profile size must match cell count");
    for(std::size_t i=0;i<states.size();++i)state_[i]=conservative(states[i],config_.eos);
}

EulerPrimitive1D CompressibleEuler1D::primitive(std::size_t i)const{
    const auto& u=state_.at(i);
    return {u.density,u.momentum/u.density,config_.eos.pressure(u.density,u.momentum,u.energy)};
}

bool CompressibleEuler1D::physical(const EulerState1D& u)const noexcept{
    if(!(u.density>=config_.density_floor)||!std::isfinite(u.density)||!std::isfinite(u.momentum)||!std::isfinite(u.energy))return false;
    const double p=raw_pressure(u,config_.eos);
    return p>=config_.pressure_floor&&std::isfinite(p);
}

void CompressibleEuler1D::require_physical(std::span<const EulerState1D> state)const{
    for(const auto& u:state)if(!physical(u))throw std::runtime_error("compressible update produced a non-physical state");
}

double CompressibleEuler1D::stable_timestep()const{
    double speed=0.0;
    for(std::size_t i=0;i<state_.size();++i){
        const auto q=primitive(i);
        speed=std::max(speed,std::abs(q.velocity)+config_.eos.sound_speed(q.density,q.pressure));
    }
    if(!(speed>0.0)||!std::isfinite(speed))throw std::runtime_error("invalid compressible characteristic speed");
    return config_.cfl*(config_.length/static_cast<double>(state_.size()))/speed;
}

std::vector<EulerState1D> CompressibleEuler1D::spatial_residual(std::span<const EulerState1D> state)const{
    const std::size_t n=state.size();
    const double dx=config_.length/static_cast<double>(n);
    std::vector<EulerState1D> face(n+1U);

    if(config_.reconstruction==CompressibleReconstruction::muscl_minmod){
        std::vector<EulerState1D> slope(n);
        for(std::size_t i=0;i<n;++i){
            if(!config_.periodic&&(i==0U||i+1U==n)){slope[i]={};continue;}
            const std::size_t im=i?i-1U:n-1U,ip=i+1U<n?i+1U:0U;
            const auto& a=state[im];const auto& b=state[i];const auto& c=state[ip];
            slope[i]={minmod(b.density-a.density,c.density-b.density),minmod(b.momentum-a.momentum,c.momentum-b.momentum),minmod(b.energy-a.energy,c.energy-b.energy)};
        }
        for(std::size_t f=0;f<=n;++f){
            const std::size_t li=f==0U?(config_.periodic?n-1U:0U):f-1U;
            const std::size_t ri=f==n?(config_.periodic?0U:n-1U):f;
            EulerState1D ul=add(state[li],slope[li],0.5),ur=add(state[ri],slope[ri],-0.5);
            if(!physical(ul))ul=state[li];
            if(!physical(ur))ur=state[ri];
            face[f]=rusanov_flux(ul,ur,config_.eos);
        }
    }else{
        for(std::size_t f=0;f<=n;++f){
            const auto sample=[&](long long index)->const EulerState1D&{return state[wrapped_index(index,n,config_.periodic)];};
            const long long ff=static_cast<long long>(f);
            const EulerState1D& center_left=sample(ff-1LL);
            const EulerState1D& center_right=sample(ff);
            EulerState1D ul=center_left,ur=center_right;
            try{
                const auto basis=roe_basis(center_left,center_right,config_.eos);
                std::array<Vector3,5> left_characteristic{};
                std::array<Vector3,5> right_characteristic{};
                for(std::size_t s=0;s<5U;++s){
                    const long long offset=static_cast<long long>(s);
                    left_characteristic[s]=multiply(basis.left,vector_of(sample(ff-3LL+offset)));
                    right_characteristic[s]=multiply(basis.left,vector_of(sample(ff+2LL-offset)));
                }
                Vector3 wl{},wr{};
                for(std::size_t k=0;k<3U;++k){
                    std::array<double,5> ql{},qr{};
                    for(std::size_t s=0;s<5U;++s){ql[s]=left_characteristic[s][k];qr[s]=right_characteristic[s][k];}
                    wl[k]=weno5_left(ql,config_.weno_epsilon);
                    wr[k]=weno5_left(qr,config_.weno_epsilon);
                }
                ul=state_of(multiply(basis.right,wl));
                ur=state_of(multiply(basis.right,wr));
            }catch(const std::exception&){
                ul=center_left;ur=center_right;
            }
            if(!physical(ul))ul=center_left;
            if(!physical(ur))ur=center_right;
            face[f]=rusanov_flux(ul,ur,config_.eos);
        }
    }

    std::vector<EulerState1D> residual(n);
    for(std::size_t i=0;i<n;++i)residual[i]=mul(add(face[i+1U],face[i],-1.0),-1.0/dx);
    return residual;
}

double CompressibleEuler1D::step(double dt){
    if(dt==0.0)dt=stable_timestep();
    if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument("invalid compressible timestep");

    const auto advance=[&](std::span<const EulerState1D> base,double scale){
        const auto rhs=spatial_residual(base);
        std::vector<EulerState1D> out(base.size());
        for(std::size_t i=0;i<base.size();++i)out[i]=add(base[i],rhs[i],scale);
        require_physical(out);
        return out;
    };

    if(config_.time_integrator==CompressibleTimeIntegrator::forward_euler){
        next_=advance(state_,dt);
    }else{
        const std::vector<EulerState1D> u0=state_;
        const auto u1=advance(u0,dt);
        const auto u1e=advance(u1,dt);
        std::vector<EulerState1D> u2(state_.size());
        for(std::size_t i=0;i<u2.size();++i)u2[i]=add(mul(u0[i],0.75),mul(u1e[i],0.25));
        require_physical(u2);
        const auto u2e=advance(u2,dt);
        for(std::size_t i=0;i<next_.size();++i)next_[i]=add(mul(u0[i],1.0/3.0),mul(u2e[i],2.0/3.0));
        require_physical(next_);
    }
    state_.swap(next_);
    return dt;
}

void CompressibleEuler1D::run(double duration){
    if(!(duration>=0.0)||!std::isfinite(duration))throw std::invalid_argument("invalid compressible duration");
    double t=0.0;
    while(t<duration){
        const double dt=std::min(stable_timestep(),duration-t);
        t+=step(dt);
    }
}

std::vector<double> CompressibleEuler1D::pressure_jump_sensor()const{
    const std::size_t n=state_.size();
    std::vector<double> pressure(n);
    for(std::size_t i=0;i<n;++i)pressure[i]=primitive(i).pressure;
    std::vector<double> sensor(n,0.0);
    for(std::size_t i=0;i<n;++i){
        const auto normalized_jump=[&](std::size_t j){
            const double denom=pressure[i]+pressure[j]+config_.pressure_floor;
            return std::abs(pressure[j]-pressure[i])/denom;
        };
        if(i>0U)sensor[i]=std::max(sensor[i],normalized_jump(i-1U));
        else if(config_.periodic)sensor[i]=std::max(sensor[i],normalized_jump(n-1U));
        if(i+1U<n)sensor[i]=std::max(sensor[i],normalized_jump(i+1U));
        else if(config_.periodic)sensor[i]=std::max(sensor[i],normalized_jump(0U));
    }
    return sensor;
}

double CompressibleEuler1D::mass()const{
    double sum=0.0;const double dx=config_.length/static_cast<double>(state_.size());
    for(const auto& u:state_)sum+=u.density*dx;
    return sum;
}

double CompressibleEuler1D::momentum()const{
    double sum=0.0;const double dx=config_.length/static_cast<double>(state_.size());
    for(const auto& u:state_)sum+=u.momentum*dx;
    return sum;
}

double CompressibleEuler1D::total_energy()const{
    double sum=0.0;const double dx=config_.length/static_cast<double>(state_.size());
    for(const auto& u:state_)sum+=u.energy*dx;
    return sum;
}

} // namespace cfd::fvm
