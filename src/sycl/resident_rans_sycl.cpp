#include "cfd/fvm/resident_rans_sycl.hpp"

#if defined(CFD_HAS_SYCL)

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace cfd::fvm {
namespace {

cfd::core::CsrMatrix cell_pattern(const PolyMesh& mesh) {
    cfd::core::CsrBuilder b(mesh.cell_count(), mesh.cell_count());
    for (std::size_t c=0;c<mesh.cell_count();++c) b.add(c,c,1.0);
    for (const auto& f:mesh.faces()) if(!f.boundary()) { b.add(f.owner,f.neighbour,-1.0); b.add(f.neighbour,f.owner,-1.0); }
    return b.build();
}

template<class T> void free_if(T*& p, sycl::queue& q) noexcept { if(p){sycl::free(p,q);p=nullptr;} }

void validate_rans_controls(const RansTransportControls& c) {
    if(!(c.dt>0.0)||!std::isfinite(c.dt)||c.linear_iterations==0U||!(c.linear_tolerance>0.0)||!std::isfinite(c.linear_tolerance))
        throw std::invalid_argument("invalid resident RANS transport controls");
    if(c.temporal_scheme!=TemporalScheme::euler)
        throw std::invalid_argument("resident RANS currently supports Euler time integration");
}

std::uint8_t boundary_kind(TurbulenceScalarBoundaryType t) {
    return static_cast<std::uint8_t>(t==TurbulenceScalarBoundaryType::fixedValue
        ? PressureBoundaryType::fixedValue : PressureBoundaryType::zeroGradient);
}

} // namespace

ResidentFvmFieldRegistrySycl::~ResidentFvmFieldRegistrySycl() noexcept {
    auto& q=mesh_.queue();try{q.wait_and_throw();}catch(...){}
    for(auto& [name,e]:fields_){(void)name;if(e.data) sycl::free(e.data,q);}
}

double* ResidentFvmFieldRegistrySycl::create(std::string name,std::size_t components) {
    if(name.empty()||components==0U) throw std::invalid_argument("resident field name/components invalid");
    if(fields_.contains(name)) throw std::invalid_argument("resident field already exists: "+name);
    auto& q=mesh_.queue();const std::size_t count=components*mesh_.cell_count();double* p=sycl::malloc_device<double>(count,q);
    if(!p) throw std::bad_alloc{};
    fields_.emplace(std::move(name),Entry{p,components});
    resident_bytes_+=count*sizeof(double);
    return p;
}

double* ResidentFvmFieldRegistrySycl::create_scalar(std::string name,double initial) {
    double* p=create(std::move(name),1U);mesh_.fill(p,mesh_.cell_count(),initial);return p;
}

double* ResidentFvmFieldRegistrySycl::create_vector(std::string name,Vec3 initial) {
    double* p=create(std::move(name),3U);const auto n=mesh_.cell_count();mesh_.fill(p,n,initial.x);mesh_.fill(p+n,n,initial.y);mesh_.fill(p+2U*n,n,initial.z);return p;
}

double* ResidentFvmFieldRegistrySycl::create_components(std::string name,std::size_t components,double initial) {
    double* p=create(std::move(name),components);mesh_.fill(p,components*mesh_.cell_count(),initial);return p;
}

bool ResidentFvmFieldRegistrySycl::contains(std::string_view name) const { return fields_.find(std::string(name))!=fields_.end(); }
double* ResidentFvmFieldRegistrySycl::scalar(std::string_view name) { auto it=fields_.find(std::string(name));if(it==fields_.end()||it->second.components!=1U) throw std::out_of_range("resident scalar field not found");return it->second.data; }
const double* ResidentFvmFieldRegistrySycl::scalar(std::string_view name) const { auto it=fields_.find(std::string(name));if(it==fields_.end()||it->second.components!=1U) throw std::out_of_range("resident scalar field not found");return it->second.data; }
double* ResidentFvmFieldRegistrySycl::vector(std::string_view name) { auto it=fields_.find(std::string(name));if(it==fields_.end()||it->second.components!=3U) throw std::out_of_range("resident vector field not found");return it->second.data; }

ResidentScalarEquationSycl::ResidentScalarEquationSycl(const PolyMesh& host_mesh,ResidentFvmFieldRegistrySycl& registry,
    std::string prefix,double dt,std::size_t linear_iterations,double linear_tolerance)
    : registry_(registry),mesh_(registry.mesh()),prefix_(std::move(prefix)),dt_(dt),linear_iterations_(linear_iterations),linear_tolerance_(linear_tolerance),
      patch_names_(host_mesh.patches().size()),patch_faces_(host_mesh.patches().size()),
      boundary_kind_host_(host_mesh.face_count(),static_cast<std::uint8_t>(PressureBoundaryType::zeroGradient)),boundary_fixed_host_(host_mesh.face_count(),0.0) {
    if(host_mesh.cell_count()!=mesh_.cell_count()||host_mesh.face_count()!=mesh_.face_count()) throw std::invalid_argument("resident scalar host/device mesh mismatch");
    if(prefix_.empty()||!(dt_>0.0)||!std::isfinite(dt_)||linear_iterations_==0U||!(linear_tolerance_>0.0)||!std::isfinite(linear_tolerance_))
        throw std::invalid_argument("invalid shared resident scalar controls");
    for(std::size_t p=0;p<host_mesh.patches().size();++p) patch_names_[p]=host_mesh.patches()[p].name;
    for(std::size_t f=0;f<host_mesh.face_count();++f) if(host_mesh.faces()[f].boundary()) patch_faces_[host_mesh.faces()[f].patch].push_back(f);
    solver_=std::make_unique<cfd::core::SyclCsrLinearAlgebra>(cell_pattern(host_mesh),mesh_.queue());
    field_=registry_.create_scalar(prefix_+".value",0.0);old_=registry_.create_scalar(prefix_+".old",0.0);
    diffusivity_=registry_.create_scalar(prefix_+".diffusivity",0.0);source_=registry_.create_scalar(prefix_+".source",0.0);
    sink_=registry_.create_scalar(prefix_+".sink",0.0);rhs_=registry_.create_scalar(prefix_+".rhs",0.0);
    auto& q=mesh_.queue();matrix_values_=sycl::malloc_device<double>(solver_->nonzeros(),q);
    boundary_kind_=sycl::malloc_device<std::uint8_t>(mesh_.face_count(),q);boundary_fixed_=sycl::malloc_device<double>(mesh_.face_count(),q);
    if(!matrix_values_||!boundary_kind_||!boundary_fixed_) throw std::bad_alloc{};
    upload_boundary_state();
}

ResidentScalarEquationSycl::~ResidentScalarEquationSycl() noexcept { auto& q=mesh_.queue();try{q.wait_and_throw();}catch(...){}free_if(matrix_values_,q);free_if(boundary_kind_,q);free_if(boundary_fixed_,q);solver_.reset(); }

void ResidentScalarEquationSycl::upload_boundary_state() {
    auto& q=mesh_.queue();const auto nf=mesh_.face_count();const auto kb=nf*sizeof(std::uint8_t),vb=nf*sizeof(double);
    q.memcpy(boundary_kind_,boundary_kind_host_.data(),kb);q.memcpy(boundary_fixed_,boundary_fixed_host_.data(),vb).wait_and_throw();
    transfer_stats_.record_host_to_device(kb+vb);transfer_stats_.record_synchronization();
}

void ResidentScalarEquationSycl::set_boundary(std::string_view patch,TurbulenceScalarBoundaryType type,double value) {
    if(!std::isfinite(value)) throw std::invalid_argument("resident scalar boundary must be finite");
    std::size_t pi=patch_names_.size();for(std::size_t p=0;p<patch_names_.size();++p) if(patch_names_[p]==patch){pi=p;break;}
    if(pi==patch_names_.size()) throw std::out_of_range("resident scalar patch not found");
    for(auto f:patch_faces_[pi]){boundary_kind_host_[f]=boundary_kind(type);boundary_fixed_host_[f]=value;}upload_boundary_state();
}

void ResidentScalarEquationSycl::initialize_uniform(double value) { if(!std::isfinite(value)) throw std::invalid_argument("resident scalar initial value must be finite");mesh_.fill(field_,mesh_.cell_count(),value);mesh_.fill(old_,mesh_.cell_count(),value);time_=0.0;steps_=0U; }

cfd::core::IterativeSolverResult ResidentScalarEquationSycl::step(const double* flux,double minimum,double maximum) {
    if(!flux) throw std::invalid_argument("resident scalar equation requires face flux");
    const auto bytes=mesh_.cell_count()*sizeof(double);
    mesh_.queue().memcpy(old_,field_,bytes);transfer_stats_.record_device_to_device(bytes);
    mesh_.assemble_scalar_transport_system_variable(solver_->row_offsets_device(),solver_->column_indices_device(),solver_->nonzeros(),old_,flux,boundary_kind_,boundary_fixed_,diffusivity_,dt_,source_,sink_,matrix_values_,rhs_);
    solver_->update_values_device(matrix_values_);auto result=solver_->bicgstab_device(rhs_,field_,linear_iterations_,linear_tolerance_);
    if(result.converged){mesh_.clamp_scalar(field_,minimum,maximum);time_+=dt_;++steps_;}else{mesh_.queue().memcpy(field_,old_,bytes);transfer_stats_.record_device_to_device(bytes);}return result;
}
void ResidentScalarEquationSycl::download(std::span<double> host) const { mesh_.download_cell_scalar(field_,host); }
void ResidentScalarEquationSycl::reset_transfer_stats() const noexcept { transfer_stats_.reset();mesh_.reset_transfer_stats();if(solver_)solver_->reset_transfer_stats(); }
std::uint64_t ResidentScalarEquationSycl::hot_loop_host_transfer_bytes() const noexcept { return transfer_stats_.host_transfer_bytes()+mesh_.transfer_stats().host_transfer_bytes()+(solver_?solver_->transfer_stats().host_transfer_bytes():0U); }

ResidentSpalartAllmarasSycl::ResidentSpalartAllmarasSycl(const PolyMesh& host,ResidentFvmFieldRegistrySycl& r,SpalartAllmarasConfig c,std::string p)
    : registry_(r),config_(c),equation_(host,r,p+".nuTilde",c.transport.dt,c.transport.linear_iterations,c.transport.linear_tolerance) {
    validate_rans_controls(c.transport);if(!(c.molecular_viscosity>0.0)||!(c.sigma_nu_tilde>0.0)||!(c.kappa>0.0)||!(c.minimum_wall_distance>0.0)) throw std::invalid_argument("invalid resident SA config");
    wall_distance_=r.create_scalar(p+".wallDistance",1.0);strain_=r.create_scalar(p+".strain",0.0);velocity_gradient_=r.create_components(p+".gradU",9U,0.0);
    scalar_gradient_=r.create_vector(p+".gradNuTilde",{});nut_=r.create_scalar(p+".nut",0.0);
}
void ResidentSpalartAllmarasSycl::initialize_uniform(double v){ if(v<config_.minimum_nu_tilde) throw std::invalid_argument("resident SA initial value below floor");equation_.initialize_uniform(v); }
void ResidentSpalartAllmarasSycl::set_boundary(std::string_view p,TurbulenceScalarBoundaryType t,double v){equation_.set_boundary(p,t,v);}
void ResidentSpalartAllmarasSycl::set_wall_distance(double d){if(!(d>=config_.minimum_wall_distance))throw std::invalid_argument("resident SA wall distance invalid");registry_.mesh().fill(wall_distance_,registry_.mesh().cell_count(),d);}

cfd::core::IterativeSolverResult ResidentSpalartAllmarasSycl::step(const double* flux,const double* velocity) {
    auto& m=registry_.mesh();m.velocity_strain_rate_magnitude(velocity,velocity_gradient_,strain_);m.gauss_gradient_scalar(equation_.field_device(),scalar_gradient_,nullptr);
    const auto n=m.cell_count();auto* nt=equation_.field_device();auto* diff=equation_.diffusivity_device();auto* src=equation_.source_device();auto* sink=equation_.sink_device();
    auto* wall=wall_distance_;auto* strain=strain_;auto* grad=scalar_gradient_;auto* nut=nut_;const auto c=config_;
    m.queue().parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const std::size_t i=id[0];const double x=nt[i]>c.minimum_nu_tilde?nt[i]:c.minimum_nu_tilde;
        const double d=wall[i]>c.minimum_wall_distance?wall[i]:c.minimum_wall_distance;const double chi=x/c.molecular_viscosity;const double chi3=chi*chi*chi,cv3=c.cv1*c.cv1*c.cv1;
        const double fv1=chi3/(chi3+cv3);const double fv2=1.0-chi/(1.0+chi*fv1);const double kd2=c.kappa*c.kappa*d*d;
        const double raw=strain[i]+x*fv2/kd2;const double st=raw>c.cs*strain[i]?raw:c.cs*strain[i];const double den=st*kd2>1.0e-30?st*kd2:1.0e-30;
        double rr=x/den;if(rr<0.0)rr=0.0;if(rr>10.0)rr=10.0;const double r2=rr*rr,r3=r2*rr,r6=r3*r3;const double g=rr+c.cw2*(r6-rr);
        const double g2=g*g,g3=g2*g,g6=g3*g3,cw32=c.cw3*c.cw3,cw36=cw32*cw32*cw32;const double ratio=(1.0+cw36)/(g6+cw36);
        double cube=ratio>1.0?ratio:1.0;for(int it=0;it<6;++it) cube=(2.0*cube+ratio/(cube*cube))/3.0;const double fw=g*sycl::sqrt(cube);
        const double cw1=c.cb1/(c.kappa*c.kappa)+(1.0+c.cb2)/c.sigma_nu_tilde;const double gx=grad[i],gy=grad[n+i],gz=grad[2U*n+i];
        diff[i]=(c.molecular_viscosity+x)/c.sigma_nu_tilde;src[i]=c.cb1*st*x+(c.cb2/c.sigma_nu_tilde)*(gx*gx+gy*gy+gz*gz);sink[i]=cw1*fw*x/(d*d);
        nut[i]=x*fv1;});
    return equation_.step(flux,config_.minimum_nu_tilde,std::numeric_limits<double>::max());
}

ResidentKEpsilonSycl::ResidentKEpsilonSycl(const PolyMesh& host,ResidentFvmFieldRegistrySycl& r,KEpsilonConfig c,std::string p)
    : registry_(r),config_(c),k_(host,r,p+".k",c.transport.dt,c.transport.linear_iterations,c.transport.linear_tolerance),epsilon_(host,r,p+".epsilon",c.transport.dt,c.transport.linear_iterations,c.transport.linear_tolerance) {
    validate_rans_controls(c.transport);if(!(c.molecular_viscosity>0.0)||!(c.c_mu>0.0)||!(c.minimum_k>0.0)||!(c.minimum_epsilon>0.0)) throw std::invalid_argument("invalid resident k-epsilon config");
    strain_=r.create_scalar(p+".strain",0.0);velocity_gradient_=r.create_components(p+".gradU",9U,0.0);nut_=r.create_scalar(p+".nut",0.0);
}
void ResidentKEpsilonSycl::initialize_uniform(double k,double e){if(k<config_.minimum_k||e<config_.minimum_epsilon)throw std::invalid_argument("resident k-epsilon initial values below floors");k_.initialize_uniform(k);epsilon_.initialize_uniform(e);}
void ResidentKEpsilonSycl::set_k_boundary(std::string_view p,TurbulenceScalarBoundaryType t,double v){k_.set_boundary(p,t,v);}void ResidentKEpsilonSycl::set_epsilon_boundary(std::string_view p,TurbulenceScalarBoundaryType t,double v){epsilon_.set_boundary(p,t,v);}
TwoEquationLinearResult ResidentKEpsilonSycl::step(const double* flux,const double* velocity){auto& m=registry_.mesh();m.velocity_strain_rate_magnitude(velocity,velocity_gradient_,strain_);const auto n=m.cell_count();
    auto* kk=k_.field_device();auto* ee=epsilon_.field_device();auto* kd=k_.diffusivity_device();auto* ed=epsilon_.diffusivity_device();auto* ks=k_.source_device();auto* es=epsilon_.source_device();auto* ksi=k_.sink_device();auto* esi=epsilon_.sink_device();auto* nut=nut_;auto* strain=strain_;const auto c=config_;
    m.queue().parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const auto i=id[0];const double k=kk[i]>c.minimum_k?kk[i]:c.minimum_k,e=ee[i]>c.minimum_epsilon?ee[i]:c.minimum_epsilon;
        double nu=c.c_mu*k*k/e;const double maxnu=c.maximum_eddy_viscosity_ratio*c.molecular_viscosity;if(nu>maxnu)nu=maxnu;nut[i]=nu;const double prod0=nu*strain[i]*strain[i],lim=c.production_limiter*e,prod=prod0<lim?prod0:lim;
        kd[i]=c.molecular_viscosity+nu/c.sigma_k;ed[i]=c.molecular_viscosity+nu/c.sigma_epsilon;ks[i]=prod;ksi[i]=e/k;es[i]=c.c1*(e/k)*prod;esi[i]=c.c2*e/k;});
    TwoEquationLinearResult out{};out.first=k_.step(flux,c.minimum_k,std::numeric_limits<double>::max());out.second=epsilon_.step(flux,c.minimum_epsilon,std::numeric_limits<double>::max());return out;}

ResidentKOmegaSSTSycl::ResidentKOmegaSSTSycl(const PolyMesh& host,ResidentFvmFieldRegistrySycl& r,KOmegaSSTConfig c,std::string p)
    : registry_(r),config_(c),k_(host,r,p+".k",c.transport.dt,c.transport.linear_iterations,c.transport.linear_tolerance),omega_(host,r,p+".omega",c.transport.dt,c.transport.linear_iterations,c.transport.linear_tolerance) {
    validate_rans_controls(c.transport);if(!(c.molecular_viscosity>0.0)||!(c.beta_star>0.0)||!(c.a1>0.0)||!(c.minimum_k>0.0)||!(c.minimum_omega>0.0)) throw std::invalid_argument("invalid resident SST config");
    wall_distance_=r.create_scalar(p+".wallDistance",1.0);grid_scale_=r.create_scalar(p+".gridScale",1.0);strain_=r.create_scalar(p+".strain",0.0);velocity_gradient_=r.create_components(p+".gradU",9U,0.0);
    grad_k_=r.create_vector(p+".gradK",{});grad_omega_=r.create_vector(p+".gradOmega",{});nut_=r.create_scalar(p+".nut",0.0);f1_=r.create_scalar(p+".F1",0.0);f2_=r.create_scalar(p+".F2",0.0);hybrid_=r.create_scalar(p+".hybrid",1.0);
    std::vector<double> delta(host.cell_count());for(std::size_t i=0;i<delta.size();++i)delta[i]=std::cbrt(host.cells()[i].volume);r.mesh().upload_cell_scalar(delta,grid_scale_);
}
void ResidentKOmegaSSTSycl::initialize_uniform(double k,double w){if(k<config_.minimum_k||w<config_.minimum_omega)throw std::invalid_argument("resident SST initial values below floors");k_.initialize_uniform(k);omega_.initialize_uniform(w);}
void ResidentKOmegaSSTSycl::set_k_boundary(std::string_view p,TurbulenceScalarBoundaryType t,double v){k_.set_boundary(p,t,v);}void ResidentKOmegaSSTSycl::set_omega_boundary(std::string_view p,TurbulenceScalarBoundaryType t,double v){omega_.set_boundary(p,t,v);}
void ResidentKOmegaSSTSycl::set_wall_distance(double d){if(!(d>=config_.minimum_wall_distance))throw std::invalid_argument("resident SST wall distance invalid");registry_.mesh().fill(wall_distance_,registry_.mesh().cell_count(),d);}void ResidentKOmegaSSTSycl::set_grid_scale(double d){if(!(d>0.0))throw std::invalid_argument("resident SST grid scale invalid");registry_.mesh().fill(grid_scale_,registry_.mesh().cell_count(),d);}

TwoEquationLinearResult ResidentKOmegaSSTSycl::step(const double* flux,const double* velocity){auto& m=registry_.mesh();m.velocity_strain_rate_magnitude(velocity,velocity_gradient_,strain_);m.gauss_gradient_scalar(k_.field_device(),grad_k_,nullptr);m.gauss_gradient_scalar(omega_.field_device(),grad_omega_,nullptr);const auto n=m.cell_count();
    auto* kk=k_.field_device();auto* ww=omega_.field_device();auto* kd=k_.diffusivity_device();auto* wd=omega_.diffusivity_device();auto* ks=k_.source_device();auto* ws=omega_.source_device();auto* ksi=k_.sink_device();auto* wsi=omega_.sink_device();auto* nut=nut_;auto* f1=f1_;auto* f2=f2_;auto* hy=hybrid_;auto* wall=wall_distance_;auto* delta=grid_scale_;auto* strain=strain_;auto* gk=grad_k_;auto* gw=grad_omega_;const auto c=config_;
    m.queue().parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const std::size_t i=id[0];const double k=kk[i]>c.minimum_k?kk[i]:c.minimum_k,w=ww[i]>c.minimum_omega?ww[i]:c.minimum_omega,d=wall[i]>c.minimum_wall_distance?wall[i]:c.minimum_wall_distance;
        const double crossraw=gk[i]*gw[i]+gk[n+i]*gw[n+i]+gk[2U*n+i]*gw[2U*n+i];const double cd0=2.0*c.alpha_omega2*crossraw/w,cd=cd0>1.0e-20?cd0:1.0e-20;
        const double a=sycl::sqrt(k)/(c.beta_star*w*d),b=500.0*c.molecular_viscosity/(d*d*w),cc=4.0*c.alpha_omega2*k/(cd*d*d);double arg=(a>b?a:b);if(arg>cc)arg=cc;if(arg<0.0)arg=0.0;if(arg>100.0)arg=100.0;const double x=arg*arg*arg*arg;double F1=1.0;if(x<20.0){const double ex=sycl::exp(2.0*x);F1=(ex-1.0)/(ex+1.0);}
        double arg2a=2.0*sycl::sqrt(k)/(c.beta_star*w*d),arg2b=500.0*c.molecular_viscosity/(d*d*w),arg2=arg2a>arg2b?arg2a:arg2b;if(arg2>100.0)arg2=100.0;const double y=arg2*arg2;double F2=1.0;if(y<20.0){const double ey=sycl::exp(2.0*y);F2=(ey-1.0)/(ey+1.0);}f1[i]=F1;f2[i]=F2;
        const double den1=c.a1*w,den2=strain[i]*F2,den=den1>den2?den1:den2;double nu=c.a1*k/(den>1.0e-30?den:1.0e-30);const double maxnu=c.maximum_eddy_viscosity_ratio*c.molecular_viscosity;if(nu>maxnu)nu=maxnu;nut[i]=nu;
        double h=1.0;if(c.des_enabled){const double lt=sycl::sqrt(k)/(c.beta_star*w),ll=c.c_des*(delta[i]>1.0e-30?delta[i]:1.0e-30),raw=lt/ll>1.0?lt/ll:1.0;double shield=1.0;if(c.des_zonal_filter==SSTDESZonalFilter::f1)shield=1.0-F1;else if(c.des_zonal_filter==SSTDESZonalFilter::f2)shield=1.0-F2;if(shield<0.0)shield=0.0;if(shield>1.0)shield=1.0;h=1.0+(raw-1.0)*shield;}hy[i]=h;
        const double ak=F1*c.alpha_k1+(1.0-F1)*c.alpha_k2,aw=F1*c.alpha_omega1+(1.0-F1)*c.alpha_omega2,beta=F1*c.beta1+(1.0-F1)*c.beta2,gamma=F1*c.gamma1+(1.0-F1)*c.gamma2;
        const double p0=nu*strain[i]*strain[i],pl=c.production_limiter*c.beta_star*k*w,p=p0<pl?p0:pl,pr=nu>1.0e-30?p/nu:0.0,cross=2.0*(1.0-F1)*c.alpha_omega2*crossraw/w;
        kd[i]=c.molecular_viscosity+ak*nu;wd[i]=c.molecular_viscosity+aw*nu;ks[i]=p;ksi[i]=c.beta_star*w*h;ws[i]=gamma*pr+(cross>0.0?cross:0.0);wsi[i]=beta*w+(cross<0.0?-cross/w:0.0);});
    TwoEquationLinearResult out{};out.first=k_.step(flux,c.minimum_k,std::numeric_limits<double>::max());out.second=omega_.step(flux,c.minimum_omega,std::numeric_limits<double>::max());return out;}

ResidentThermoSpeciesSourceSycl::ResidentThermoSpeciesSourceSycl(ResidentFvmFieldRegistrySycl& r,std::string p):registry_(r){heat_source_=r.create_scalar(p+".heatSource",0.0);species_source_=r.create_scalar(p+".speciesSource",0.0);rate_=r.create_scalar(p+".rate",0.0);}
void ResidentThermoSpeciesSourceSycl::arrhenius_one_step(const double* temperature,const double* fuel,double A,double Ta,double Q){if(!temperature||!fuel||!(A>=0.0)||!(Ta>=0.0)||!std::isfinite(A)||!std::isfinite(Ta)||!std::isfinite(Q))throw std::invalid_argument("invalid resident Arrhenius coupling controls");const auto n=registry_.mesh().cell_count();auto* rate=rate_;auto* hs=heat_source_;auto* ys=species_source_;
    registry_.mesh().queue().parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const auto i=id[0];const double T=temperature[i]>1.0?temperature[i]:1.0,Y=fuel[i]>0.0?fuel[i]:0.0;const double w=A*Y*sycl::exp(-Ta/T);rate[i]=w;ys[i]=-w;hs[i]=Q*w;});}
void ResidentThermoSpeciesSourceSycl::apply_sources(ResidentScalarEquationSycl& temperature_equation,ResidentScalarEquationSycl& species_equation){
    const auto bytes=registry_.mesh().cell_count()*sizeof(double);auto& q=registry_.mesh().queue();q.memcpy(temperature_equation.source_device(),heat_source_,bytes);q.memcpy(species_equation.source_device(),species_source_,bytes);
}

} // namespace cfd::fvm
#endif
