#include "cfd/multibody/resident_cfd_dem_sycl.hpp"

#if defined(CFD_HAS_SYCL)

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::multibody {
namespace {
using AtomicFloat = sycl::atomic_ref<float, sycl::memory_order::relaxed,
    sycl::memory_scope::device, sycl::access::address_space::global_space>;

template<class T> void free_if(T*& p, sycl::queue& q) noexcept { if(p){sycl::free(p,q);p=nullptr;} }
}

ResidentCfdDemSyclCoupler::ResidentCfdDemSyclCoupler(cfd::fvm::ResidentPolyMeshSycl& mesh,
                                                       ResidentDemSycl& dem,
                                                       ResidentCfdDemSyclConfig config)
    : mesh_(mesh), dem_(dem), config_(config) {
    if(!(config_.fluid_density>0.0)||!(config_.dynamic_viscosity>0.0)||
       !(config_.minimum_void_fraction>0.0&&config_.minimum_void_fraction<=1.0)||
       !(config_.void_drag_exponent>=0.0)||!(config_.saffman_coefficient>=0.0)) {
        throw std::invalid_argument("invalid resident CFD/DEM coupling configuration");
    }
    allocate();
}

ResidentCfdDemSyclCoupler::~ResidentCfdDemSyclCoupler() noexcept { release(); }

void ResidentCfdDemSyclCoupler::allocate() {
    auto& q=mesh_.queue();const auto n=mesh_.cell_count();
    alpha_accum_=sycl::malloc_device<float>(n,q);
    reaction_x_=sycl::malloc_device<float>(n,q);reaction_y_=sycl::malloc_device<float>(n,q);reaction_z_=sycl::malloc_device<float>(n,q);
    solid_fraction_=sycl::malloc_device<double>(n,q);void_fraction_=sycl::malloc_device<double>(n,q);
    if(!alpha_accum_||!reaction_x_||!reaction_y_||!reaction_z_||!solid_fraction_||!void_fraction_){release();throw std::bad_alloc{};}
}

void ResidentCfdDemSyclCoupler::release() noexcept {
    auto& q=mesh_.queue();try{q.wait_and_throw();}catch(...){}
    free_if(alpha_accum_,q);free_if(reaction_x_,q);free_if(reaction_y_,q);free_if(reaction_z_,q);
    free_if(solid_fraction_,q);free_if(void_fraction_,q);
}

std::size_t ResidentCfdDemSyclCoupler::resident_bytes() const noexcept {
    return mesh_.cell_count()*(4U*sizeof(float)+2U*sizeof(double));
}

void ResidentCfdDemSyclCoupler::couple(const double* fluid_velocity,
                                        const double* pressure_gradient,
                                        const double* vorticity,
                                        double* fluid_acceleration) {
    if(!fluid_velocity) throw std::invalid_argument("resident CFD/DEM coupling requires fluid velocity");
    auto& q=mesh_.queue();const auto mv=mesh_.device_view();const auto dv=dem_.device_view();
    const std::size_t nc=mv.cell_count,np=dv.count;
    q.memset(alpha_accum_,0,nc*sizeof(float));q.memset(reaction_x_,0,nc*sizeof(float));q.memset(reaction_y_,0,nc*sizeof(float));q.memset(reaction_z_,0,nc*sizeof(float));
    const auto* cg=mv.cell_geometry;float* alpha=alpha_accum_;const double pi=std::numbers::pi;
    const float* px=dv.x;const float* py=dv.y;const float* pz=dv.z;const float* pr=dv.radius;
    q.parallel_for(sycl::range<1>(np),[=](sycl::id<1> id){
        const std::size_t i=id[0];std::size_t best=0U;double bestd=1.0e300;
        for(std::size_t c=0;c<nc;++c){const double dx=cg[c]-px[i],dy=cg[nc+c]-py[i],dz=cg[2U*nc+c]-pz[i];const double d2=dx*dx+dy*dy+dz*dz;if(d2<bestd){bestd=d2;best=c;}}
        const double r=pr[i],volume=(4.0/3.0)*pi*r*r*r,cell_volume=cg[3U*nc+best];
        AtomicFloat(alpha[best]).fetch_add(static_cast<float>(volume/cell_volume));
    });
    const double min_void=config_.minimum_void_fraction;double* solid=solid_fraction_;double* eps=void_fraction_;
    q.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){const auto c=id[0];double a=alpha[c];const double amax=1.0-min_void;if(a<0.0)a=0.0;if(a>amax)a=amax;solid[c]=a;eps[c]=1.0-a;});

    const double rho=config_.fluid_density,mu=config_.dynamic_viscosity,nu=mu/rho,void_exp=config_.void_drag_exponent,saff=config_.saffman_coefficient;
    const bool with_pg=config_.include_pressure_gradient&&pressure_gradient;const bool with_lift=config_.include_saffman_lift&&vorticity;
    float *pfx=dv.fx,*pfy=dv.fy,*pfz=dv.fz;const float *pvx=dv.vx,*pvy=dv.vy,*pvz=dv.vz;
    float *rx=reaction_x_,*ry=reaction_y_,*rz=reaction_z_;
    q.parallel_for(sycl::range<1>(np),[=](sycl::id<1> id){
        const std::size_t i=id[0];std::size_t c0=0U;double bestd=1.0e300;
        for(std::size_t c=0;c<nc;++c){const double dx=cg[c]-px[i],dy=cg[nc+c]-py[i],dz=cg[2U*nc+c]-pz[i];const double d2=dx*dx+dy*dy+dz*dz;if(d2<bestd){bestd=d2;c0=c;}}
        const double ux=fluid_velocity[c0],uy=fluid_velocity[nc+c0],uz=fluid_velocity[2U*nc+c0];
        const double rvx=ux-pvx[i],rvy=uy-pvy[i],rvz=uz-pvz[i],speed=sycl::sqrt(rvx*rvx+rvy*rvy+rvz*rvz);
        const double d=2.0*pr[i],re=rho*speed*d/mu;
        double cd=0.0;if(speed>0.0){cd=re<1000.0?(24.0/(re>1.0e-30?re:1.0e-30))*(1.0+0.15*sycl::pow(re,0.687)):0.44;}
        const double area=0.25*pi*d*d,ev=eps[c0]>1.0e-12?eps[c0]:1.0e-12,hindrance=sycl::pow(ev,-void_exp);
        double fx=0.5*cd*rho*area*speed*hindrance*rvx,fy=0.5*cd*rho*area*speed*hindrance*rvy,fz=0.5*cd*rho*area*speed*hindrance*rvz;
        const double volume=(pi/6.0)*d*d*d;
        if(with_pg){fx-=volume*pressure_gradient[c0];fy-=volume*pressure_gradient[nc+c0];fz-=volume*pressure_gradient[2U*nc+c0];}
        if(with_lift){
            const double ox=vorticity[c0],oy=vorticity[nc+c0],oz=vorticity[2U*nc+c0],om=sycl::sqrt(ox*ox+oy*oy+oz*oz);
            if(om>0.0){const double scale=saff*rho*d*d*sycl::sqrt(nu*om)/om;fx+=(rvy*oz-rvz*oy)*scale;fy+=(rvz*ox-rvx*oz)*scale;fz+=(rvx*oy-rvy*ox)*scale;}
        }
        pfx[i]+=static_cast<float>(fx);pfy[i]+=static_cast<float>(fy);pfz[i]+=static_cast<float>(fz);
        AtomicFloat(rx[c0]).fetch_add(static_cast<float>(-fx));AtomicFloat(ry[c0]).fetch_add(static_cast<float>(-fy));AtomicFloat(rz[c0]).fetch_add(static_cast<float>(-fz));
    });

    if(fluid_acceleration){
        const double invrho=1.0/rho;
        q.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){const auto c=id[0];const double scale=invrho/cg[3U*nc+c];fluid_acceleration[c]+=rx[c]*scale;fluid_acceleration[nc+c]+=ry[c]*scale;fluid_acceleration[2U*nc+c]+=rz[c]*scale;});
    }
}

void ResidentCfdDemSyclCoupler::download_solid_volume_fraction(std::vector<double>& out) const {
    out.resize(mesh_.cell_count());mesh_.download_cell_scalar(solid_fraction_,out);
}
void ResidentCfdDemSyclCoupler::download_void_fraction(std::vector<double>& out) const {
    out.resize(mesh_.cell_count());mesh_.download_cell_scalar(void_fraction_,out);
}

} // namespace cfd::multibody
#endif
