#include "cfd/multibody/resident_dem_sycl.hpp"

#if defined(CFD_HAS_SYCL)

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>

namespace cfd::multibody {
namespace {
using GlobalAtomicFloat = sycl::atomic_ref<float, sycl::memory_order::relaxed,
    sycl::memory_scope::device, sycl::access::address_space::global_space>;
using GlobalAtomicU32 = sycl::atomic_ref<std::uint32_t, sycl::memory_order::relaxed,
    sycl::memory_scope::device, sycl::access::address_space::global_space>;

std::size_t grid_extent(double lo,double hi,double h) {
    if(!(hi>lo) || !(h>0.0) || !std::isfinite(lo) || !std::isfinite(hi) || !std::isfinite(h)) {
        throw std::invalid_argument("invalid resident DEM domain/grid");
    }
    return std::max<std::size_t>(1U,static_cast<std::size_t>(std::ceil((hi-lo)/h)));
}
}

ResidentDemSycl::ResidentDemSycl(sycl::queue& queue, ResidentDemSyclConfig config)
    : queue_(queue), config_(config) {
    if(config.capacity==0U || config.history_slots_per_particle==0U || !(config.dt>0.0) || !std::isfinite(config.dt)) {
        throw std::invalid_argument("invalid resident DEM capacity/history/dt");
    }
    if(!(config.normal_stiffness>=0.0) || !(config.normal_damping>=0.0) ||
       !(config.tangential_stiffness>=0.0) || !(config.tangential_damping>=0.0) ||
       !(config.friction>=0.0) || !(config.rolling_resistance>=0.0) || !(config.cohesion_force>=0.0)) {
        throw std::invalid_argument("invalid resident DEM contact coefficients");
    }
    nx_=grid_extent(config.domain_min.x,config.domain_max.x,config.cell_size);
    ny_=grid_extent(config.domain_min.y,config.domain_max.y,config.cell_size);
    nz_=grid_extent(config.domain_min.z,config.domain_max.z,config.cell_size);
    if(nx_>std::numeric_limits<std::size_t>::max()/ny_ || nx_*ny_>std::numeric_limits<std::size_t>::max()/nz_) {
        throw std::overflow_error("resident DEM grid too large");
    }
    cell_count_=nx_*ny_*nz_;
    allocate();
}

void ResidentDemSycl::allocate() {
    const auto cap=config_.capacity;
    const auto hist=cap*config_.history_slots_per_particle;
    auto allocf=[&](std::size_t n){ auto* p=sycl::malloc_device<float>(n,queue_); if(!p)throw std::bad_alloc{}; return p; };
    auto allocu=[&](std::size_t n){ auto* p=sycl::malloc_device<std::uint32_t>(n,queue_); if(!p)throw std::bad_alloc{}; return p; };
    try {
        x_=allocf(cap);y_=allocf(cap);z_=allocf(cap);vx_=allocf(cap);vy_=allocf(cap);vz_=allocf(cap);
        wx_=allocf(cap);wy_=allocf(cap);wz_=allocf(cap);radius_=allocf(cap);inverse_mass_=allocf(cap);inverse_inertia_=allocf(cap);
        fx_=allocf(cap);fy_=allocf(cap);fz_=allocf(cap);tx_=allocf(cap);ty_=allocf(cap);tz_=allocf(cap);
        particle_cell_=allocu(cap);sorted_indices_=allocu(cap);
        cell_counts_=allocu(cell_count_);cell_offsets_=allocu(cell_count_+1U);cell_cursor_=allocu(cell_count_);
        history_neighbor_=allocu(hist);history_seen_=allocu(hist);history_tx_=allocf(hist);history_ty_=allocf(hist);history_tz_=allocf(hist);
        queue_.memset(history_neighbor_,0xff,hist*sizeof(std::uint32_t));
        queue_.memset(history_seen_,0,hist*sizeof(std::uint32_t));
        queue_.memset(history_tx_,0,hist*sizeof(float));queue_.memset(history_ty_,0,hist*sizeof(float));queue_.memset(history_tz_,0,hist*sizeof(float));
    } catch(...) { release(); throw; }
}

ResidentDemSycl::~ResidentDemSycl() noexcept { release(); }
void ResidentDemSycl::release() noexcept {
    try { queue_.wait_and_throw(); } catch(...) {}
    auto freep=[&](auto*& p){ if(p){ sycl::free(p,queue_);p=nullptr; } };
    freep(x_);freep(y_);freep(z_);freep(vx_);freep(vy_);freep(vz_);freep(wx_);freep(wy_);freep(wz_);
    freep(radius_);freep(inverse_mass_);freep(inverse_inertia_);freep(fx_);freep(fy_);freep(fz_);freep(tx_);freep(ty_);freep(tz_);
    freep(particle_cell_);freep(cell_counts_);freep(cell_offsets_);freep(cell_cursor_);freep(sorted_indices_);
    freep(history_neighbor_);freep(history_seen_);freep(history_tx_);freep(history_ty_);freep(history_tz_);freep(reduction_scalar_);
}

void ResidentDemSycl::upload(std::span<const ResidentDemParticle> particles) {
    if(particles.size()>config_.capacity) throw std::invalid_argument("resident DEM particle capacity exceeded");
    count_=particles.size();
    std::vector<float> x(count_),y(count_),z(count_),vx(count_),vy(count_),vz(count_),wx(count_),wy(count_),wz(count_);
    std::vector<float> radius(count_),inv_mass(count_),inv_inertia(count_);
    for(std::size_t i=0;i<count_;++i){
        const auto& p=particles[i];
        if(!(p.radius>0.0) || !(p.mass>0.0) || !std::isfinite(p.radius) || !std::isfinite(p.mass) || 2.0*p.radius>config_.cell_size) {
            throw std::invalid_argument("resident DEM requires finite particles with diameter <= cell size");
        }
        if(p.position.x<config_.domain_min.x || p.position.x>=config_.domain_max.x ||
           p.position.y<config_.domain_min.y || p.position.y>=config_.domain_max.y ||
           p.position.z<config_.domain_min.z || p.position.z>=config_.domain_max.z) {
            throw std::invalid_argument("resident DEM particle outside configured domain");
        }
        x[i]=static_cast<float>(p.position.x);y[i]=static_cast<float>(p.position.y);z[i]=static_cast<float>(p.position.z);
        vx[i]=static_cast<float>(p.linear_velocity.x);vy[i]=static_cast<float>(p.linear_velocity.y);vz[i]=static_cast<float>(p.linear_velocity.z);
        wx[i]=static_cast<float>(p.angular_velocity.x);wy[i]=static_cast<float>(p.angular_velocity.y);wz[i]=static_cast<float>(p.angular_velocity.z);
        radius[i]=static_cast<float>(p.radius);inv_mass[i]=static_cast<float>(1.0/p.mass);
        const double inertia=0.4*p.mass*p.radius*p.radius;inv_inertia[i]=static_cast<float>(1.0/inertia);
    }
    auto copy=[&](float* dst,const std::vector<float>& src){if(!src.empty()){queue_.memcpy(dst,src.data(),src.size()*sizeof(float));transfer_stats_.record_host_to_device(src.size()*sizeof(float));}};
    copy(x_,x);copy(y_,y);copy(z_,z);copy(vx_,vx);copy(vy_,vy);copy(vz_,vz);copy(wx_,wx);copy(wy_,wy);copy(wz_,wz);copy(radius_,radius);copy(inverse_mass_,inv_mass);copy(inverse_inertia_,inv_inertia);
    const auto hist=count_*config_.history_slots_per_particle;
    if(hist){ queue_.memset(history_neighbor_,0xff,hist*sizeof(std::uint32_t));queue_.memset(history_seen_,0,hist*sizeof(std::uint32_t));queue_.memset(history_tx_,0,hist*sizeof(float));queue_.memset(history_ty_,0,hist*sizeof(float));queue_.memset(history_tz_,0,hist*sizeof(float)); }
    history_stamp_=1U;
    queue_.wait_and_throw();transfer_stats_.record_synchronization();
}

std::vector<ResidentDemParticle> ResidentDemSycl::download() const {
    std::vector<ResidentDemParticle> out(count_);
    std::vector<float> x(count_),y(count_),z(count_),vx(count_),vy(count_),vz(count_),wx(count_),wy(count_),wz(count_),r(count_),im(count_);
    auto copy=[&](std::vector<float>& dst,const float* src){if(!dst.empty()){queue_.memcpy(dst.data(),src,dst.size()*sizeof(float));transfer_stats_.record_device_to_host(dst.size()*sizeof(float));}};
    copy(x,x_);copy(y,y_);copy(z,z_);copy(vx,vx_);copy(vy,vy_);copy(vz,vz_);copy(wx,wx_);copy(wy,wy_);copy(wz,wz_);copy(r,radius_);copy(im,inverse_mass_);
    queue_.wait_and_throw();transfer_stats_.record_synchronization();
    for(std::size_t i=0;i<count_;++i){ out[i].position={x[i],y[i],z[i]};out[i].linear_velocity={vx[i],vy[i],vz[i]};out[i].angular_velocity={wx[i],wy[i],wz[i]};out[i].radius=r[i];out[i].mass=1.0/static_cast<double>(im[i]); }
    return out;
}

void ResidentDemSycl::build_neighbor_buckets() {
    if(count_==0U) return;
    queue_.memset(cell_counts_,0,cell_count_*sizeof(std::uint32_t));
    const float* x=x_;const float* y=y_;const float* z=z_;std::uint32_t* particle_cell=particle_cell_;std::uint32_t* counts=cell_counts_;
    const float x0=static_cast<float>(config_.domain_min.x),y0=static_cast<float>(config_.domain_min.y),z0=static_cast<float>(config_.domain_min.z),h=static_cast<float>(config_.cell_size);
    const std::size_t nx=nx_,ny=ny_,nz=nz_,n=count_;
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){
        const std::size_t i=id[0];
        auto ix=static_cast<long long>((x[i]-x0)/h);auto iy=static_cast<long long>((y[i]-y0)/h);auto iz=static_cast<long long>((z[i]-z0)/h);
        if(ix<0) ix=0;
        if(iy<0) iy=0;
        if(iz<0) iz=0;
        if(ix>=static_cast<long long>(nx)) ix=static_cast<long long>(nx)-1;
        if(iy>=static_cast<long long>(ny)) iy=static_cast<long long>(ny)-1;
        if(iz>=static_cast<long long>(nz)) iz=static_cast<long long>(nz)-1;
        const std::uint32_t cell=static_cast<std::uint32_t>((static_cast<std::size_t>(iz)*ny+static_cast<std::size_t>(iy))*nx+static_cast<std::size_t>(ix));
        particle_cell[i]=cell; GlobalAtomicU32(counts[cell]).fetch_add(1U);
    });
    std::uint32_t* offsets=cell_offsets_;std::uint32_t* cursor=cell_cursor_;const std::size_t cells=cell_count_;
    queue_.parallel_for(sycl::range<1>(1U),[=](sycl::id<1>){
        std::uint32_t sum=0U;
        for(std::size_t c=0;c<cells;++c){ offsets[c]=sum;cursor[c]=sum;sum+=counts[c]; }
        offsets[cells]=sum;
    });
    std::uint32_t* sorted=sorted_indices_;
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){
        const std::uint32_t i=static_cast<std::uint32_t>(id[0]);const std::uint32_t c=particle_cell[i];
        const std::uint32_t slot=GlobalAtomicU32(cursor[c]).fetch_add(1U);sorted[slot]=i;
    });
}

void ResidentDemSycl::clear_forces_with_gravity() {
    const std::size_t n=count_;float* fx=fx_;float* fy=fy_;float* fz=fz_;float* tx=tx_;float* ty=ty_;float* tz=tz_;const float* im=inverse_mass_;
    const float gx=static_cast<float>(config_.gravity.x),gy=static_cast<float>(config_.gravity.y),gz=static_cast<float>(config_.gravity.z);
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const std::size_t i=id[0];const float m=1.0F/im[i];fx[i]=m*gx;fy[i]=m*gy;fz[i]=m*gz;tx[i]=0.0F;ty[i]=0.0F;tz[i]=0.0F;});
}

void ResidentDemSycl::apply_contacts() {
    const std::size_t n=count_,nx=nx_,ny=ny_,nz=nz_,slots=config_.history_slots_per_particle;
    const float *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*r=radius_;
    float *fx=fx_,*fy=fy_,*fz=fz_,*tx=tx_,*ty=ty_,*tz=tz_;
    const std::uint32_t *pc=particle_cell_,*offsets=cell_offsets_,*sorted=sorted_indices_;
    std::uint32_t *hneighbor=history_neighbor_,*hseen=history_seen_;float *htx=history_tx_,*hty=history_ty_,*htz=history_tz_;
    const float kn=static_cast<float>(config_.normal_stiffness),cn=static_cast<float>(config_.normal_damping),kt0=static_cast<float>(config_.tangential_stiffness),ct0=static_cast<float>(config_.tangential_damping),mu=static_cast<float>(config_.friction),mur=static_cast<float>(config_.rolling_resistance),cohesion=static_cast<float>(config_.cohesion_force),dt=static_cast<float>(config_.dt);
    const std::uint32_t stamp=history_stamp_;
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){
        const std::size_t i=id[0];const std::size_t c=pc[i];const std::size_t cx=c%nx;const std::size_t yz=c/nx;const std::size_t cy=yz%ny;const std::size_t cz=yz/ny;
        for(int dz=-1;dz<=1;++dz){const long long zz=static_cast<long long>(cz)+dz;if(zz<0||zz>=static_cast<long long>(nz))continue;
            for(int dy=-1;dy<=1;++dy){const long long yy=static_cast<long long>(cy)+dy;if(yy<0||yy>=static_cast<long long>(ny))continue;
                for(int dx=-1;dx<=1;++dx){const long long xx=static_cast<long long>(cx)+dx;if(xx<0||xx>=static_cast<long long>(nx))continue;
                    const std::size_t nc=(static_cast<std::size_t>(zz)*ny+static_cast<std::size_t>(yy))*nx+static_cast<std::size_t>(xx);
                    for(std::uint32_t s=offsets[nc];s<offsets[nc+1U];++s){const std::size_t j=sorted[s];if(j<=i)continue;
                        const float dxp=x[i]-x[j],dyp=y[i]-y[j],dzp=z[i]-z[j];const float d2=dxp*dxp+dyp*dyp+dzp*dzp;const float rs=r[i]+r[j];if(!(d2<rs*rs))continue;
                        const float d=sycl::sqrt(d2>1.0e-20F?d2:1.0e-20F);const float nxv=dxp/d,nyv=dyp/d,nzv=dzp/d;const float pen=rs-d;const float root=sycl::sqrt(pen>0.0F?pen:0.0F);
                        const float cpx=x[i]-nxv*r[i],cpy=y[i]-nyv*r[i],cpz=z[i]-nzv*r[i];
                        const float rix=cpx-x[i],riy=cpy-y[i],riz=cpz-z[i],rjx=cpx-x[j],rjy=cpy-y[j],rjz=cpz-z[j];
                        const float vix=vx[i]+wy[i]*riz-wz[i]*riy,viy=vy[i]+wz[i]*rix-wx[i]*riz,viz=vz[i]+wx[i]*riy-wy[i]*rix;
                        const float vjx=vx[j]+wy[j]*rjz-wz[j]*rjy,vjy=vy[j]+wz[j]*rjx-wx[j]*rjz,vjz=vz[j]+wx[j]*rjy-wy[j]*rjx;
                        const float rvx=vix-vjx,rvy=viy-vjy,rvz=viz-vjz;const float vn=rvx*nxv+rvy*nyv+rvz*nzv;
                        float fn=kn*pen*root-cn*root*vn-cohesion;if(fn < -cohesion)fn=-cohesion;
                        const float vtx=rvx-vn*nxv,vty=rvy-vn*nyv,vtz=rvz-vn*nzv;
                        std::size_t hslot=slots;std::size_t empty=slots;const std::size_t base=i*slots;
                        for(std::size_t k=0;k<slots;++k){const auto nb=hneighbor[base+k];if(nb==static_cast<std::uint32_t>(j)){hslot=k;break;}if(nb==std::numeric_limits<std::uint32_t>::max()&&empty==slots)empty=k;}
                        if(hslot==slots && empty<slots){hslot=empty;hneighbor[base+hslot]=static_cast<std::uint32_t>(j);htx[base+hslot]=0.0F;hty[base+hslot]=0.0F;htz[base+hslot]=0.0F;}
                        float hx=0.0F,hy=0.0F,hz=0.0F;
                        if(hslot<slots){const std::size_t hs=base+hslot;hx=htx[hs]+vtx*dt;hy=hty[hs]+vty*dt;hz=htz[hs]+vtz*dt;const float hn=hx*nxv+hy*nyv+hz*nzv;hx-=hn*nxv;hy-=hn*nyv;hz-=hn*nzv;hseen[hs]=stamp;}
                        const float kt=kt0*root,ct=ct0*root;float ftx=-kt*hx-ct*vtx,fty=-kt*hy-ct*vty,ftz=-kt*hz-ct*vtz;
                        const float ftm=sycl::sqrt(ftx*ftx+fty*fty+ftz*ftz);const float maxft=mu*(fn<0.0F?-fn:fn);
                        if(ftm>maxft && ftm>1.0e-20F){const float scale=maxft/ftm;ftx*=scale;fty*=scale;ftz*=scale;if(hslot<slots&&kt>1.0e-20F){hx=-(ftx+ct*vtx)/kt;hy=-(fty+ct*vty)/kt;hz=-(ftz+ct*vtz)/kt;}}
                        if(hslot<slots){const std::size_t hs=base+hslot;htx[hs]=hx;hty[hs]=hy;htz[hs]=hz;}
                        const float fpx=nxv*fn+ftx,fpy=nyv*fn+fty,fpz=nzv*fn+ftz;
                        GlobalAtomicFloat(fx[i]).fetch_add(fpx);GlobalAtomicFloat(fy[i]).fetch_add(fpy);GlobalAtomicFloat(fz[i]).fetch_add(fpz);GlobalAtomicFloat(fx[j]).fetch_add(-fpx);GlobalAtomicFloat(fy[j]).fetch_add(-fpy);GlobalAtomicFloat(fz[j]).fetch_add(-fpz);
                        const float tix=riy*fpz-riz*fpy,tiy=riz*fpx-rix*fpz,tiz=rix*fpy-riy*fpx;const float tjx=-(rjy*fpz-rjz*fpy),tjy=-(rjz*fpx-rjx*fpz),tjz=-(rjx*fpy-rjy*fpx);
                        GlobalAtomicFloat(tx[i]).fetch_add(tix);GlobalAtomicFloat(ty[i]).fetch_add(tiy);GlobalAtomicFloat(tz[i]).fetch_add(tiz);GlobalAtomicFloat(tx[j]).fetch_add(tjx);GlobalAtomicFloat(ty[j]).fetch_add(tjy);GlobalAtomicFloat(tz[j]).fetch_add(tjz);
                        const float wrx=wx[i]-wx[j],wry=wy[i]-wy[j],wrz=wz[i]-wz[j],wm=sycl::sqrt(wrx*wrx+wry*wry+wrz*wrz);
                        if(wm>1.0e-20F&&mur>0.0F){const float reff=r[i]*r[j]/(r[i]+r[j]);const float trscale=-mur*(fn<0.0F?-fn:fn)*reff/wm;const float trx=wrx*trscale,tryv=wry*trscale,trz=wrz*trscale;GlobalAtomicFloat(tx[i]).fetch_add(trx);GlobalAtomicFloat(ty[i]).fetch_add(tryv);GlobalAtomicFloat(tz[i]).fetch_add(trz);GlobalAtomicFloat(tx[j]).fetch_add(-trx);GlobalAtomicFloat(ty[j]).fetch_add(-tryv);GlobalAtomicFloat(tz[j]).fetch_add(-trz);}
                    }
                }
            }
        }
    });
}

void ResidentDemSycl::cleanup_history(std::uint32_t stamp) {
    const std::size_t total=count_*config_.history_slots_per_particle;std::uint32_t* nb=history_neighbor_;std::uint32_t* seen=history_seen_;float* hx=history_tx_;float* hy=history_ty_;float* hz=history_tz_;
    queue_.parallel_for(sycl::range<1>(total),[=](sycl::id<1> id){const std::size_t i=id[0];if(nb[i]!=std::numeric_limits<std::uint32_t>::max()&&seen[i]!=stamp){nb[i]=std::numeric_limits<std::uint32_t>::max();hx[i]=0.0F;hy[i]=0.0F;hz[i]=0.0F;}});
}

void ResidentDemSycl::integrate() {
    const std::size_t n=count_;float *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_;const float *im=inverse_mass_,*ii=inverse_inertia_,*fx=fx_,*fy=fy_,*fz=fz_,*tx=tx_,*ty=ty_,*tz=tz_;const float dt=static_cast<float>(config_.dt);
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const std::size_t i=id[0];vx[i]+=fx[i]*im[i]*dt;vy[i]+=fy[i]*im[i]*dt;vz[i]+=fz[i]*im[i]*dt;wx[i]+=tx[i]*ii[i]*dt;wy[i]+=ty[i]*ii[i]*dt;wz[i]+=tz[i]*ii[i]*dt;x[i]+=vx[i]*dt;y[i]+=vy[i]*dt;z[i]+=vz[i]*dt;});
}

void ResidentDemSycl::step(std::size_t steps) {
    for(std::size_t s=0;s<steps;++s){
        build_neighbor_buckets();clear_forces_with_gravity();
        ++history_stamp_;if(history_stamp_==0U)history_stamp_=1U;
        apply_contacts();cleanup_history(history_stamp_);integrate();
    }
}

double ResidentDemSycl::total_kinetic_energy() const {
    if(!reduction_scalar_){reduction_scalar_=sycl::malloc_shared<float>(1U,queue_);if(!reduction_scalar_)throw std::bad_alloc{};}
    *reduction_scalar_=0.0F;const std::size_t n=count_;const float *vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*im=inverse_mass_,*ii=inverse_inertia_;
    auto reduction=sycl::reduction(reduction_scalar_,sycl::plus<float>());
    queue_.parallel_for(sycl::range<1>(n),reduction,[=](sycl::id<1> id,auto& sum){const std::size_t i=id[0];const float m=1.0F/im[i],inertia=1.0F/ii[i];sum.combine(0.5F*m*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])+0.5F*inertia*(wx[i]*wx[i]+wy[i]*wy[i]+wz[i]*wz[i]));}).wait_and_throw();
    transfer_stats_.record_synchronization();return static_cast<double>(*reduction_scalar_);
}

std::size_t ResidentDemSycl::resident_bytes() const noexcept {
    const std::size_t cap=config_.capacity,hist=cap*config_.history_slots_per_particle;
    return cap*(18U*sizeof(float)+2U*sizeof(std::uint32_t))+(cell_count_*2U+cell_count_+1U)*sizeof(std::uint32_t)+hist*(3U*sizeof(float)+2U*sizeof(std::uint32_t))+(reduction_scalar_?sizeof(float):0U);
}
ResidentDemDeviceView ResidentDemSycl::device_view() noexcept { return {count_,x_,y_,z_,vx_,vy_,vz_,wx_,wy_,wz_,radius_,inverse_mass_}; }

} // namespace cfd::multibody
#endif
