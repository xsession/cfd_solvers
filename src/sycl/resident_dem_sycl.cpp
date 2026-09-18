#include "cfd/multibody/resident_dem_sycl.hpp"

#if defined(CFD_HAS_SYCL)

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

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
    auto allocu64=[&](std::size_t n){ auto* p=sycl::malloc_device<std::uint64_t>(n,queue_); if(!p)throw std::bad_alloc{}; return p; };
    auto alloci32=[&](std::size_t n){ auto* p=sycl::malloc_device<std::int32_t>(n,queue_); if(!p)throw std::bad_alloc{}; return p; };
    try {
        global_id_=allocu64(cap);owner_rank_=alloci32(cap);
        x_=allocf(cap);y_=allocf(cap);z_=allocf(cap);vx_=allocf(cap);vy_=allocf(cap);vz_=allocf(cap);
        wx_=allocf(cap);wy_=allocf(cap);wz_=allocf(cap);radius_=allocf(cap);inverse_mass_=allocf(cap);inverse_inertia_=allocf(cap);
        fx_=allocf(cap);fy_=allocf(cap);fz_=allocf(cap);tx_=allocf(cap);ty_=allocf(cap);tz_=allocf(cap);
        particle_cell_=allocu(cap);sorted_indices_=allocu(cap);
        cell_counts_=allocu(cell_count_);cell_offsets_=allocu(cell_count_+1U);cell_cursor_=allocu(cell_count_);
        history_neighbor_=allocu64(hist);history_seen_=allocu(hist);history_tx_=allocf(hist);history_ty_=allocf(hist);history_tz_=allocf(hist);
        queue_.memset(history_neighbor_,0xff,hist*sizeof(std::uint64_t));
        queue_.memset(history_seen_,0,hist*sizeof(std::uint32_t));
        queue_.memset(history_tx_,0,hist*sizeof(float));queue_.memset(history_ty_,0,hist*sizeof(float));queue_.memset(history_tz_,0,hist*sizeof(float));
    } catch(...) { release(); throw; }
}

ResidentDemSycl::~ResidentDemSycl() noexcept { release(); }
void ResidentDemSycl::release() noexcept {
    try { queue_.wait_and_throw(); } catch(...) {}
    auto freep=[&](auto*& p){ if(p){ sycl::free(p,queue_);p=nullptr; } };
    freep(global_id_);freep(owner_rank_);freep(x_);freep(y_);freep(z_);freep(vx_);freep(vy_);freep(vz_);freep(wx_);freep(wy_);freep(wz_);
    freep(radius_);freep(inverse_mass_);freep(inverse_inertia_);freep(fx_);freep(fy_);freep(fz_);freep(tx_);freep(ty_);freep(tz_);
    freep(particle_cell_);freep(cell_counts_);freep(cell_offsets_);freep(cell_cursor_);freep(sorted_indices_);
    freep(history_neighbor_);freep(history_seen_);freep(history_tx_);freep(history_ty_);freep(history_tz_);freep(reduction_scalar_);freep(compaction_count_);
}

void ResidentDemSycl::upload(std::span<const ResidentDemParticle> particles) {
    if(particles.size()>config_.capacity) throw std::invalid_argument("resident DEM particle capacity exceeded");
    count_=particles.size();
    std::vector<std::uint64_t> gids(count_);
    std::vector<std::int32_t> owners(count_);
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
        gids[i]=p.global_id==invalid_dem_global_id?static_cast<std::uint64_t>(i):p.global_id;
        owners[i]=static_cast<std::int32_t>(p.owner_rank);
        x[i]=static_cast<float>(p.position.x);y[i]=static_cast<float>(p.position.y);z[i]=static_cast<float>(p.position.z);
        vx[i]=static_cast<float>(p.linear_velocity.x);vy[i]=static_cast<float>(p.linear_velocity.y);vz[i]=static_cast<float>(p.linear_velocity.z);
        wx[i]=static_cast<float>(p.angular_velocity.x);wy[i]=static_cast<float>(p.angular_velocity.y);wz[i]=static_cast<float>(p.angular_velocity.z);
        radius[i]=static_cast<float>(p.radius);inv_mass[i]=static_cast<float>(1.0/p.mass);
        const double inertia=0.4*p.mass*p.radius*p.radius;inv_inertia[i]=static_cast<float>(1.0/inertia);
    }
    auto sorted_ids=gids;
    std::sort(sorted_ids.begin(),sorted_ids.end());
    if(std::adjacent_find(sorted_ids.begin(),sorted_ids.end())!=sorted_ids.end()) throw std::invalid_argument("resident DEM global IDs must be unique");
    auto copyf=[&](float* dst,const std::vector<float>& src){if(!src.empty()){queue_.memcpy(dst,src.data(),src.size()*sizeof(float));transfer_stats_.record_host_to_device(src.size()*sizeof(float));}};
    if(!gids.empty()){queue_.memcpy(global_id_,gids.data(),gids.size()*sizeof(std::uint64_t));transfer_stats_.record_host_to_device(gids.size()*sizeof(std::uint64_t));}
    if(!owners.empty()){queue_.memcpy(owner_rank_,owners.data(),owners.size()*sizeof(std::int32_t));transfer_stats_.record_host_to_device(owners.size()*sizeof(std::int32_t));}
    copyf(x_,x);copyf(y_,y);copyf(z_,z);copyf(vx_,vx);copyf(vy_,vy);copyf(vz_,vz);copyf(wx_,wx);copyf(wy_,wy);copyf(wz_,wz);copyf(radius_,radius);copyf(inverse_mass_,inv_mass);copyf(inverse_inertia_,inv_inertia);
    const auto hist=count_*config_.history_slots_per_particle;
    if(hist){ queue_.memset(history_neighbor_,0xff,hist*sizeof(std::uint64_t));queue_.memset(history_seen_,0,hist*sizeof(std::uint32_t));queue_.memset(history_tx_,0,hist*sizeof(float));queue_.memset(history_ty_,0,hist*sizeof(float));queue_.memset(history_tz_,0,hist*sizeof(float)); }
    history_stamp_=1U;
    queue_.wait_and_throw();transfer_stats_.record_synchronization();
}

std::vector<ResidentDemParticle> ResidentDemSycl::download() const {
    std::vector<ResidentDemParticle> out(count_);
    std::vector<std::uint64_t> gids(count_);
    std::vector<std::int32_t> owners(count_);
    std::vector<float> x(count_),y(count_),z(count_),vx(count_),vy(count_),vz(count_),wx(count_),wy(count_),wz(count_),r(count_),im(count_);
    auto copyf=[&](std::vector<float>& dst,const float* src){if(!dst.empty()){queue_.memcpy(dst.data(),src,dst.size()*sizeof(float));transfer_stats_.record_device_to_host(dst.size()*sizeof(float));}};
    if(!gids.empty()){queue_.memcpy(gids.data(),global_id_,gids.size()*sizeof(std::uint64_t));transfer_stats_.record_device_to_host(gids.size()*sizeof(std::uint64_t));}
    if(!owners.empty()){queue_.memcpy(owners.data(),owner_rank_,owners.size()*sizeof(std::int32_t));transfer_stats_.record_device_to_host(owners.size()*sizeof(std::int32_t));}
    copyf(x,x_);copyf(y,y_);copyf(z,z_);copyf(vx,vx_);copyf(vy,vy_);copyf(vz,vz_);copyf(wx,wx_);copyf(wy,wy_);copyf(wz,wz_);copyf(r,radius_);copyf(im,inverse_mass_);
    queue_.wait_and_throw();transfer_stats_.record_synchronization();
    for(std::size_t i=0;i<count_;++i){ out[i].global_id=gids[i];out[i].owner_rank=static_cast<int>(owners[i]);out[i].position={x[i],y[i],z[i]};out[i].linear_velocity={vx[i],vy[i],vz[i]};out[i].angular_velocity={wx[i],wy[i],wz[i]};out[i].radius=r[i];out[i].mass=1.0/static_cast<double>(im[i]); }
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

void ResidentDemSycl::apply_contacts() { apply_contacts_partition(nullptr,0); }

void ResidentDemSycl::apply_contacts_partition(const std::uint8_t* boundary_mask,int partition_mode) {
    const std::size_t n=count_,nx=nx_,ny=ny_,nz=nz_,slots=config_.history_slots_per_particle;
    const std::uint64_t* gid=global_id_;
    const float *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*r=radius_;
    float *fx=fx_,*fy=fy_,*fz=fz_,*tx=tx_,*ty=ty_,*tz=tz_;
    const std::uint32_t *pc=particle_cell_,*offsets=cell_offsets_,*sorted=sorted_indices_;
    std::uint64_t *hneighbor=history_neighbor_;std::uint32_t* hseen=history_seen_;float *htx=history_tx_,*hty=history_ty_,*htz=history_tz_;
    const float kn=static_cast<float>(config_.normal_stiffness),cn=static_cast<float>(config_.normal_damping),kt0=static_cast<float>(config_.tangential_stiffness),ct0=static_cast<float>(config_.tangential_damping),mu=static_cast<float>(config_.friction),mur=static_cast<float>(config_.rolling_resistance),cohesion=static_cast<float>(config_.cohesion_force),dt=static_cast<float>(config_.dt);
    const std::uint32_t stamp=history_stamp_;
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){
        const std::size_t i=id[0];const std::size_t c=pc[i];const std::size_t cx=c%nx;const std::size_t yz=c/nx;const std::size_t cy=yz%ny;const std::size_t cz=yz/ny;
        for(int dz=-1;dz<=1;++dz){const long long zz=static_cast<long long>(cz)+dz;if(zz<0||zz>=static_cast<long long>(nz))continue;
            for(int dy=-1;dy<=1;++dy){const long long yy=static_cast<long long>(cy)+dy;if(yy<0||yy>=static_cast<long long>(ny))continue;
                for(int dx=-1;dx<=1;++dx){const long long xx=static_cast<long long>(cx)+dx;if(xx<0||xx>=static_cast<long long>(nx))continue;
                    const std::size_t nc=(static_cast<std::size_t>(zz)*ny+static_cast<std::size_t>(yy))*nx+static_cast<std::size_t>(xx);
                    for(std::uint32_t s=offsets[nc];s<offsets[nc+1U];++s){const std::size_t j=sorted[s];if(gid[j]<=gid[i])continue;
                        const bool boundary_pair=boundary_mask && (boundary_mask[i]!=0U || boundary_mask[j]!=0U);
                        if(partition_mode==1 && boundary_pair)continue;
                        if(partition_mode==2 && !boundary_pair)continue;
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
                        for(std::size_t k=0;k<slots;++k){const auto nb=hneighbor[base+k];if(nb==gid[j]){hslot=k;break;}if(nb==invalid_dem_global_id&&empty==slots)empty=k;}
                        if(hslot==slots && empty<slots){hslot=empty;hneighbor[base+hslot]=gid[j];htx[base+hslot]=0.0F;hty[base+hslot]=0.0F;htz[base+hslot]=0.0F;}
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
    const std::size_t total=count_*config_.history_slots_per_particle;std::uint64_t* nb=history_neighbor_;std::uint32_t* seen=history_seen_;float* hx=history_tx_;float* hy=history_ty_;float* hz=history_tz_;
    queue_.parallel_for(sycl::range<1>(total),[=](sycl::id<1> id){const std::size_t i=id[0];if(nb[i]!=invalid_dem_global_id&&seen[i]!=stamp){nb[i]=invalid_dem_global_id;hx[i]=0.0F;hy[i]=0.0F;hz[i]=0.0F;}});
}

void ResidentDemSycl::integrate() {
    const std::size_t n=count_;float *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_;const float *im=inverse_mass_,*ii=inverse_inertia_,*fx=fx_,*fy=fy_,*fz=fz_,*tx=tx_,*ty=ty_,*tz=tz_;const float dt=static_cast<float>(config_.dt);
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){const std::size_t i=id[0];vx[i]+=fx[i]*im[i]*dt;vy[i]+=fy[i]*im[i]*dt;vz[i]+=fz[i]*im[i]*dt;wx[i]+=tx[i]*ii[i]*dt;wy[i]+=ty[i]*ii[i]*dt;wz[i]+=tz[i]*ii[i]*dt;x[i]+=vx[i]*dt;y[i]+=vy[i]*dt;z[i]+=vz[i]*dt;});
}

void ResidentDemSycl::step(std::size_t steps) {
    for(std::size_t s=0;s<steps;++s){
        begin_distributed_contact_step();
        apply_contacts();
        finish_distributed_contact_step();
    }
}

void ResidentDemSycl::begin_distributed_contact_step() {
    build_neighbor_buckets();
    clear_forces_with_gravity();
    ++history_stamp_;
    if(history_stamp_==0U) history_stamp_=1U;
}

void ResidentDemSycl::apply_owned_contacts_interior_device(const std::uint8_t* boundary_mask) {
    if(count_!=0U && !boundary_mask) throw std::invalid_argument("resident DEM interior partition requires a boundary mask");
    apply_contacts_partition(boundary_mask,1);
}

void ResidentDemSycl::apply_owned_contacts_boundary_device(const std::uint8_t* boundary_mask) {
    if(count_!=0U && !boundary_mask) throw std::invalid_argument("resident DEM boundary partition requires a boundary mask");
    apply_contacts_partition(boundary_mask,2);
}

void ResidentDemSycl::apply_ghost_contacts_device(
    const ResidentDemPackedParticle* device_ghosts,
    std::size_t ghost_count,
    const std::uint8_t* boundary_mask,
    int local_rank,
    ResidentDemPackedForce* device_remote_forces,
    std::uint32_t* shared_remote_force_count,
    std::size_t remote_force_capacity) {
    if(ghost_count==0U || count_==0U) return;
    if(!device_ghosts || !boundary_mask || !device_remote_forces || !shared_remote_force_count || remote_force_capacity==0U) {
        throw std::invalid_argument("resident DEM ghost contact requires valid device buffers");
    }
    const std::size_t n=count_,slots=config_.history_slots_per_particle;
    const auto* ghosts=device_ghosts;
    const std::uint64_t* gid=global_id_;
    const std::int32_t* owner=owner_rank_;
    const float *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*r=radius_;
    float *fx=fx_,*fy=fy_,*fz=fz_,*tx=tx_,*ty=ty_,*tz=tz_;
    std::uint64_t *hneighbor=history_neighbor_;std::uint32_t* hseen=history_seen_;float *htx=history_tx_,*hty=history_ty_,*htz=history_tz_;
    const float kn=static_cast<float>(config_.normal_stiffness),cn=static_cast<float>(config_.normal_damping),kt0=static_cast<float>(config_.tangential_stiffness),ct0=static_cast<float>(config_.tangential_damping),mu=static_cast<float>(config_.friction),mur=static_cast<float>(config_.rolling_resistance),cohesion=static_cast<float>(config_.cohesion_force),dt=static_cast<float>(config_.dt);
    const std::uint32_t stamp=history_stamp_;
    auto* remote=device_remote_forces;auto* remote_count=shared_remote_force_count;const auto remote_cap=remote_force_capacity;
    const auto local_owner=static_cast<std::int32_t>(local_rank);
    queue_.parallel_for(sycl::range<1>(n),[=](sycl::id<1> id){
        const std::size_t i=id[0];
        if(boundary_mask[i]==0U || owner[i]!=local_owner) return;
        for(std::size_t g=0;g<ghost_count;++g){
            const auto ghost=ghosts[g];
            if(ghost.global_id==invalid_dem_global_id || gid[i]>=ghost.global_id) continue;
            const float dxp=x[i]-ghost.x,dyp=y[i]-ghost.y,dzp=z[i]-ghost.z;
            const float d2=dxp*dxp+dyp*dyp+dzp*dzp,rs=r[i]+ghost.radius;
            if(!(d2<rs*rs)) continue;
            const float d=sycl::sqrt(d2>1.0e-20F?d2:1.0e-20F),nxv=dxp/d,nyv=dyp/d,nzv=dzp/d;
            const float pen=rs-d,root=sycl::sqrt(pen>0.0F?pen:0.0F);
            const float cpx=x[i]-nxv*r[i],cpy=y[i]-nyv*r[i],cpz=z[i]-nzv*r[i];
            const float rix=cpx-x[i],riy=cpy-y[i],riz=cpz-z[i],rjx=cpx-ghost.x,rjy=cpy-ghost.y,rjz=cpz-ghost.z;
            const float vix=vx[i]+wy[i]*riz-wz[i]*riy,viy=vy[i]+wz[i]*rix-wx[i]*riz,viz=vz[i]+wx[i]*riy-wy[i]*rix;
            const float vjx=ghost.vx+ghost.wy*rjz-ghost.wz*rjy,vjy=ghost.vy+ghost.wz*rjx-ghost.wx*rjz,vjz=ghost.vz+ghost.wx*rjy-ghost.wy*rjx;
            const float rvx=vix-vjx,rvy=viy-vjy,rvz=viz-vjz,vn=rvx*nxv+rvy*nyv+rvz*nzv;
            float fn=kn*pen*root-cn*root*vn-cohesion;if(fn < -cohesion) fn=-cohesion;
            const float vtx=rvx-vn*nxv,vty=rvy-vn*nyv,vtz=rvz-vn*nzv;
            std::size_t hslot=slots,empty=slots;const std::size_t base=i*slots;
            for(std::size_t k=0;k<slots;++k){const auto nb=hneighbor[base+k];if(nb==ghost.global_id){hslot=k;break;}if(nb==invalid_dem_global_id&&empty==slots)empty=k;}
            if(hslot==slots&&empty<slots){hslot=empty;hneighbor[base+hslot]=ghost.global_id;htx[base+hslot]=0.0F;hty[base+hslot]=0.0F;htz[base+hslot]=0.0F;}
            float hx=0.0F,hy=0.0F,hz=0.0F;
            if(hslot<slots){const std::size_t hs=base+hslot;hx=htx[hs]+vtx*dt;hy=hty[hs]+vty*dt;hz=htz[hs]+vtz*dt;const float hn=hx*nxv+hy*nyv+hz*nzv;hx-=hn*nxv;hy-=hn*nyv;hz-=hn*nzv;hseen[hs]=stamp;}
            const float kt=kt0*root,ct=ct0*root;float ftx=-kt*hx-ct*vtx,fty=-kt*hy-ct*vty,ftz=-kt*hz-ct*vtz;
            const float ftm=sycl::sqrt(ftx*ftx+fty*fty+ftz*ftz),maxft=mu*(fn<0.0F?-fn:fn);
            if(ftm>maxft&&ftm>1.0e-20F){const float scale=maxft/ftm;ftx*=scale;fty*=scale;ftz*=scale;if(hslot<slots&&kt>1.0e-20F){hx=-(ftx+ct*vtx)/kt;hy=-(fty+ct*vty)/kt;hz=-(ftz+ct*vtz)/kt;}}
            if(hslot<slots){const auto hs=base+hslot;htx[hs]=hx;hty[hs]=hy;htz[hs]=hz;}
            const float fpx=nxv*fn+ftx,fpy=nyv*fn+fty,fpz=nzv*fn+ftz;
            const float tix=riy*fpz-riz*fpy,tiy=riz*fpx-rix*fpz,tiz=rix*fpy-riy*fpx;
            float tjx=-(rjy*fpz-rjz*fpy),tjy=-(rjz*fpx-rjx*fpz),tjz=-(rjx*fpy-rjy*fpx);
            float trxi=0.0F,tryi=0.0F,trzi=0.0F;
            const float wrx=wx[i]-ghost.wx,wry=wy[i]-ghost.wy,wrz=wz[i]-ghost.wz,wm=sycl::sqrt(wrx*wrx+wry*wry+wrz*wrz);
            if(wm>1.0e-20F&&mur>0.0F){const float reff=r[i]*ghost.radius/(r[i]+ghost.radius),trscale=-mur*(fn<0.0F?-fn:fn)*reff/wm;trxi=wrx*trscale;tryi=wry*trscale;trzi=wrz*trscale;tjx-=trxi;tjy-=tryi;tjz-=trzi;}
            fx[i]+=fpx;fy[i]+=fpy;fz[i]+=fpz;tx[i]+=tix+trxi;ty[i]+=tiy+tryi;tz[i]+=tiz+trzi;
            const auto slot=GlobalAtomicU32(*remote_count).fetch_add(1U);
            if(static_cast<std::size_t>(slot)<remote_cap){auto& out=remote[slot];out.global_id=ghost.global_id;out.owner_rank=ghost.owner_rank;out.fx=-fpx;out.fy=-fpy;out.fz=-fpz;out.tx=tjx;out.ty=tjy;out.tz=tjz;}
        }
    });
}

void ResidentDemSycl::apply_bonds_device(
    ResidentDemPackedBond* device_bonds,
    std::size_t bond_count,
    const ResidentDemPackedParticle* device_ghosts,
    std::size_t ghost_count,
    int local_rank,
    ResidentDemPackedForce* device_remote_forces,
    std::uint32_t* shared_remote_force_count,
    std::size_t remote_force_capacity) {
    if(bond_count==0U || count_==0U) return;
    if(!device_bonds || !device_remote_forces || !shared_remote_force_count || remote_force_capacity==0U) throw std::invalid_argument("resident DEM bond kernel requires valid buffers");
    const auto n=count_;const auto local_owner=static_cast<std::int32_t>(local_rank);const float dt=static_cast<float>(config_.dt);
    const auto* gid=global_id_;const auto* owner=owner_rank_;const auto *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_;
    float *fx=fx_,*fy=fy_,*fz=fz_;auto* bonds=device_bonds;const auto* ghosts=device_ghosts;auto* remote=device_remote_forces;auto* remote_count=shared_remote_force_count;const auto remote_cap=remote_force_capacity;
    queue_.parallel_for(sycl::range<1>(bond_count),[=](sycl::id<1> id){
        auto& bond=bonds[id[0]];if(bond.broken!=0U || bond.particle_a==invalid_dem_global_id || !(bond.particle_a<bond.particle_b)) return;
        std::size_t ia=n;for(std::size_t i=0;i<n;++i)if(gid[i]==bond.particle_a){ia=i;break;}if(ia==n||owner[ia]!=local_owner)return;
        bool b_remote=false;std::size_t ib=n,gb=ghost_count;
        for(std::size_t i=0;i<n;++i)if(gid[i]==bond.particle_b){ib=i;break;}
        if(ib==n&&ghosts){for(std::size_t g=0;g<ghost_count;++g)if(ghosts[g].global_id==bond.particle_b){gb=g;b_remote=true;break;}}
        if(ib==n&&!b_remote)return;
        const float bx=b_remote?ghosts[gb].x:x[ib],by=b_remote?ghosts[gb].y:y[ib],bz=b_remote?ghosts[gb].z:z[ib];
        const float bvx=b_remote?ghosts[gb].vx:vx[ib],bvy=b_remote?ghosts[gb].vy:vy[ib],bvz=b_remote?ghosts[gb].vz:vz[ib];
        const float dx=bx-x[ia],dy=by-y[ia],dz=bz-z[ia],d2=dx*dx+dy*dy+dz*dz,d=sycl::sqrt(d2>1.0e-20F?d2:1.0e-20F),nxv=dx/d,nyv=dy/d,nzv=dz/d;
        const float rvx=bvx-vx[ia],rvy=bvy-vy[ia],rvz=bvz-vz[ia],vn=rvx*nxv+rvy*nyv+rvz*nzv,vtx=rvx-vn*nxv,vty=rvy-vn*nyv,vtz=rvz-vn*nzv;
        float hx=bond.tangential_x+vtx*dt,hy=bond.tangential_y+vty*dt,hz=bond.tangential_z+vtz*dt;const float hn=hx*nxv+hy*nyv+hz*nzv;hx-=hn*nxv;hy-=hn*nyv;hz-=hn*nzv;
        const float extension=d-bond.rest_length,normal_elastic=bond.normal_stiffness*extension;
        const float sx=bond.shear_stiffness*hx,syv=bond.shear_stiffness*hy,sz=bond.shear_stiffness*hz,shear_mag=sycl::sqrt(sx*sx+syv*syv+sz*sz);
        const float tensile_ratio=normal_elastic>0.0F?normal_elastic/bond.tensile_failure_force:0.0F,shear_ratio=shear_mag/bond.shear_failure_force,failure=tensile_ratio>shear_ratio?tensile_ratio:shear_ratio;
        if(failure>bond.damage_onset_ratio){const float denom=1.0F-bond.damage_onset_ratio;float target=denom>1.0e-12F?(failure-bond.damage_onset_ratio)/denom:1.0F;if(target<0.0F)target=0.0F;if(target>1.0F)target=1.0F;if(target>bond.damage)bond.damage=target;}
        bond.tangential_x=hx;bond.tangential_y=hy;bond.tangential_z=hz;
        if(failure>=1.0F||bond.damage>=1.0F){bond.damage=1.0F;bond.broken=1U;return;}
        const float intact=1.0F-bond.damage,normal_force=intact*normal_elastic+bond.normal_damping*vn;
        const float fsx=intact*sx+bond.shear_damping*vtx,fsy=intact*syv+bond.shear_damping*vty,fsz=intact*sz+bond.shear_damping*vtz;
        const float fax=nxv*normal_force+fsx,fay=nyv*normal_force+fsy,faz=nzv*normal_force+fsz;
        GlobalAtomicFloat(fx[ia]).fetch_add(fax);GlobalAtomicFloat(fy[ia]).fetch_add(fay);GlobalAtomicFloat(fz[ia]).fetch_add(faz);
        if(!b_remote){GlobalAtomicFloat(fx[ib]).fetch_add(-fax);GlobalAtomicFloat(fy[ib]).fetch_add(-fay);GlobalAtomicFloat(fz[ib]).fetch_add(-faz);return;}
        const auto slot=GlobalAtomicU32(*remote_count).fetch_add(1U);if(static_cast<std::size_t>(slot)<remote_cap){auto& out=remote[slot];out.global_id=bond.particle_b;out.owner_rank=ghosts[gb].owner_rank;out.fx=-fax;out.fy=-fay;out.fz=-faz;out.tx=0.0F;out.ty=0.0F;out.tz=0.0F;}
    });
}

void ResidentDemSycl::apply_remote_forces_device(const ResidentDemPackedForce* device_forces,std::size_t force_count) {
    if(force_count==0U)return;
    if(!device_forces)throw std::invalid_argument("resident DEM remote force buffer is null");
    const auto n=count_;const auto* gid=global_id_;float *fx=fx_,*fy=fy_,*fz=fz_,*tx=tx_,*ty=ty_,*tz=tz_;const auto* forces=device_forces;
    queue_.parallel_for(sycl::range<1>(force_count),[=](sycl::id<1> id){const auto f=forces[id[0]];for(std::size_t i=0;i<n;++i){if(gid[i]!=f.global_id)continue;GlobalAtomicFloat(fx[i]).fetch_add(f.fx);GlobalAtomicFloat(fy[i]).fetch_add(f.fy);GlobalAtomicFloat(fz[i]).fetch_add(f.fz);GlobalAtomicFloat(tx[i]).fetch_add(f.tx);GlobalAtomicFloat(ty[i]).fetch_add(f.ty);GlobalAtomicFloat(tz[i]).fetch_add(f.tz);break;}});
}

void ResidentDemSycl::finish_distributed_contact_step() {
    cleanup_history(history_stamp_);
    integrate();
}

double ResidentDemSycl::total_kinetic_energy() const {
    if(!reduction_scalar_){reduction_scalar_=sycl::malloc_shared<float>(1U,queue_);if(!reduction_scalar_)throw std::bad_alloc{};}
    *reduction_scalar_=0.0F;const std::size_t n=count_;const float *vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*im=inverse_mass_,*ii=inverse_inertia_;
    auto reduction=sycl::reduction(reduction_scalar_,sycl::plus<float>());
    queue_.parallel_for(sycl::range<1>(n),reduction,[=](sycl::id<1> id,auto& sum){const std::size_t i=id[0];const float m=1.0F/im[i],inertia=1.0F/ii[i];sum.combine(0.5F*m*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])+0.5F*inertia*(wx[i]*wx[i]+wy[i]*wy[i]+wz[i]*wz[i]));}).wait_and_throw();
    transfer_stats_.record_synchronization();return static_cast<double>(*reduction_scalar_);
}

void ResidentDemSycl::pack_indices_device(
    const std::uint32_t* device_indices,
    std::size_t index_count,
    ResidentDemPackedParticle* device_particles,
    ResidentDemPackedHistory* device_history) const {
    if(index_count==0U) return;
    if(!device_indices || !device_particles || !device_history) throw std::invalid_argument("resident DEM pack requires non-null device pointers");
    const auto slots=config_.history_slots_per_particle;
    const auto n=count_;
    const auto* gid=global_id_;const auto* owner=owner_rank_;
    const auto *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*radius=radius_,*im=inverse_mass_,*ii=inverse_inertia_;
    const auto* hneighbor=history_neighbor_;const auto* htx=history_tx_;const auto* hty=history_ty_;const auto* htz=history_tz_;
    queue_.parallel_for(sycl::range<1>(index_count),[=](sycl::id<1> id){
        const std::size_t p=id[0];const std::size_t i=device_indices[p];
        if(i>=n) return;
        auto& out=device_particles[p];out.global_id=gid[i];out.owner_rank=owner[i];out.x=x[i];out.y=y[i];out.z=z[i];out.vx=vx[i];out.vy=vy[i];out.vz=vz[i];out.wx=wx[i];out.wy=wy[i];out.wz=wz[i];out.radius=radius[i];out.inverse_mass=im[i];out.inverse_inertia=ii[i];
        for(std::size_t k=0;k<slots;++k){const std::size_t hs=i*slots+k;auto& h=device_history[p*slots+k];h.owner_global_id=gid[i];h.neighbor_global_id=hneighbor[hs];h.tx=htx[hs];h.ty=hty[hs];h.tz=htz[hs];}
    });
}

void ResidentDemSycl::append_packed_device(
    const ResidentDemPackedParticle* device_particles,
    const ResidentDemPackedHistory* device_history,
    std::size_t packed_count,
    int owner_rank_override) {
    if(packed_count==0U) return;
    if(!device_particles || !device_history) throw std::invalid_argument("resident DEM append requires non-null device pointers");
    if(count_+packed_count>config_.capacity) throw std::overflow_error("resident DEM append exceeds capacity");
    const auto base=count_,slots=config_.history_slots_per_particle;
    auto* gid=global_id_;auto* owner=owner_rank_;auto *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*radius=radius_,*im=inverse_mass_,*ii=inverse_inertia_;
    auto* hneighbor=history_neighbor_;auto* hseen=history_seen_;auto* htx=history_tx_;auto* hty=history_ty_;auto* htz=history_tz_;
    const auto stamp=history_stamp_;
    queue_.parallel_for(sycl::range<1>(packed_count),[=](sycl::id<1> id){
        const std::size_t p=id[0],i=base+p;const auto in=device_particles[p];
        gid[i]=in.global_id;owner[i]=owner_rank_override>=0?static_cast<std::int32_t>(owner_rank_override):in.owner_rank;
        x[i]=in.x;y[i]=in.y;z[i]=in.z;vx[i]=in.vx;vy[i]=in.vy;vz[i]=in.vz;wx[i]=in.wx;wy[i]=in.wy;wz[i]=in.wz;radius[i]=in.radius;im[i]=in.inverse_mass;ii[i]=in.inverse_inertia;
        for(std::size_t k=0;k<slots;++k){const auto h=device_history[p*slots+k];const std::size_t hs=i*slots+k;hneighbor[hs]=h.neighbor_global_id;htx[hs]=h.tx;hty[hs]=h.ty;htz[hs]=h.tz;hseen[hs]=h.neighbor_global_id==invalid_dem_global_id?0U:stamp;}
    });
    count_+=packed_count;
}

void ResidentDemSycl::retain_x_slab_device(double x_min,double x_max,int owner_rank) {
    if(!(x_max>x_min)) throw std::invalid_argument("resident DEM retain slab requires x_max>x_min");
    if(!compaction_count_){compaction_count_=sycl::malloc_shared<std::uint32_t>(1U,queue_);if(!compaction_count_)throw std::bad_alloc{};}
    *compaction_count_=0U;
    const auto n=count_,slots=config_.history_slots_per_particle;
    auto* gid=global_id_;auto* owner=owner_rank_;auto *x=x_,*y=y_,*z=z_,*vx=vx_,*vy=vy_,*vz=vz_,*wx=wx_,*wy=wy_,*wz=wz_,*radius=radius_,*im=inverse_mass_,*ii=inverse_inertia_;
    auto* fx=fx_;auto* fy=fy_;auto* fz=fz_;auto* tx=tx_;auto* ty=ty_;auto* tz=tz_;
    auto* hneighbor=history_neighbor_;auto* hseen=history_seen_;auto* htx=history_tx_;auto* hty=history_ty_;auto* htz=history_tz_;auto* out_count=compaction_count_;
    const float lo=static_cast<float>(x_min),hi=static_cast<float>(x_max);const auto local_owner=static_cast<std::int32_t>(owner_rank);
    queue_.parallel_for(sycl::range<1>(1U),[=](sycl::id<1>){
        std::uint32_t dst=0U;
        for(std::size_t src=0;src<n;++src){
            if(!(x[src]>=lo && x[src]<hi)) continue;
            if(dst!=src){gid[dst]=gid[src];owner[dst]=local_owner;x[dst]=x[src];y[dst]=y[src];z[dst]=z[src];vx[dst]=vx[src];vy[dst]=vy[src];vz[dst]=vz[src];wx[dst]=wx[src];wy[dst]=wy[src];wz[dst]=wz[src];radius[dst]=radius[src];im[dst]=im[src];ii[dst]=ii[src];fx[dst]=fx[src];fy[dst]=fy[src];fz[dst]=fz[src];tx[dst]=tx[src];ty[dst]=ty[src];tz[dst]=tz[src];
                for(std::size_t k=0;k<slots;++k){const auto hs=src*slots+k,hd=static_cast<std::size_t>(dst)*slots+k;hneighbor[hd]=hneighbor[hs];hseen[hd]=hseen[hs];htx[hd]=htx[hs];hty[hd]=hty[hs];htz[hd]=htz[hs];}
            }else owner[dst]=local_owner;
            ++dst;
        }
        *out_count=dst;
    }).wait_and_throw();
    count_=static_cast<std::size_t>(*compaction_count_);transfer_stats_.record_synchronization();
}

std::size_t ResidentDemSycl::resident_bytes() const noexcept {
    const std::size_t cap=config_.capacity,hist=cap*config_.history_slots_per_particle;
    return cap*(18U*sizeof(float)+2U*sizeof(std::uint32_t)+sizeof(std::uint64_t)+sizeof(std::int32_t))+
        (cell_count_*2U+cell_count_+1U)*sizeof(std::uint32_t)+
        hist*(3U*sizeof(float)+sizeof(std::uint64_t)+sizeof(std::uint32_t))+
        (reduction_scalar_?sizeof(float):0U)+(compaction_count_?sizeof(std::uint32_t):0U);
}
ResidentDemDeviceView ResidentDemSycl::device_view() noexcept { return {count_,global_id_,owner_rank_,x_,y_,z_,vx_,vy_,vz_,wx_,wy_,wz_,radius_,inverse_mass_,fx_,fy_,fz_,tx_,ty_,tz_}; }
ResidentDemExchangeDeviceView ResidentDemSycl::exchange_device_view() noexcept { return {count_,config_.capacity,config_.history_slots_per_particle,global_id_,owner_rank_,x_,y_,z_,vx_,vy_,vz_,wx_,wy_,wz_,radius_,inverse_mass_,inverse_inertia_,history_neighbor_,history_seen_,history_tx_,history_ty_,history_tz_}; }

} // namespace cfd::multibody
#endif
