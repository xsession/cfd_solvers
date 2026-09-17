#pragma once

#include "cfd/core/device_residency.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfd::fdtd {

enum class ResidentBoundary3D : unsigned char { pec, pmc, cpml };
enum class HaloFace3D : unsigned char { x_minus, x_plus, y_minus, y_plus, z_minus, z_plus };

struct ResidentMaxwell3DConfig {
    std::size_t nx{32U};
    std::size_t ny{32U};
    std::size_t nz{32U};
    float dx{1.0e-3F};
    float dy{1.0e-3F};
    float dz{1.0e-3F};
    float courant{0.8F};
    float epsilon_r{1.0F};
    float mu_r{1.0F};
    ResidentBoundary3D boundary{ResidentBoundary3D::cpml};
    std::size_t cpml_cells{8U};
    float cpml_order{3.0F};
    float cpml_target_reflection{1.0e-7F};
    float cpml_kappa_max{5.0F};
    float cpml_alpha_fraction{0.05F};
    // Set a face false for an MPI/internal subdomain interface. CPML and PEC/PMC
    // enforcement are then omitted on that face so device halo exchange can own it.
    std::array<bool,6> physical_faces{true,true,true,true,true,true};
};

struct ResidentFieldView3D {
    float* ex{}; float* ey{}; float* ez{};
    float* hx{}; float* hy{}; float* hz{};
    std::size_t cells{};
};

// Device-resident 3-D Yee FDTD baseline. E/H, material coefficients and the
// CPML auxiliary fields stay in USM device allocations across the hot loop.
// Host transfers are restricted to setup, explicit output and scalar API
// boundaries; halo pack/unpack accepts caller-owned device buffers.
class ResidentMaxwell3DSycl {
public:
    static constexpr float epsilon0 = 8.8541878128e-12F;
    static constexpr float mu0 = 1.25663706212e-6F;

    explicit ResidentMaxwell3DSycl(ResidentMaxwell3DConfig config)
        : config_(validate(config)),
          cells_(config_.nx * config_.ny * config_.nz),
          queue_(sycl::default_selector_v, sycl::property::queue::in_order{}) {
        const float epsilon_ref = epsilon0 * config_.epsilon_r;
        const float mu_ref = mu0 * config_.mu_r;
        const float wave = 1.0F / sycl::sqrt(epsilon_ref * mu_ref);
        const float metric = sycl::sqrt(1.0F/(config_.dx*config_.dx) +
                                        1.0F/(config_.dy*config_.dy) +
                                        1.0F/(config_.dz*config_.dz));
        dt_ = config_.courant / (wave * metric);

        fields_ = sycl::malloc_device<float>(6U * cells_, queue_);
        material_ = sycl::malloc_device<float>(5U * cells_, queue_);
        psi_e_ = sycl::malloc_device<float>(6U * cells_, queue_);
        psi_h_ = sycl::malloc_device<float>(6U * cells_, queue_);
        profile_x_ = sycl::malloc_device<float>(6U * config_.nx, queue_);
        profile_y_ = sycl::malloc_device<float>(6U * config_.ny, queue_);
        profile_z_ = sycl::malloc_device<float>(6U * config_.nz, queue_);
        reduction_scalar_ = sycl::malloc_shared<float>(1U, queue_);
        if (!fields_ || !material_ || !psi_e_ || !psi_h_ || !profile_x_ || !profile_y_ ||
            !profile_z_ || !reduction_scalar_) throw std::bad_alloc{};

        queue_.memset(fields_, 0, 6U * cells_ * sizeof(float));
        queue_.memset(psi_e_, 0, 6U * cells_ * sizeof(float));
        queue_.memset(psi_h_, 0, 6U * cells_ * sizeof(float));
        initialize_uniform_material(config_.epsilon_r, config_.epsilon_r, config_.epsilon_r,
                                    config_.mu_r, 0.0F);
        upload_cpml_profiles();
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
    }

    ~ResidentMaxwell3DSycl() noexcept {
        try { queue_.wait_and_throw(); } catch (...) {}
        if (fields_) sycl::free(fields_, queue_);
        if (material_) sycl::free(material_, queue_);
        if (psi_e_) sycl::free(psi_e_, queue_);
        if (psi_h_) sycl::free(psi_h_, queue_);
        if (profile_x_) sycl::free(profile_x_, queue_);
        if (profile_y_) sycl::free(profile_y_, queue_);
        if (profile_z_) sycl::free(profile_z_, queue_);
        if (reduction_scalar_) sycl::free(reduction_scalar_, queue_);
    }

    ResidentMaxwell3DSycl(const ResidentMaxwell3DSycl&) = delete;
    ResidentMaxwell3DSycl& operator=(const ResidentMaxwell3DSycl&) = delete;

    [[nodiscard]] sycl::queue& queue() noexcept { return queue_; }
    [[nodiscard]] const sycl::queue& queue() const noexcept { return queue_; }
    [[nodiscard]] float dt() const noexcept { return dt_; }
    [[nodiscard]] std::uint64_t time_step() const noexcept { return time_step_; }
    [[nodiscard]] std::size_t cells() const noexcept { return cells_; }
    [[nodiscard]] const ResidentMaxwell3DConfig& config() const noexcept { return config_; }
    [[nodiscard]] const cfd::core::DeviceTransferStats& transfer_stats() const noexcept { return transfer_stats_; }
    void reset_transfer_stats() const noexcept { transfer_stats_.reset(); }

    [[nodiscard]] ResidentFieldView3D device_fields() noexcept {
        return {fields_, fields_ + cells_, fields_ + 2U*cells_, fields_ + 3U*cells_,
                fields_ + 4U*cells_, fields_ + 5U*cells_, cells_};
    }

    void initialize_gaussian_ez(float amplitude = 1.0F, float width_fraction = 0.12F) {
        if (!std::isfinite(amplitude) || !(width_fraction > 0.0F))
            throw std::invalid_argument("invalid resident FDTD Gaussian parameters");
        const std::size_t nx=config_.nx, ny=config_.ny, nz=config_.nz, cells=cells_;
        const float cx=0.5F*static_cast<float>(nx-1U), cy=0.5F*static_cast<float>(ny-1U),
                    cz=0.5F*static_cast<float>(nz-1U);
        const float width=width_fraction*static_cast<float>(std::min({nx,ny,nz}));
        float* state=fields_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
            const std::size_t n=gid[0];
            const std::size_t i=n%nx, j=(n/nx)%ny, k=n/(nx*ny);
            const float x=(static_cast<float>(i)-cx)/width;
            const float y=(static_cast<float>(j)-cy)/width;
            const float z=(static_cast<float>(k)-cz)/width;
            for (std::size_t c=0;c<6U;++c) state[c*cells+n]=0.0F;
            state[2U*cells+n]=amplitude*sycl::exp(-0.5F*(x*x+y*y+z*z));
        });
        enforce_outer_boundary();
    }

    // Material update is a device kernel, so changing an internal box does not
    // stage full material fields through the CPU. The CPML layer should remain
    // in the reference background material for impedance matching.
    void set_material_box(std::size_t i0,std::size_t i1,std::size_t j0,std::size_t j1,
                          std::size_t k0,std::size_t k1,float epsilon_rx,float epsilon_ry,
                          float epsilon_rz,float mu_r,float conductivity_s_per_m=0.0F) {
        if (i0>=i1||j0>=j1||k0>=k1||i1>config_.nx||j1>config_.ny||k1>config_.nz||
            !(epsilon_rx>0.0F)||!(epsilon_ry>0.0F)||!(epsilon_rz>0.0F)||!(mu_r>0.0F)||
            !(conductivity_s_per_m>=0.0F)||!std::isfinite(epsilon_rx)||!std::isfinite(epsilon_ry)||
            !std::isfinite(epsilon_rz)||!std::isfinite(mu_r)||!std::isfinite(conductivity_s_per_m))
            throw std::invalid_argument("invalid resident FDTD material box");
        if(config_.boundary==ResidentBoundary3D::cpml){
            const auto l=config_.cpml_cells;
            if((config_.physical_faces[0]&&i0<l)||(config_.physical_faces[1]&&i1>config_.nx-l)||
               (config_.physical_faces[2]&&j0<l)||(config_.physical_faces[3]&&j1>config_.ny-l)||
               (config_.physical_faces[4]&&k0<l)||(config_.physical_faces[5]&&k1>config_.nz-l))
                throw std::invalid_argument("resident FDTD material box overlaps a physical CPML layer");
        }
        const std::size_t nx=config_.nx,ny=config_.ny,cells=cells_;
        float* material=material_;
        const float ix=1.0F/(epsilon0*epsilon_rx), iy=1.0F/(epsilon0*epsilon_ry),
                    iz=1.0F/(epsilon0*epsilon_rz), im=1.0F/(mu0*mu_r), sigma=conductivity_s_per_m;
        const std::size_t count=(i1-i0)*(j1-j0)*(k1-k0);
        queue_.parallel_for(sycl::range<1>(count), [=](sycl::id<1> gid) {
            const std::size_t q=gid[0];
            const std::size_t ni=i1-i0,nj=j1-j0;
            const std::size_t i=i0+q%ni, j=j0+(q/ni)%nj, k=k0+q/(ni*nj);
            const std::size_t n=(k*ny+j)*nx+i;
            material[n]=ix; material[cells+n]=iy; material[2U*cells+n]=iz;
            material[3U*cells+n]=im; material[4U*cells+n]=sigma;
        });
    }

    void add_soft_ez_source(std::size_t i,std::size_t j,std::size_t k,float value) {
        if(i>=config_.nx||j>=config_.ny||k>=config_.nz||!std::isfinite(value))
            throw std::out_of_range("invalid resident FDTD source");
        const std::size_t n=(k*config_.ny+j)*config_.nx+i,cells=cells_;
        float* state=fields_;
        queue_.parallel_for(sycl::range<1>(1U), [=](sycl::id<1>) { state[2U*cells+n]+=value; });
    }

    void step(std::size_t count=1U) {
        for(std::size_t s=0;s<count;++s) {
            update_magnetic();
            update_electric();
            enforce_outer_boundary();
            ++time_step_;
        }
    }

    void wait() const {
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
    }

    // Scalar reduction stays on-device except for one shared scalar used as the
    // API result. This intentionally does not count as a bulk field transfer.
    [[nodiscard]] double energy() const {
        *reduction_scalar_=0.0F;
        const float* state=fields_; const float* material=material_;
        const std::size_t cells=cells_; const float dv=config_.dx*config_.dy*config_.dz;
        auto reduction=sycl::reduction(reduction_scalar_,sycl::plus<float>());
        queue_.parallel_for(sycl::range<1>(cells),reduction,[=](sycl::id<1> gid,auto& sum){
            const std::size_t n=gid[0];
            const float ex=state[n],ey=state[cells+n],ez=state[2U*cells+n];
            const float hx=state[3U*cells+n],hy=state[4U*cells+n],hz=state[5U*cells+n];
            const float density=0.5F*(ex*ex/material[n]+ey*ey/material[cells+n]+ez*ez/material[2U*cells+n]+
                                      (hx*hx+hy*hy+hz*hz)/material[3U*cells+n]);
            sum.combine(density*dv);
        }).wait_and_throw();
        transfer_stats_.record_synchronization();
        return static_cast<double>(*reduction_scalar_);
    }

    void capture_ez_to_device(std::size_t i,std::size_t j,std::size_t k,float* device_samples,std::size_t slot) {
        if(!device_samples||i>=config_.nx||j>=config_.ny||k>=config_.nz)
            throw std::invalid_argument("invalid resident FDTD probe destination");
        const std::size_t n=(k*config_.ny+j)*config_.nx+i,cells=cells_;
        const float* state=fields_;
        queue_.parallel_for(sycl::range<1>(1U), [=](sycl::id<1>) { device_samples[slot]=state[2U*cells+n]; });
    }

    [[nodiscard]] std::size_t halo_cells(HaloFace3D face) const noexcept {
        switch(face) {
            case HaloFace3D::x_minus: case HaloFace3D::x_plus: return config_.ny*config_.nz;
            case HaloFace3D::y_minus: case HaloFace3D::y_plus: return config_.nx*config_.nz;
            case HaloFace3D::z_minus: case HaloFace3D::z_plus: return config_.nx*config_.ny;
        }
        return 0U;
    }
    [[nodiscard]] std::size_t halo_value_count(HaloFace3D face) const noexcept { return 6U*halo_cells(face); }

    // Pack/unpack caller-owned device buffers. This is the seam for MPI-aware
    // GPU-direct or staged halo exchange; no host field buffer is required.
    void pack_halo(HaloFace3D face,float* device_buffer) {
        if(!device_buffer) throw std::invalid_argument("null FDTD halo buffer");
        const std::size_t nx=config_.nx,ny=config_.ny,nz=config_.nz,cells=cells_,count=halo_cells(face);
        const float* state=fields_; const unsigned f=static_cast<unsigned>(face);
        queue_.parallel_for(sycl::range<1>(count), [=](sycl::id<1> gid) {
            const std::size_t q=gid[0]; std::size_t i=0,j=0,k=0;
            if(f<=1U){j=q%ny;k=q/ny;i=(f==0U?0U:nx-1U);} 
            else if(f<=3U){i=q%nx;k=q/nx;j=(f==2U?0U:ny-1U);} 
            else {i=q%nx;j=q/nx;k=(f==4U?0U:nz-1U);} 
            const std::size_t n=(k*ny+j)*nx+i;
            for(std::size_t c=0;c<6U;++c) device_buffer[c*count+q]=state[c*cells+n];
        });
    }

    void unpack_halo(HaloFace3D face,const float* device_buffer) {
        if(!device_buffer) throw std::invalid_argument("null FDTD halo buffer");
        const std::size_t nx=config_.nx,ny=config_.ny,nz=config_.nz,cells=cells_,count=halo_cells(face);
        float* state=fields_; const unsigned f=static_cast<unsigned>(face);
        queue_.parallel_for(sycl::range<1>(count), [=](sycl::id<1> gid) {
            const std::size_t q=gid[0]; std::size_t i=0,j=0,k=0;
            if(f<=1U){j=q%ny;k=q/ny;i=(f==0U?0U:nx-1U);} 
            else if(f<=3U){i=q%nx;k=q/nx;j=(f==2U?0U:ny-1U);} 
            else {i=q%nx;j=q/nx;k=(f==4U?0U:nz-1U);} 
            const std::size_t n=(k*ny+j)*nx+i;
            for(std::size_t c=0;c<6U;++c) state[c*cells+n]=device_buffer[c*count+q];
        });
    }

    [[nodiscard]] std::size_t resident_bytes() const noexcept {
        return (6U+5U+6U+6U)*cells_*sizeof(float) +
               6U*(config_.nx+config_.ny+config_.nz)*sizeof(float) + sizeof(float);
    }

private:
    ResidentMaxwell3DConfig config_{};
    std::size_t cells_{};
    mutable sycl::queue queue_;
    float dt_{};
    std::uint64_t time_step_{};
    float* fields_{nullptr};
    float* material_{nullptr};
    float* psi_e_{nullptr};
    float* psi_h_{nullptr};
    float* profile_x_{nullptr};
    float* profile_y_{nullptr};
    float* profile_z_{nullptr};
    mutable float* reduction_scalar_{nullptr};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};

    static ResidentMaxwell3DConfig validate(ResidentMaxwell3DConfig c) {
        if(c.nx<6U||c.ny<6U||c.nz<6U||!(c.dx>0.0F)||!(c.dy>0.0F)||!(c.dz>0.0F)||
           !(c.courant>0.0F&&c.courant<1.0F)||!(c.epsilon_r>0.0F)||!(c.mu_r>0.0F))
            throw std::invalid_argument("invalid resident Maxwell3D configuration");
        if(c.boundary==ResidentBoundary3D::cpml &&
           (c.cpml_cells<2U || 2U*c.cpml_cells+2U>std::min({c.nx,c.ny,c.nz}) ||
            !(c.cpml_order>=1.0F) || !(c.cpml_target_reflection>0.0F&&c.cpml_target_reflection<1.0F) ||
            !(c.cpml_kappa_max>=1.0F) || !(c.cpml_alpha_fraction>=0.0F)))
            throw std::invalid_argument("invalid resident Maxwell3D CPML configuration");
        return c;
    }

    void initialize_uniform_material(float erx,float ery,float erz,float mur,float sigma) {
        const float ix=1.0F/(epsilon0*erx),iy=1.0F/(epsilon0*ery),iz=1.0F/(epsilon0*erz),im=1.0F/(mu0*mur);
        const std::size_t cells=cells_; float* material=material_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid){
            const std::size_t n=gid[0]; material[n]=ix; material[cells+n]=iy; material[2U*cells+n]=iz;
            material[3U*cells+n]=im; material[4U*cells+n]=sigma;
        });
    }

    static double depth_e(std::size_t i,std::size_t n,std::size_t layer,bool minus_face,bool plus_face) noexcept {
        if(minus_face&&i<layer) return (static_cast<double>(layer)-static_cast<double>(i)-0.5)/static_cast<double>(layer);
        if(plus_face&&i+layer>=n) return (static_cast<double>(i)-static_cast<double>(n-layer)+0.5)/static_cast<double>(layer);
        return 0.0;
    }
    static double depth_h(std::size_t i,std::size_t n,std::size_t layer,bool minus_face,bool plus_face) noexcept {
        const double x=static_cast<double>(i)+0.5;
        if(minus_face&&x<static_cast<double>(layer)) return (static_cast<double>(layer)-x)/static_cast<double>(layer);
        if(plus_face&&x>static_cast<double>(n-layer)) return (x-static_cast<double>(n-layer))/static_cast<double>(layer);
        return 0.0;
    }

    std::vector<float> make_profile(std::size_t n,float spacing,bool minus_face,bool plus_face) const {
        std::vector<float> out(6U*n,1.0F);
        for(std::size_t i=0;i<n;++i){out[2U*n+i]=0.0F;out[5U*n+i]=0.0F;}
        if(config_.boundary!=ResidentBoundary3D::cpml) return out;
        const double eps=static_cast<double>(epsilon0*config_.epsilon_r),mu=static_cast<double>(mu0*config_.mu_r);
        const double eta=std::sqrt(mu/eps),thickness=static_cast<double>(config_.cpml_cells)*spacing;
        const double sigma_max=-(static_cast<double>(config_.cpml_order)+1.0)*std::log(config_.cpml_target_reflection)/(2.0*eta*thickness);
        const double alpha_max=static_cast<double>(config_.cpml_alpha_fraction)*sigma_max;
        auto fill=[&](std::size_t i,double depth,bool magnetic){
            depth=std::clamp(depth,0.0,1.0); if(!(depth>0.0)) return;
            const double grade=std::pow(depth,config_.cpml_order);
            const double sigma_e=sigma_max*grade,kappa=1.0+(config_.cpml_kappa_max-1.0)*grade,alpha_e=alpha_max*(1.0-depth);
            const double sigma=magnetic?sigma_e*mu/eps:sigma_e,alpha=magnetic?alpha_e*mu/eps:alpha_e,den=magnetic?mu:eps;
            const double b=std::exp(-(sigma/kappa+alpha)*static_cast<double>(dt_)/den);
            const double denom=sigma*kappa+kappa*kappa*alpha;
            const double c=denom>0.0?sigma*(b-1.0)/denom:0.0;
            const std::size_t base=magnetic?3U*n:0U;
            out[base+i]=static_cast<float>(kappa);out[base+n+i]=static_cast<float>(b);out[base+2U*n+i]=static_cast<float>(c);
        };
        for(std::size_t i=0;i<n;++i){fill(i,depth_e(i,n,config_.cpml_cells,minus_face,plus_face),false);fill(i,depth_h(i,n,config_.cpml_cells,minus_face,plus_face),true);}
        return out;
    }

    void upload_cpml_profiles() {
        const auto px=make_profile(config_.nx,config_.dx,config_.physical_faces[0],config_.physical_faces[1]),
                   py=make_profile(config_.ny,config_.dy,config_.physical_faces[2],config_.physical_faces[3]),
                   pz=make_profile(config_.nz,config_.dz,config_.physical_faces[4],config_.physical_faces[5]);
        queue_.memcpy(profile_x_,px.data(),px.size()*sizeof(float));
        queue_.memcpy(profile_y_,py.data(),py.size()*sizeof(float));
        queue_.memcpy(profile_z_,pz.data(),pz.size()*sizeof(float));
        transfer_stats_.record_host_to_device((px.size()+py.size()+pz.size())*sizeof(float));
    }

    void update_magnetic() {
        const std::size_t nx=config_.nx,ny=config_.ny,nz=config_.nz,cells=cells_;
        const float dx=config_.dx,dy=config_.dy,dz=config_.dz,dt=dt_;
        float* s=fields_; const float* mat=material_; float* psi=psi_h_;
        const float* px=profile_x_;const float* py=profile_y_;const float* pz=profile_z_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid){
            const std::size_t n=gid[0],i=n%nx,j=(n/nx)%ny,k=n/(nx*ny);
            if(i+1U>=nx||j+1U>=ny||k+1U>=nz) return;
            auto cp=[&](float d,float* p,const float* profile,std::size_t coord,std::size_t extent,bool magnetic){
                const std::size_t base=magnetic?3U*extent:0U;
                const float kap=profile[base+coord],b=profile[base+extent+coord],c=profile[base+2U*extent+coord];
                *p=b*(*p)+c*d; return d/kap+(*p);
            };
            const std::size_t xp=n+1U,yp=n+nx,zp=n+nx*ny;
            float dezdy=(s[2U*cells+yp]-s[2U*cells+n])/dy;
            float deydz=(s[cells+zp]-s[cells+n])/dz;
            dezdy=cp(dezdy,&psi[n],py,j,ny,true);deydz=cp(deydz,&psi[cells+n],pz,k,nz,true);
            float dexdz=(s[zp]-s[n])/dz;
            float dezdx=(s[2U*cells+xp]-s[2U*cells+n])/dx;
            dexdz=cp(dexdz,&psi[2U*cells+n],pz,k,nz,true);dezdx=cp(dezdx,&psi[3U*cells+n],px,i,nx,true);
            float deydx=(s[cells+xp]-s[cells+n])/dx;
            float dexdy=(s[yp]-s[n])/dy;
            deydx=cp(deydx,&psi[4U*cells+n],px,i,nx,true);dexdy=cp(dexdy,&psi[5U*cells+n],py,j,ny,true);
            const float im=mat[3U*cells+n];
            s[3U*cells+n]-=dt*im*(dezdy-deydz);
            s[4U*cells+n]-=dt*im*(dexdz-dezdx);
            s[5U*cells+n]-=dt*im*(deydx-dexdy);
        });
    }

    void update_electric() {
        const std::size_t nx=config_.nx,ny=config_.ny,nz=config_.nz,cells=cells_;
        const float dx=config_.dx,dy=config_.dy,dz=config_.dz,dt=dt_;
        float* s=fields_; const float* mat=material_; float* psi=psi_e_;
        const float* px=profile_x_;const float* py=profile_y_;const float* pz=profile_z_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid){
            const std::size_t n=gid[0],i=n%nx,j=(n/nx)%ny,k=n/(nx*ny);
            if(i==0U||j==0U||k==0U||i+1U>=nx||j+1U>=ny||k+1U>=nz) return;
            auto cp=[&](float d,float* p,const float* profile,std::size_t coord,std::size_t extent){
                const float kap=profile[coord],b=profile[extent+coord],c=profile[2U*extent+coord];
                *p=b*(*p)+c*d; return d/kap+(*p);
            };
            const std::size_t xm=n-1U,ym=n-nx,zm=n-nx*ny;
            float dhzdy=(s[5U*cells+n]-s[5U*cells+ym])/dy;
            float dhydz=(s[4U*cells+n]-s[4U*cells+zm])/dz;
            dhzdy=cp(dhzdy,&psi[n],py,j,ny);dhydz=cp(dhydz,&psi[cells+n],pz,k,nz);
            float dhxdz=(s[3U*cells+n]-s[3U*cells+zm])/dz;
            float dhzdx=(s[5U*cells+n]-s[5U*cells+xm])/dx;
            dhxdz=cp(dhxdz,&psi[2U*cells+n],pz,k,nz);dhzdx=cp(dhzdx,&psi[3U*cells+n],px,i,nx);
            float dhydx=(s[4U*cells+n]-s[4U*cells+xm])/dx;
            float dhxdy=(s[3U*cells+n]-s[3U*cells+ym])/dy;
            dhydx=cp(dhydx,&psi[4U*cells+n],px,i,nx);dhxdy=cp(dhxdy,&psi[5U*cells+n],py,j,ny);
            const float sigma=mat[4U*cells+n];
            const float invx=mat[n],invy=mat[cells+n],invz=mat[2U*cells+n];
            const float lx=0.5F*sigma*dt*invx,ly=0.5F*sigma*dt*invy,lz=0.5F*sigma*dt*invz;
            s[n]=((1.0F-lx)*s[n]+dt*invx*(dhzdy-dhydz))/(1.0F+lx);
            s[cells+n]=((1.0F-ly)*s[cells+n]+dt*invy*(dhxdz-dhzdx))/(1.0F+ly);
            s[2U*cells+n]=((1.0F-lz)*s[2U*cells+n]+dt*invz*(dhydx-dhxdy))/(1.0F+lz);
        });
    }

    void enforce_outer_boundary() {
        const std::size_t nx=config_.nx,ny=config_.ny,nz=config_.nz,cells=cells_; float* s=fields_;
        const bool pec=config_.boundary!=ResidentBoundary3D::pmc;
        const bool xm=config_.physical_faces[0],xp=config_.physical_faces[1],ym=config_.physical_faces[2],
                   yp=config_.physical_faces[3],zm=config_.physical_faces[4],zp=config_.physical_faces[5];
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid){
            const std::size_t n=gid[0],i=n%nx,j=(n/nx)%ny,k=n/(nx*ny);
            if((xm&&i==0U)||(xp&&i+1U==nx)||(ym&&j==0U)||(yp&&j+1U==ny)||(zm&&k==0U)||(zp&&k+1U==nz)){
                const std::size_t base=pec?0U:3U;
                s[(base+0U)*cells+n]=0.0F;s[(base+1U)*cells+n]=0.0F;s[(base+2U)*cells+n]=0.0F;
            }
        });
    }
};

} // namespace cfd::fdtd
#endif
