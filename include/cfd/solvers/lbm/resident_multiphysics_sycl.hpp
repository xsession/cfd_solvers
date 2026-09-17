#pragma once

#include "cfd/core/device_residency.hpp"
#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"
#include "cfd/solvers/lbm/particles.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <new>
#include <span>
#include <stdexcept>
#include <vector>

namespace cfd::lbm {

namespace resident_detail {
inline std::size_t shift(std::size_t v, int delta, std::size_t n) noexcept {
    if (delta > 0) return v + 1U == n ? 0U : v + 1U;
    if (delta < 0) return v == 0U ? n - 1U : v - 1U;
    return v;
}
inline std::size_t index(std::size_t x, std::size_t y, std::size_t z,
                         std::size_t nx, std::size_t ny) noexcept {
    return (z * ny + y) * nx + x;
}
inline long long floor_to_ll(float value) noexcept {
    long long i = static_cast<long long>(value);
    if (static_cast<float>(i) > value) --i;
    return i;
}
inline float periodic_coordinate(float value, std::size_t extent) noexcept {
    const float e = static_cast<float>(extent);
    const long long turns = floor_to_ll(value / e);
    value -= static_cast<float>(turns) * e;
    if (value >= e) value -= e;
    if (value < 0.0F) value += e;
    return value;
}
inline std::size_t wrap_cell(long long v, std::size_t extent) noexcept {
    const long long e = static_cast<long long>(extent);
    v %= e;
    if (v < 0) v += e;
    return static_cast<std::size_t>(v);
}
} // namespace resident_detail

struct ResidentThermalConfig {
    std::size_t nx{64U}, ny{64U}, nz{64U};
    float diffusivity{0.01F};
};

// Device-resident D3Q7 passive-scalar / temperature lattice. The velocity is
// consumed directly from EsotericPullSyclSolver macroscopic scratch; no host
// field transfer is required between hydrodynamics and thermal transport.
class ResidentThermalD3Q7Sycl {
public:
    explicit ResidentThermalD3Q7Sycl(sycl::queue& queue, ResidentThermalConfig config)
        : queue_(queue), config_(config), cells_(config.nx * config.ny * config.nz),
          tau_(0.5F + 4.0F * config.diffusivity) {
        if (config.nx < 2U || config.ny < 2U || config.nz < 2U || !(config.diffusivity > 0.0F)) {
            throw std::invalid_argument("invalid resident thermal lattice");
        }
        g_ = sycl::malloc_device<float>(7U * cells_, queue_);
        next_ = sycl::malloc_device<float>(7U * cells_, queue_);
        if (!g_ || !next_) throw std::bad_alloc{};
        initialize_uniform(0.0F);
    }

    ~ResidentThermalD3Q7Sycl() noexcept {
        try { queue_.wait_and_throw(); } catch (...) {}
        if (g_) sycl::free(g_, queue_);
        if (next_) sycl::free(next_, queue_);
        if (scalar_scratch_) sycl::free(scalar_scratch_, queue_);
        if (reduction_scalar_) sycl::free(reduction_scalar_, queue_);
    }

    ResidentThermalD3Q7Sycl(const ResidentThermalD3Q7Sycl&) = delete;
    ResidentThermalD3Q7Sycl& operator=(const ResidentThermalD3Q7Sycl&) = delete;

    void initialize_uniform(float temperature) {
        if (!std::isfinite(temperature)) throw std::invalid_argument("thermal state must be finite");
        float* g = g_;
        const std::size_t cells = cells_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n = id[0];
            for (int d = 0; d < 7; ++d) g[static_cast<std::size_t>(d) * cells + n] = weight(d) * temperature;
        });
    }

    void step(const D3Q19SyclSolver::DeviceMacroscopicView& flow, std::size_t count = 1U) {
        if (flow.cells != cells_) throw std::invalid_argument("thermal/flow lattice size mismatch");
        const std::size_t nx = config_.nx, ny = config_.ny, nz = config_.nz, cells = cells_;
        const float omega = 1.0F / tau_;
        const float* ux = flow.ux; const float* uy = flow.uy; const float* uz = flow.uz;
        for (std::size_t iteration = 0; iteration < count; ++iteration) {
            const float* in = g_;
            float* out = next_;
            queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
                const std::size_t n = id[0];
                const std::size_t x = n % nx;
                const std::size_t yz = n / nx;
                const std::size_t y = yz % ny;
                const std::size_t z = yz / ny;
                float local[7];
                float temperature = 0.0F;
                for (int d = 0; d < 7; ++d) {
                    const std::size_t xs = resident_detail::shift(x, -cx(d), nx);
                    const std::size_t ys = resident_detail::shift(y, -cy(d), ny);
                    const std::size_t zs = resident_detail::shift(z, -cz(d), nz);
                    const std::size_t src = resident_detail::index(xs, ys, zs, nx, ny);
                    local[d] = in[static_cast<std::size_t>(d) * cells + src];
                    temperature += local[d];
                }
                for (int d = 0; d < 7; ++d) {
                    const float cu = static_cast<float>(cx(d)) * ux[n] +
                                     static_cast<float>(cy(d)) * uy[n] +
                                     static_cast<float>(cz(d)) * uz[n];
                    const float geq = weight(d) * temperature * (1.0F + 4.0F * cu);
                    out[static_cast<std::size_t>(d) * cells + n] = local[d] - omega * (local[d] - geq);
                }
            });
            float* old = g_; g_ = next_; next_ = old;
        }
    }

    [[nodiscard]] float* refresh_device_scalar() const {
        ensure_scalar_scratch();
        const float* g = g_; float* scalar = scalar_scratch_; const std::size_t cells = cells_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n = id[0]; float value = 0.0F;
            for (int d = 0; d < 7; ++d) value += g[static_cast<std::size_t>(d) * cells + n];
            scalar[n] = value;
        });
        return scalar_scratch_;
    }

    void add_boussinesq_acceleration(D3Q19SyclSolver::DeviceAccelerationView force,
                                     float reference_temperature, float beta,
                                     ParticleVector3 gravity) const {
        if (force.cells != cells_ || !std::isfinite(reference_temperature) || !std::isfinite(beta)) {
            throw std::invalid_argument("invalid resident Boussinesq coupling");
        }
        const float* temperature = refresh_device_scalar();
        float* ax = force.ax; float* ay = force.ay; float* az = force.az;
        const float gx = static_cast<float>(gravity.x), gy = static_cast<float>(gravity.y), gz = static_cast<float>(gravity.z);
        const std::size_t cells = cells_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n = id[0];
            const float scale = -beta * (temperature[n] - reference_temperature);
            ax[n] += scale * gx; ay[n] += scale * gy; az[n] += scale * gz;
        });
    }

    [[nodiscard]] double total_scalar() const {
        ensure_reduction_scalar();
        *reduction_scalar_ = 0.0F;
        const float* g = g_; const std::size_t count = 7U * cells_;
        auto reduction = sycl::reduction(reduction_scalar_, sycl::plus<float>());
        queue_.parallel_for(sycl::range<1>(count), reduction, [=](sycl::id<1> id, auto& sum) {
            sum.combine(g[id[0]]);
        }).wait_and_throw();
        transfer_stats_.record_synchronization();
        return static_cast<double>(*reduction_scalar_);
    }

    [[nodiscard]] std::size_t resident_bytes() const noexcept {
        return 14U * cells_ * sizeof(float) + (scalar_scratch_ ? cells_ * sizeof(float) : 0U) +
               (reduction_scalar_ ? sizeof(float) : 0U);
    }
    [[nodiscard]] const cfd::core::DeviceTransferStats& transfer_stats() const noexcept { return transfer_stats_; }
    void reset_transfer_stats() const noexcept { transfer_stats_.reset(); }

private:
    sycl::queue& queue_;
    ResidentThermalConfig config_;
    std::size_t cells_{};
    float tau_{};
    float* g_{nullptr};
    float* next_{nullptr};
    mutable float* scalar_scratch_{nullptr};
    mutable float* reduction_scalar_{nullptr};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};

    static constexpr int cx(int d) noexcept { return d == 1 ? 1 : d == 2 ? -1 : 0; }
    static constexpr int cy(int d) noexcept { return d == 3 ? 1 : d == 4 ? -1 : 0; }
    static constexpr int cz(int d) noexcept { return d == 5 ? 1 : d == 6 ? -1 : 0; }
    static constexpr float weight(int d) noexcept { return d == 0 ? 0.25F : 0.125F; }
    void ensure_scalar_scratch() const {
        if (!scalar_scratch_) scalar_scratch_ = sycl::malloc_device<float>(cells_, queue_);
        if (!scalar_scratch_) throw std::bad_alloc{};
    }
    void ensure_reduction_scalar() const {
        if (!reduction_scalar_) reduction_scalar_ = sycl::malloc_shared<float>(1U, queue_);
        if (!reduction_scalar_) throw std::bad_alloc{};
    }
};

class ResidentFreeSurfaceSycl3D {
public:
    ResidentFreeSurfaceSycl3D(sycl::queue& queue, std::size_t nx, std::size_t ny, std::size_t nz)
        : queue_(queue), nx_(nx), ny_(ny), nz_(nz), cells_(nx * ny * nz) {
        if (nx < 2U || ny < 2U || nz < 2U) throw std::invalid_argument("invalid resident free-surface lattice");
        fill_ = sycl::malloc_device<float>(cells_, queue_);
        next_ = sycl::malloc_device<float>(cells_, queue_);
        geometry_ = sycl::malloc_device<float>(4U * cells_, queue_);
        if (!fill_ || !next_ || !geometry_) throw std::bad_alloc{};
        initialize_uniform(0.0F);
    }
    ~ResidentFreeSurfaceSycl3D() noexcept {
        try { queue_.wait_and_throw(); } catch (...) {}
        if (fill_) sycl::free(fill_, queue_);
        if (next_) sycl::free(next_, queue_);
        if (geometry_) sycl::free(geometry_, queue_);
        if (reduction_scalar_) sycl::free(reduction_scalar_, queue_);
    }
    ResidentFreeSurfaceSycl3D(const ResidentFreeSurfaceSycl3D&) = delete;
    ResidentFreeSurfaceSycl3D& operator=(const ResidentFreeSurfaceSycl3D&) = delete;

    void initialize_uniform(float fill) {
        if (!(fill >= 0.0F && fill <= 1.0F)) throw std::invalid_argument("free-surface fill must be in [0,1]");
        float* field = fill_; const std::size_t cells = cells_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) { field[id[0]] = fill; });
    }
    void initialize_half_space_x(float interface_x) {
        if (!std::isfinite(interface_x)) throw std::invalid_argument("invalid free-surface interface position");
        float* field = fill_; const std::size_t nx = nx_, cells = cells_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n = id[0], x = n % nx;
            const float left = static_cast<float>(x), right = left + 1.0F;
            float value = interface_x <= left ? 0.0F : interface_x >= right ? 1.0F : interface_x - left;
            field[n] = value;
        });
    }

    void advect(const D3Q19SyclSolver::DeviceMacroscopicView& flow, float dt) {
        if (flow.cells != cells_ || !(dt > 0.0F)) throw std::invalid_argument("invalid resident free-surface advection");
        const std::size_t nx = nx_, ny = ny_, nz = nz_, cells = cells_;
        const float* phi = fill_; float* out = next_;
        const float* ux = flow.ux; const float* uy = flow.uy; const float* uz = flow.uz;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n = id[0], x = n % nx, yz = n / nx, y = yz % ny, z = yz / ny;
            const std::size_t xp = resident_detail::index(resident_detail::shift(x,1,nx),y,z,nx,ny);
            const std::size_t xm = resident_detail::index(resident_detail::shift(x,-1,nx),y,z,nx,ny);
            const std::size_t yp = resident_detail::index(x,resident_detail::shift(y,1,ny),z,nx,ny);
            const std::size_t ym = resident_detail::index(x,resident_detail::shift(y,-1,ny),z,nx,ny);
            const std::size_t zp = resident_detail::index(x,y,resident_detail::shift(z,1,nz),nx,ny);
            const std::size_t zm = resident_detail::index(x,y,resident_detail::shift(z,-1,nz),nx,ny);
            const float uxp = 0.5F * (ux[n] + ux[xp]), uxm = 0.5F * (ux[xm] + ux[n]);
            const float uyp = 0.5F * (uy[n] + uy[yp]), uym = 0.5F * (uy[ym] + uy[n]);
            const float uzp = 0.5F * (uz[n] + uz[zp]), uzm = 0.5F * (uz[zm] + uz[n]);
            const float fxp = uxp >= 0.0F ? uxp * phi[n] : uxp * phi[xp];
            const float fxm = uxm >= 0.0F ? uxm * phi[xm] : uxm * phi[n];
            const float fyp = uyp >= 0.0F ? uyp * phi[n] : uyp * phi[yp];
            const float fym = uym >= 0.0F ? uym * phi[ym] : uym * phi[n];
            const float fzp = uzp >= 0.0F ? uzp * phi[n] : uzp * phi[zp];
            const float fzm = uzm >= 0.0F ? uzm * phi[zm] : uzm * phi[n];
            float value = phi[n] - dt * ((fxp-fxm) + (fyp-fym) + (fzp-fzm));
            value = value < 0.0F ? 0.0F : value > 1.0F ? 1.0F : value;
            out[n] = value;
        });
        float* old = fill_; fill_ = next_; next_ = old;
    }

    [[nodiscard]] float* fill_device() noexcept { return fill_; }
    [[nodiscard]] const float* fill_device() const noexcept { return fill_; }

    void refresh_geometry() {
        const std::size_t nx = nx_, ny = ny_, nz = nz_, cells = cells_;
        const float* phi = fill_; float* nxv = geometry_; float* nyv = geometry_ + cells;
        float* nzv = geometry_ + 2U*cells; float* curvature = geometry_ + 3U*cells;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n=id[0],x=n%nx,y=(n/nx)%ny,z=n/(nx*ny);
            const auto xp=resident_detail::index(resident_detail::shift(x,1,nx),y,z,nx,ny);
            const auto xm=resident_detail::index(resident_detail::shift(x,-1,nx),y,z,nx,ny);
            const auto yp=resident_detail::index(x,resident_detail::shift(y,1,ny),z,nx,ny);
            const auto ym=resident_detail::index(x,resident_detail::shift(y,-1,ny),z,nx,ny);
            const auto zp=resident_detail::index(x,y,resident_detail::shift(z,1,nz),nx,ny);
            const auto zm=resident_detail::index(x,y,resident_detail::shift(z,-1,nz),nx,ny);
            const float gx=0.5F*(phi[xp]-phi[xm]), gy=0.5F*(phi[yp]-phi[ym]), gz=0.5F*(phi[zp]-phi[zm]);
            const float mag=sycl::sqrt(gx*gx+gy*gy+gz*gz);
            if (mag > 1.0e-12F) { nxv[n]=gx/mag; nyv[n]=gy/mag; nzv[n]=gz/mag; }
            else { nxv[n]=0.0F; nyv[n]=0.0F; nzv[n]=0.0F; }
        });
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n=id[0],x=n%nx,y=(n/nx)%ny,z=n/(nx*ny);
            const auto xp=resident_detail::index(resident_detail::shift(x,1,nx),y,z,nx,ny);
            const auto xm=resident_detail::index(resident_detail::shift(x,-1,nx),y,z,nx,ny);
            const auto yp=resident_detail::index(x,resident_detail::shift(y,1,ny),z,nx,ny);
            const auto ym=resident_detail::index(x,resident_detail::shift(y,-1,ny),z,nx,ny);
            const auto zp=resident_detail::index(x,y,resident_detail::shift(z,1,nz),nx,ny);
            const auto zm=resident_detail::index(x,y,resident_detail::shift(z,-1,nz),nx,ny);
            curvature[n] = -0.5F*((nxv[xp]-nxv[xm])+(nyv[yp]-nyv[ym])+(nzv[zp]-nzv[zm]));
        });
    }

    void add_surface_tension_acceleration(D3Q19SyclSolver::DeviceAccelerationView force,
                                          float surface_tension, float density = 1.0F) {
        if (force.cells != cells_ || !(surface_tension >= 0.0F) || !(density > 0.0F)) {
            throw std::invalid_argument("invalid resident surface-tension coupling");
        }
        refresh_geometry();
        const std::size_t nx=nx_,ny=ny_,nz=nz_,cells=cells_; const float* phi=fill_;
        const float* curvature=geometry_+3U*cells; float* ax=force.ax; float* ay=force.ay; float* az=force.az;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> id) {
            const std::size_t n=id[0],x=n%nx,y=(n/nx)%ny,z=n/(nx*ny);
            const auto xp=resident_detail::index(resident_detail::shift(x,1,nx),y,z,nx,ny);
            const auto xm=resident_detail::index(resident_detail::shift(x,-1,nx),y,z,nx,ny);
            const auto yp=resident_detail::index(x,resident_detail::shift(y,1,ny),z,nx,ny);
            const auto ym=resident_detail::index(x,resident_detail::shift(y,-1,ny),z,nx,ny);
            const auto zp=resident_detail::index(x,y,resident_detail::shift(z,1,nz),nx,ny);
            const auto zm=resident_detail::index(x,y,resident_detail::shift(z,-1,nz),nx,ny);
            const float gx=0.5F*(phi[xp]-phi[xm]), gy=0.5F*(phi[yp]-phi[ym]), gz=0.5F*(phi[zp]-phi[zm]);
            const float c=surface_tension*curvature[n]/density;
            ax[n]+=c*gx; ay[n]+=c*gy; az[n]+=c*gz;
        });
    }

    [[nodiscard]] double volume() const {
        ensure_reduction_scalar(); *reduction_scalar_=0.0F; const float* phi=fill_; const std::size_t cells=cells_;
        auto reduction=sycl::reduction(reduction_scalar_,sycl::plus<float>());
        queue_.parallel_for(sycl::range<1>(cells),reduction,[=](sycl::id<1> id,auto& sum){sum.combine(phi[id[0]]);}).wait_and_throw();
        return static_cast<double>(*reduction_scalar_);
    }
    [[nodiscard]] std::size_t resident_bytes() const noexcept { return 6U*cells_*sizeof(float)+(reduction_scalar_?sizeof(float):0U); }

private:
    sycl::queue& queue_; std::size_t nx_{},ny_{},nz_{},cells_{}; float* fill_{nullptr}; float* next_{nullptr};
    float* geometry_{nullptr}; mutable float* reduction_scalar_{nullptr};
    void ensure_reduction_scalar() const { if(!reduction_scalar_)reduction_scalar_=sycl::malloc_shared<float>(1U,queue_); if(!reduction_scalar_)throw std::bad_alloc{}; }
};

class ResidentParticlesSycl3D {
public:
    ResidentParticlesSycl3D(sycl::queue& queue,std::size_t nx,std::size_t ny,std::size_t nz,std::size_t count)
        :queue_(queue),nx_(nx),ny_(ny),nz_(nz),count_(count){
        if(nx<2U||ny<2U||nz<2U)throw std::invalid_argument("invalid resident particle lattice");
        state_=sycl::malloc_device<float>(8U*count_,queue_); if(count_&&!state_)throw std::bad_alloc{};
    }
    ~ResidentParticlesSycl3D() noexcept { try{queue_.wait_and_throw();}catch(...){} if(state_)sycl::free(state_,queue_); }
    ResidentParticlesSycl3D(const ResidentParticlesSycl3D&)=delete; ResidentParticlesSycl3D& operator=(const ResidentParticlesSycl3D&)=delete;

    void upload(std::span<const ImmersedBoundaryParticle> particles){
        if(particles.size()!=count_)throw std::invalid_argument("resident particle count mismatch");
        std::vector<float> host(8U*count_);
        for(std::size_t i=0;i<count_;++i){const auto&p=particles[i];if(!(p.mass>0.0)||!(p.response_time>0.0))throw std::invalid_argument("particle mass and response time must be positive");
            host[i]=static_cast<float>(p.position.x);host[count_+i]=static_cast<float>(p.position.y);host[2U*count_+i]=static_cast<float>(p.position.z);
            host[3U*count_+i]=static_cast<float>(p.velocity.x);host[4U*count_+i]=static_cast<float>(p.velocity.y);host[5U*count_+i]=static_cast<float>(p.velocity.z);
            host[6U*count_+i]=static_cast<float>(p.mass);host[7U*count_+i]=static_cast<float>(p.response_time);}
        if(count_){queue_.memcpy(state_,host.data(),host.size()*sizeof(float)).wait_and_throw();transfer_stats_.record_host_to_device(host.size()*sizeof(float));transfer_stats_.record_synchronization();}
    }

    void advance(const D3Q19SyclSolver::DeviceMacroscopicView& flow,
                 D3Q19SyclSolver::DeviceAccelerationView force,float dt,bool two_way,
                 ParticleVector3 gravity={}){
        if(flow.cells!=nx_*ny_*nz_||force.cells!=flow.cells||!(dt>0.0F))throw std::invalid_argument("invalid resident particle coupling");
        float* s=state_;const std::size_t count=count_,nx=nx_,ny=ny_,nz=nz_;const float* ux=flow.ux;const float* uy=flow.uy;const float* uz=flow.uz;
        float* ax=force.ax;float* ay=force.ay;float* az=force.az;const float gx=static_cast<float>(gravity.x),gy=static_cast<float>(gravity.y),gz=static_cast<float>(gravity.z);
        queue_.parallel_for(sycl::range<1>(count),[=](sycl::id<1> id){
            const std::size_t p=id[0];float px=resident_detail::periodic_coordinate(s[p],nx),py=resident_detail::periodic_coordinate(s[count+p],ny),pz=resident_detail::periodic_coordinate(s[2U*count+p],nz);
            float vx=s[3U*count+p],vy=s[4U*count+p],vz=s[5U*count+p],mass=s[6U*count+p],response=s[7U*count+p];
            const long long x0=resident_detail::floor_to_ll(px),y0=resident_detail::floor_to_ll(py),z0=resident_detail::floor_to_ll(pz);
            const float tx=px-static_cast<float>(x0),ty=py-static_cast<float>(y0),tz=pz-static_cast<float>(z0);float fu=0.0F,fv=0.0F,fw=0.0F;
            std::size_t ids[8];float weights[8];int k=0;
            for(int dz=0;dz<=1;++dz)for(int dy=0;dy<=1;++dy)for(int dx=0;dx<=1;++dx){const float wx=dx?tx:1.0F-tx,wy=dy?ty:1.0F-ty,wz=dz?tz:1.0F-tz;
                const auto xi=resident_detail::wrap_cell(x0+dx,nx),yi=resident_detail::wrap_cell(y0+dy,ny),zi=resident_detail::wrap_cell(z0+dz,nz);const auto n=resident_detail::index(xi,yi,zi,nx,ny);const float w=wx*wy*wz;ids[k]=n;weights[k]=w;++k;fu+=w*ux[n];fv+=w*uy[n];fw+=w*uz[n];}
            const float dragx=mass*(fu-vx)/response,dragy=mass*(fv-vy)/response,dragz=mass*(fw-vz)/response;
            if(two_way){for(int j=0;j<8;++j){using Atomic=sycl::atomic_ref<float,sycl::memory_order::relaxed,sycl::memory_scope::device,sycl::access::address_space::global_space>;
                    Atomic(ax[ids[j]]).fetch_add(-dragx*weights[j]);Atomic(ay[ids[j]]).fetch_add(-dragy*weights[j]);Atomic(az[ids[j]]).fetch_add(-dragz*weights[j]);}}
            vx+=dt*(dragx/mass+gx);vy+=dt*(dragy/mass+gy);vz+=dt*(dragz/mass+gz);px=resident_detail::periodic_coordinate(px+dt*vx,nx);py=resident_detail::periodic_coordinate(py+dt*vy,ny);pz=resident_detail::periodic_coordinate(pz+dt*vz,nz);
            s[p]=px;s[count+p]=py;s[2U*count+p]=pz;s[3U*count+p]=vx;s[4U*count+p]=vy;s[5U*count+p]=vz;
        });
    }

    [[nodiscard]] std::vector<ImmersedBoundaryParticle> download() const {
        std::vector<float> host(8U*count_);if(count_)queue_.memcpy(host.data(),state_,host.size()*sizeof(float)).wait_and_throw();
        transfer_stats_.record_device_to_host(host.size()*sizeof(float));transfer_stats_.record_synchronization();std::vector<ImmersedBoundaryParticle> out(count_);
        for (std::size_t i = 0; i < count_; ++i) {
            out[i] = {{host[i], host[count_ + i], host[2U * count_ + i]},
                      {host[3U * count_ + i], host[4U * count_ + i], host[5U * count_ + i]},
                      host[6U * count_ + i], host[7U * count_ + i]};
        }
        return out;
    }
    [[nodiscard]] std::size_t resident_bytes()const noexcept{return 8U*count_*sizeof(float);} [[nodiscard]] std::size_t size()const noexcept{return count_;}
    [[nodiscard]] const cfd::core::DeviceTransferStats& transfer_stats()const noexcept{return transfer_stats_;}void reset_transfer_stats()const noexcept{transfer_stats_.reset();}
private: sycl::queue&queue_;std::size_t nx_{},ny_{},nz_{},count_{};float*state_{nullptr};mutable cfd::core::DeviceTransferStats transfer_stats_{};
};

class ResidentFlowDiagnosticsSycl3D {
public:
    ResidentFlowDiagnosticsSycl3D(sycl::queue&queue,std::size_t nx,std::size_t ny,std::size_t nz):queue_(queue),nx_(nx),ny_(ny),nz_(nz),cells_(nx*ny*nz){
        if (nx < 3U || ny < 3U || nz < 3U) {
            throw std::invalid_argument("Q-criterion lattice must be at least 3x3x3");
        }
        q_ = sycl::malloc_device<float>(cells_, queue_);
        if (!q_) throw std::bad_alloc{};
    }
    ~ResidentFlowDiagnosticsSycl3D()noexcept{try{queue_.wait_and_throw();}catch(...){}if(q_)sycl::free(q_,queue_);} ResidentFlowDiagnosticsSycl3D(const ResidentFlowDiagnosticsSycl3D&)=delete;ResidentFlowDiagnosticsSycl3D&operator=(const ResidentFlowDiagnosticsSycl3D&)=delete;
    [[nodiscard]] float* q_criterion(const D3Q19SyclSolver::DeviceMacroscopicView&flow){if(flow.cells!=cells_)throw std::invalid_argument("diagnostic/flow lattice size mismatch");
        const std::size_t nx=nx_,ny=ny_,nz=nz_,cells=cells_;const float*ux=flow.ux,*uy=flow.uy,*uz=flow.uz;float*out=q_;
        queue_.parallel_for(sycl::range<1>(cells),[=](sycl::id<1>id){const std::size_t n=id[0],x=n%nx,y=(n/nx)%ny,z=n/(nx*ny);float g[3][3]{};
            for(int axis=0;axis<3;++axis){std::size_t xp=x,xm=x,yp=y,ym=y,zp=z,zm=z;if(axis==0){xp=resident_detail::shift(x,1,nx);xm=resident_detail::shift(x,-1,nx);}else if(axis==1){yp=resident_detail::shift(y,1,ny);ym=resident_detail::shift(y,-1,ny);}else{zp=resident_detail::shift(z,1,nz);zm=resident_detail::shift(z,-1,nz);}const auto p=resident_detail::index(xp,yp,zp,nx,ny),m=resident_detail::index(xm,ym,zm,nx,ny);g[0][axis]=0.5F*(ux[p]-ux[m]);g[1][axis]=0.5F*(uy[p]-uy[m]);g[2][axis]=0.5F*(uz[p]-uz[m]);}
            float s2=0.0F,o2=0.0F;for(int i=0;i<3;++i)for(int j=0;j<3;++j){const float sij=0.5F*(g[i][j]+g[j][i]),oij=0.5F*(g[i][j]-g[j][i]);s2+=sij*sij;o2+=oij*oij;}out[n]=0.5F*(o2-s2);});return q_;}
    [[nodiscard]] std::size_t resident_bytes()const noexcept{return cells_*sizeof(float);}private:sycl::queue&queue_;std::size_t nx_{},ny_{},nz_{},cells_{};float*q_{nullptr};
};

} // namespace cfd::lbm
#endif
