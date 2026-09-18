#include "cfd/battery/thermal3d.hpp"
#include "cfd/core/conjugate_gradient.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace cfd::battery {
namespace {
bool positive(double x) { return std::isfinite(x) && x>0; }
bool nonnegative(double x) { return std::isfinite(x) && x>=0; }
}
ThermalGrid3D::ThermalGrid3D(ThermalGridConfig c, std::vector<ThermalVoxel> m):config_(c) {
    std::size_t n=1;
    double volume=1;
    for(std::size_t a=0;a<3;++a) {
        if(c.cells[a]==0 || c.cells[a]>std::numeric_limits<std::size_t>::max()/n || !positive(c.spacing_m[a]))
            throw std::invalid_argument("invalid thermal grid dimensions");
        n*=c.cells[a]; volume*=c.spacing_m[a];
    }
    if(!positive(volume) || !positive(c.initial_temperature_k) || c.maximum_iterations==0
       || !positive(c.relative_tolerance) || c.relative_tolerance>=1)
        throw std::invalid_argument("invalid thermal solver configuration");
    for(const auto& b:c.boundaries)
        if(!nonnegative(b.convection_w_per_m2_k) || !positive(b.ambient_temperature_k))
            throw std::invalid_argument("invalid thermal boundary");
    if(m.size()==1) m.resize(n,m.front());
    if(m.size()!=n) throw std::invalid_argument("thermal material count mismatch");
    temperature_.assign(n,c.initial_temperature_k); capacity_.resize(n);
    cooling_.assign(n,0); ambient_source_.assign(n,0);
    for(std::size_t i=0;i<n;++i) {
        capacity_[i]=volume*m[i].volumetric_heat_capacity_j_per_m3_k;
        if(!positive(capacity_[i])) throw std::invalid_argument("invalid voxel heat capacity");
        for(double k:m[i].conductivity_w_per_m_k)
            if(!nonnegative(k)) throw std::invalid_argument("invalid voxel conductivity");
    }
    const std::array<std::size_t,3> stride{1,c.cells[0],c.cells[0]*c.cells[1]};
    for(std::size_t i=0;i<n;++i) for(std::size_t a=0;a<3;++a) {
        const auto coordinate=(i/stride[a])%c.cells[a];
        const double area=volume/c.spacing_m[a], k=m[i].conductivity_w_per_m_k[a];
        if(coordinate+1<c.cells[a]) {
            const auto j=i+stride[a]; const double kj=m[j].conductivity_w_per_m_k[a];
            const double g=(k>0 && kj>0)?area/(0.5*c.spacing_m[a]/k+0.5*c.spacing_m[a]/kj):0;
            if(!nonnegative(g)) throw std::invalid_argument("thermal face conductance overflow");
            faces_.push_back({i,j,g});
        }
        for(std::size_t side=0;side<2;++side) if(side==0?coordinate==0:coordinate+1==c.cells[a]) {
            const auto b=c.boundaries[2*a+side]; const double h=b.convection_w_per_m2_k;
            const double g=(h>0 && k>0)?area/(0.5*c.spacing_m[a]/k+1/h):0;
            cooling_[i]+=g; ambient_source_[i]+=g*b.ambient_temperature_k;
        }
        if(!nonnegative(cooling_[i]) || !nonnegative(ambient_source_[i]))
            throw std::invalid_argument("thermal boundary conductance overflow");
    }
}
void ThermalGrid3D::reset(std::span<const double> t) {
    if(t.size()!=temperature_.size() || !std::all_of(t.begin(),t.end(),positive))
        throw std::invalid_argument("invalid thermal reset field");
    temperature_.assign(t.begin(),t.end());
}
ThermalStepResult ThermalGrid3D::step(double dt,std::span<const double> q) {
    const auto n=temperature_.size();
    if(!positive(dt) || q.size()!=n) throw std::invalid_argument("invalid thermal step");
    auto conduction=[&](std::span<const double> t,std::span<double> out) {
        for(std::size_t i=0;i<n;++i) out[i]=cooling_[i]*t[i];
        for(const auto& f:faces_) {
            const double flux=f.conductance*(t[f.a]-t[f.b]);
            out[f.a]+=flux; out[f.b]-=flux;
        }
    };
    std::vector<double> rhs(n), delta(n,0), loss(n);
    conduction(temperature_,loss);
    for(std::size_t i=0;i<n;++i) {
        rhs[i]=dt*(q[i]+ambient_source_[i]-loss[i]);
        if(!std::isfinite(q[i]) || !std::isfinite(rhs[i])) throw std::invalid_argument("nonfinite thermal power");
    }
    auto apply=[&](std::span<const double> x,std::span<double> y) {
        conduction(x,y);
        for(std::size_t i=0;i<n;++i) y[i]=capacity_[i]*x[i]+dt*y[i];
    };
    core::ConjugateGradientWorkspace workspace;
    const auto solved=core::conjugate_gradient(rhs,delta,apply,workspace,config_.maximum_iterations,config_.relative_tolerance);
    if(!solved.converged) throw std::runtime_error("3D thermal CG did not converge");
    auto next=temperature_; ThermalStepResult out; out.iterations=solved.iterations;
    for(std::size_t i=0;i<n;++i) {
        next[i]+=delta[i];
        if(!positive(next[i])) throw std::domain_error("invalid thermal temperature");
        out.supplied_power_w+=q[i];
        out.outward_cooling_power_w+=cooling_[i]*next[i]-ambient_source_[i];
        out.stored_energy_change_j+=capacity_[i]*delta[i];
        out.maximum_temperature_k=std::max(out.maximum_temperature_k,next[i]);
    }
    out.energy_balance_error_j=out.stored_energy_change_j-dt*(out.supplied_power_w-out.outward_cooling_power_w);
    if(!std::isfinite(out.energy_balance_error_j)) throw std::overflow_error("thermal energy overflow");
    temperature_=std::move(next);
    return out;
}
ElectrothermalBatteryPack::ElectrothermalBatteryPack(std::span<const LithiumIonCellConfig> cells,PackControl control,
    ThermalGrid3D thermal,std::vector<std::vector<std::size_t>> regions)
    :pack_(cells,control),thermal_(std::move(thermal)),regions_(std::move(regions)),maximum_temperature_k_(control.maximum_temperature_k) {
    if(regions_.size()!=cells.size()) throw std::invalid_argument("one thermal region required per cell");
    std::vector<bool> used(thermal_.temperature().size(),false);
    for(std::size_t c=0;c<regions_.size();++c) {
        if(regions_[c].empty()) throw std::invalid_argument("empty cell thermal region");
        double capacity=0;
        for(auto i:regions_[c]) {
            if(i>=used.size() || used[i]) throw std::invalid_argument("overlapping or invalid cell thermal voxel");
            used[i]=true; capacity+=thermal_.heat_capacity_j_per_k()[i];
        }
        const double expected=cells[c].thermal.mass_kg*cells[c].thermal.heat_capacity_j_per_kg_k;
        if(!positive(expected) || !positive(capacity) || std::abs(capacity-expected)>1e-8*expected)
            throw std::invalid_argument("cell and mapped thermal heat capacities must agree");
    }
    if(*std::max_element(thermal_.temperature().begin(),thermal_.temperature().end())>maximum_temperature_k_)
        throw std::domain_error("initial pack thermal hotspot cutoff");
    pack_.set_external_temperatures(cell_temperatures(thermal_));
}
std::vector<double> ElectrothermalBatteryPack::cell_temperatures(const ThermalGrid3D& field) const {
    std::vector<double> t;
    for(const auto& region:regions_) {
        double energy=0,capacity=0;
        for(auto i:region) { const double c=field.heat_capacity_j_per_k()[i]; energy+=c*field.temperature()[i]; capacity+=c; }
        t.push_back(energy/capacity);
    }
    return t;
}
ElectrothermalPackResult ElectrothermalBatteryPack::step(double current,double dt) {
    auto trial_pack=pack_; auto trial_thermal=thermal_;
    auto electrical=trial_pack.step(current,dt);
    std::vector<double> q(thermal_.temperature().size(),0);
    for(std::size_t c=0;c<regions_.size();++c) {
        // Uniform volumetric deposition: this regular grid has equal voxel volumes.
        const double power=electrical.cells[c].heat_generation_w/static_cast<double>(regions_[c].size());
        for(auto i:regions_[c]) q[i]+=power;
    }
    const auto thermal=trial_thermal.step(dt,q);
    if(thermal.maximum_temperature_k>maximum_temperature_k_) throw std::domain_error("pack thermal hotspot cutoff");
    const auto temperatures=cell_temperatures(trial_thermal);
    trial_pack.set_external_temperatures(temperatures);
    for(std::size_t c=0;c<temperatures.size();++c) electrical.cells[c].temperature_k=temperatures[c];
    pack_=std::move(trial_pack); thermal_=std::move(trial_thermal);
    return {std::move(electrical),thermal};
}
}
