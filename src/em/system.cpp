#include "cfd/em/system.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>

namespace cfd::em {
namespace {
void validate_frequency(double f){if(!(f>=0.0)||!std::isfinite(f))throw std::invalid_argument("invalid frequency");}

double mean(std::span<const double> values){double s=0.0;for(double v:values)s+=v;return s/static_cast<double>(values.size());}

double standard_deviation(std::span<const double> values,double m){if(values.size()<2U)return 0.0;double s=0.0;for(double v:values){const double d=v-m;s+=d*d;}return std::sqrt(s/static_cast<double>(values.size()-1U));}

double qfunc(double x){return 0.5*std::erfc(x/std::sqrt(2.0));}

double worst_ratio(std::span<const double> frequencies,std::span<const double> targets,
                   std::span<const DecouplingCapacitor> caps,Complex plane){
    double worst=0.0;
    for(std::size_t i=0;i<frequencies.size();++i){
        const double target=targets[i];
        if(!(target>0.0)||!std::isfinite(target))throw std::invalid_argument("invalid PDN target impedance");
        worst=std::max(worst,std::abs(pdn_parallel_impedance_ohm(caps,frequencies[i],plane))/target);
    }
    return worst;
}
}

EyeDiagramMetrics analyze_nrz_eye(std::span<const double> waveform,std::span<const int> symbols,
                                  std::size_t samples_per_symbol,std::size_t decision_sample){
    if(samples_per_symbol<2U||symbols.empty())throw std::invalid_argument("invalid eye-diagram controls");
    if(decision_sample==0U)decision_sample=samples_per_symbol/2U;
    if(decision_sample>=samples_per_symbol)throw std::invalid_argument("invalid eye decision sample");
    if(waveform.size()<symbols.size()*samples_per_symbol)throw std::invalid_argument("waveform shorter than bit stream");
    std::vector<double> zeros,ones;zeros.reserve(symbols.size());ones.reserve(symbols.size());
    for(std::size_t k=0;k<symbols.size();++k){
        const double v=waveform[k*samples_per_symbol+decision_sample];
        if(!std::isfinite(v))throw std::invalid_argument("non-finite eye waveform");
        if(symbols[k]==0)zeros.push_back(v);else if(symbols[k]==1)ones.push_back(v);else throw std::invalid_argument("eye symbols must be 0/1");
    }
    if(zeros.empty()||ones.empty())throw std::invalid_argument("eye analysis requires both zero and one symbols");
    EyeDiagramMetrics result;result.symbol_count=symbols.size();result.samples_per_symbol=samples_per_symbol;
    result.zero_mean=mean(zeros);result.one_mean=mean(ones);result.zero_sigma=standard_deviation(zeros,result.zero_mean);result.one_sigma=standard_deviation(ones,result.one_mean);
    result.decision_threshold=0.5*(result.zero_mean+result.one_mean);
    result.eye_height=(result.one_mean-3.0*result.one_sigma)-(result.zero_mean+3.0*result.zero_sigma);
    const double sigma0=std::max(result.zero_sigma,1.0e-30),sigma1=std::max(result.one_sigma,1.0e-30);
    result.estimated_ber=0.5*(qfunc((result.one_mean-result.decision_threshold)/sigma1)+qfunc((result.decision_threshold-result.zero_mean)/sigma0));
    std::size_t open_samples=0U;
    for(std::size_t phase=0;phase<samples_per_symbol;++phase){
        double max_zero=-std::numeric_limits<double>::infinity();
        double min_one=std::numeric_limits<double>::infinity();
        for(std::size_t k=0;k<symbols.size();++k){
            const double v=waveform[k*samples_per_symbol+phase];
            if(symbols[k]==0)max_zero=std::max(max_zero,v);else min_one=std::min(min_one,v);
        }
        if(min_one>max_zero)++open_samples;
    }
    result.eye_width_ui=static_cast<double>(open_samples)/static_cast<double>(samples_per_symbol);
    return result;
}

Complex decoupling_impedance_ohm(const DecouplingCapacitor& capacitor,double frequency){
    validate_frequency(frequency);
    if(!(capacitor.capacitance_f>0.0)||!(capacitor.esr_ohm>=0.0)||!(capacitor.esl_h>=0.0)||capacitor.count==0U)
        throw std::invalid_argument("invalid decoupling capacitor");
    const double omega=2.0*std::numbers::pi*frequency;
    Complex z=Complex{capacitor.esr_ohm,omega*capacitor.esl_h};
    if(omega==0.0)z+=Complex{0.0,-std::numeric_limits<double>::infinity()};
    else z+=Complex{0.0,-1.0/(omega*capacitor.capacitance_f)};
    return z/static_cast<double>(capacitor.count);
}

Complex pdn_parallel_impedance_ohm(std::span<const DecouplingCapacitor> capacitors,double frequency,Complex plane){
    validate_frequency(frequency);Complex admittance{};
    if(std::abs(plane)>0.0)admittance+=1.0/plane;
    for(const auto& cap:capacitors){const Complex z=decoupling_impedance_ohm(cap,frequency);if(std::abs(z)>0.0&&std::isfinite(z.real())&&std::isfinite(z.imag()))admittance+=1.0/z;}
    if(std::abs(admittance)==0.0)return {std::numeric_limits<double>::infinity(),0.0};
    return 1.0/admittance;
}

DecouplingOptimizationResult optimize_decoupling_greedy(std::span<const double> frequencies,std::span<const double> targets,
                                                        std::span<const DecouplingCapacitor> library,
                                                        std::size_t maximum_parts,Complex plane){
    if(frequencies.empty()||frequencies.size()!=targets.size()||library.empty())throw std::invalid_argument("invalid decoupling optimization inputs");
    DecouplingOptimizationResult result;result.worst_ratio_before=worst_ratio(frequencies,targets,{},plane);result.worst_ratio_after=result.worst_ratio_before;
    for(std::size_t iteration=0;iteration<maximum_parts;++iteration){
        double best=result.worst_ratio_after;std::size_t best_index=library.size();std::vector<DecouplingCapacitor> best_selection;
        for(std::size_t candidate=0;candidate<library.size();++candidate){
            std::vector<DecouplingCapacitor> trial=result.selected;trial.push_back(library[candidate]);
            const double ratio=worst_ratio(frequencies,targets,trial,plane);
            if(ratio<best){best=ratio;best_index=candidate;best_selection=std::move(trial);}
        }
        if(best_index==library.size())break;
        result.selected=std::move(best_selection);result.worst_ratio_after=best;++result.iterations;
        if(result.worst_ratio_after<=1.0)break;
    }
    return result;
}

PdnIrDropResult solve_pdn_ir_drop_grid(const PdnGridConfig& config,std::span<const double> loads,
                                       std::span<const PdnVoltageSource> sources){
    if(config.nx<2U||config.ny<2U||!(config.horizontal_resistance_ohm>0.0)||!(config.vertical_resistance_ohm>0.0)||config.max_iterations==0U||!(config.relative_tolerance>0.0))throw std::invalid_argument("invalid PDN grid");
    const std::size_t n=config.nx*config.ny;if(loads.size()!=n||sources.empty())throw std::invalid_argument("PDN load/source size mismatch");
    std::vector<double> fixed(n,std::numeric_limits<double>::quiet_NaN()),voltage(n,0.0);
    for(const auto& source:sources){if(source.node>=n||!std::isfinite(source.voltage_v))throw std::invalid_argument("invalid PDN voltage source");fixed[source.node]=source.voltage_v;voltage[source.node]=source.voltage_v;}
    double source_mean=0.0;std::size_t source_count=0U;for(double v:fixed)if(!std::isnan(v)){source_mean+=v;++source_count;}source_mean/=static_cast<double>(source_count);for(std::size_t i=0;i<n;++i)if(std::isnan(fixed[i]))voltage[i]=source_mean;
    const double gx=1.0/config.horizontal_resistance_ohm,gy=1.0/config.vertical_resistance_ohm;
    double residual0=0.0;
    PdnIrDropResult result;result.converged=false;
    for(std::size_t iteration=0;iteration<config.max_iterations;++iteration){
        double residual=0.0;
        for(std::size_t y=0;y<config.ny;++y){
            for(std::size_t x=0;x<config.nx;++x){
                const std::size_t node=y*config.nx+x;if(!std::isnan(fixed[node])){voltage[node]=fixed[node];continue;}
                double conductance=0.0,rhs=-loads[node];
                if(x>0U){conductance+=gx;rhs+=gx*voltage[node-1U];}
                if(x+1U<config.nx){conductance+=gx;rhs+=gx*voltage[node+1U];}
                if(y>0U){conductance+=gy;rhs+=gy*voltage[node-config.nx];}
                if(y+1U<config.ny){conductance+=gy;rhs+=gy*voltage[node+config.nx];}
                const double next=rhs/conductance;residual=std::max(residual,std::abs(next-voltage[node]));voltage[node]=next;
            }
        }
        if(iteration==0U)residual0=std::max(residual,1.0e-30);
        result.iterations=iteration+1U;
        if(residual<=config.relative_tolerance*residual0){result.converged=true;break;}
    }
    result.voltage_v=std::move(voltage);result.minimum_voltage_v=*std::min_element(result.voltage_v.begin(),result.voltage_v.end());
    double maximum_source=*std::max_element(fixed.begin(),fixed.end(),[](double a,double b){return (std::isnan(a)?-std::numeric_limits<double>::infinity():a)<(std::isnan(b)?-std::numeric_limits<double>::infinity():b);});
    result.maximum_drop_v=maximum_source-result.minimum_voltage_v;return result;
}

double normalized_double_exponential_pulse(double time,double amplitude,double rise_tau,double fall_tau){
    if(time<0.0)return 0.0;
    if(!(amplitude>=0.0)||!(rise_tau>0.0)||!(fall_tau>rise_tau))throw std::invalid_argument("invalid double exponential pulse");
    const double peak_time=rise_tau*fall_tau*std::log(fall_tau/rise_tau)/(fall_tau-rise_tau);
    const double peak=std::exp(-peak_time/fall_tau)-std::exp(-peak_time/rise_tau);
    return amplitude*(std::exp(-time/fall_tau)-std::exp(-time/rise_tau))/peak;
}

std::vector<double> sample_double_exponential_pulse(double duration,double dt,double amplitude,double rise_tau,double fall_tau){
    if(!(duration>0.0)||!(dt>0.0))throw std::invalid_argument("invalid waveform sampling");
    const std::size_t n=static_cast<std::size_t>(std::floor(duration/dt))+1U;std::vector<double> values(n);
    for(std::size_t i=0;i<n;++i)values[i]=normalized_double_exponential_pulse(static_cast<double>(i)*dt,amplitude,rise_tau,fall_tau);
    return values;
}

std::vector<double> sample_damped_sine(double duration,double dt,double amplitude,double frequency,double decay_tau){
    if(!(duration>0.0)||!(dt>0.0)||!(frequency>=0.0)||!(decay_tau>0.0))throw std::invalid_argument("invalid damped sine sampling");
    const std::size_t n=static_cast<std::size_t>(std::floor(duration/dt))+1U;std::vector<double> values(n);
    for(std::size_t i=0;i<n;++i){const double t=static_cast<double>(i)*dt;values[i]=amplitude*std::exp(-t/decay_tau)*std::sin(2.0*std::numbers::pi*frequency*t);}return values;
}

ProbeMetrics probe_waveform(std::span<const double> samples,double dt){
    if(samples.empty()||!(dt>0.0))throw std::invalid_argument("invalid probe waveform");
    ProbeMetrics out;
    double sum2=0.0;
    for(double v:samples){if(!std::isfinite(v))throw std::invalid_argument("non-finite waveform sample");out.peak=std::max(out.peak,std::abs(v));sum2+=v*v;out.impulse+=v*dt;out.energy+=v*v*dt;}
    out.rms=std::sqrt(sum2/static_cast<double>(samples.size()));return out;
}

} // namespace cfd::em
