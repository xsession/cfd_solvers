#include "cfd/circuit/analysis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>

namespace cfd::circuit {

cfd::rf::Matrix2C two_port_s_parameters(const Circuit& circuit,const AcPort& first,const AcPort& second,
                                         double frequency,double reference_impedance) {
    if (!(frequency > 0.0) || !std::isfinite(frequency) || !(reference_impedance > 0.0))
        throw std::invalid_argument("invalid two-port S-parameter controls");
    const auto excite=[&](const AcPort& driven,const AcPort& measured1,const AcPort& measured2){
        Circuit copy=circuit;
        const auto dp=copy.find_node(driven.positive),dn=copy.find_node(driven.negative);
        copy.add_current_source("__PORT_EXCITATION__",dp,dn,0.0,Complex{-1.0,0.0});
        const auto ac=copy.ac(frequency);
        const auto voltage=[&](const AcPort& port){
            const auto p=copy.find_node(port.positive),n=copy.find_node(port.negative);
            return ac.node_voltage[p]-ac.node_voltage[n];
        };
        return std::array<Complex,2>{voltage(measured1),voltage(measured2)};
    };
    const auto column1=excite(first,first,second);
    const auto column2=excite(second,first,second);
    return cfd::rf::z_to_s({column1[0],column2[0],column1[1],column2[1]},reference_impedance);
}

SensitivityResult voltage_source_dc_sensitivity(Circuit circuit,std::string_view source,
                                                std::string_view output,double perturbation,
                                                const NewtonConfig& config) {
    if (!(perturbation>0.0)||!std::isfinite(perturbation)) throw std::invalid_argument("sensitivity perturbation must be positive");
    const double nominal=circuit.voltage_source_dc(source);
    const auto base=circuit.dc_operating_point_homotopy(config);
    if(!base.converged)throw std::runtime_error("sensitivity baseline operating point failed");
    circuit.set_voltage_source_dc(source,nominal+perturbation);const auto plus=circuit.dc_operating_point_homotopy(config);
    circuit.set_voltage_source_dc(source,nominal-perturbation);const auto minus=circuit.dc_operating_point_homotopy(config);
    if(!plus.converged||!minus.converged)throw std::runtime_error("sensitivity perturbed operating point failed");
    return {circuit.voltage(base,output),(circuit.voltage(plus,output)-circuit.voltage(minus,output))/(2.0*perturbation)};
}

FourierMeasurement fourier_measurement(std::span<const TransientPoint> samples,std::size_t node,double frequency) {
    if(samples.size()<2U||!(frequency>=0.0)||!std::isfinite(frequency))throw std::invalid_argument("invalid Fourier measurement");
    Complex integral{};double duration=samples.back().time-samples.front().time;
    if(!(duration>0.0))throw std::invalid_argument("transient samples need increasing time");
    for(std::size_t i=1;i<samples.size();++i){
        if(node>=samples[i-1].node_voltage.size()||node>=samples[i].node_voltage.size())throw std::out_of_range("Fourier node index out of range");
        const double t0=samples[i-1].time,t1=samples[i].time,dt=t1-t0;if(!(dt>0.0))throw std::invalid_argument("transient sample times must increase");
        const Complex e0=std::exp(Complex{0.0,-2.0*std::numbers::pi*frequency*t0});
        const Complex e1=std::exp(Complex{0.0,-2.0*std::numbers::pi*frequency*t1});
        integral+=0.5*dt*(samples[i-1].node_voltage[node]*e0+samples[i].node_voltage[node]*e1);
    }
    const Complex phasor=(frequency==0.0?1.0:2.0)*integral/duration;
    return {frequency,phasor,std::abs(phasor),std::arg(phasor)};
}

double total_harmonic_distortion(std::span<const TransientPoint> samples,std::size_t node,double fundamental,std::size_t harmonics) {
    if(!(fundamental>0.0)||harmonics<2U)throw std::invalid_argument("invalid THD controls");
    const double base=fourier_measurement(samples,node,fundamental).magnitude;if(!(base>0.0))return std::numeric_limits<double>::infinity();
    double sum=0.0;for(std::size_t h=2;h<=harmonics;++h){const double amplitude=fourier_measurement(samples,node,fundamental*static_cast<double>(h)).magnitude;sum+=amplitude*amplitude;}
    return std::sqrt(sum)/base;
}

MonteCarloResult monte_carlo_dc_resistors(const Circuit& nominal,std::string_view output,
                                           std::span<const ResistorTolerance> tolerances,
                                           std::size_t samples,std::uint64_t seed,const NewtonConfig& config) {
    if(samples==0U)throw std::invalid_argument("Monte Carlo requires at least one sample");
    for(const auto&t:tolerances)if(!(t.relative_sigma>=0.0)||!std::isfinite(t.relative_sigma))throw std::invalid_argument("invalid resistor tolerance");
    std::mt19937_64 rng(seed);std::normal_distribution<double> normal;
    MonteCarloResult result;result.samples=samples;result.minimum=std::numeric_limits<double>::infinity();result.maximum=-result.minimum;double mean=0.0,m2=0.0;
    for(std::size_t sample=0;sample<samples;++sample){Circuit circuit=nominal;for(const auto&t:tolerances){const double r=nominal.resistance(t.element);double varied=r*(1.0+t.relative_sigma*normal(rng));if(!(varied>0.0))varied=std::max(r*1.0e-6,1.0e-18);circuit.set_resistance(t.element,varied);}const auto op=circuit.dc_operating_point_homotopy(config);if(!op.converged)throw std::runtime_error("Monte Carlo operating point failed");const double value=circuit.voltage(op,output);const double delta=value-mean;mean+=delta/static_cast<double>(sample+1U);m2+=delta*(value-mean);result.minimum=std::min(result.minimum,value);result.maximum=std::max(result.maximum,value);}
    result.mean=mean;result.standard_deviation=std::sqrt(std::max(0.0,m2/static_cast<double>(samples)));return result;
}


cfd::rf::ComplexMatrix extract_s_parameters(Circuit circuit,std::span<const CircuitPort> ports,
                                              double frequency,const NewtonConfig& config) {
    if(ports.empty()||!(frequency>0.0)||!std::isfinite(frequency))
        throw std::invalid_argument("invalid circuit S-parameter extraction controls");
    std::vector<double> reference(ports.size());
    for(std::size_t i=0;i<ports.size();++i){
        if(ports[i].positive>=circuit.node_count()||ports[i].negative>=circuit.node_count()
           ||!(ports[i].reference_impedance>0.0)||!std::isfinite(ports[i].reference_impedance))
            throw std::invalid_argument("invalid circuit port");
        reference[i]=ports[i].reference_impedance;
    }
    cfd::rf::ComplexMatrix y(ports.size());
    for(std::size_t excitation=0;excitation<ports.size();++excitation){
        Circuit driven=circuit;
        std::vector<std::string> names(ports.size());
        for(std::size_t p=0;p<ports.size();++p){
            names[p]="__CFD_PORT_"+std::to_string(p);
            driven.add_voltage_source(names[p],ports[p].positive,ports[p].negative,0.0,
                                      p==excitation?Complex{1.0,0.0}:Complex{});
        }
        const auto operating=driven.dc_operating_point_homotopy(config);
        if(!operating.converged)throw std::runtime_error("circuit S-parameter operating point failed");
        const auto solution=driven.ac(frequency,&operating,config);
        for(std::size_t response=0;response<ports.size();++response){
            // MNA source current is positive from + to - through the ideal source;
            // current entering the DUT is its negative.
            y(response,excitation)=-driven.branch_current(solution,names[response]);
        }
    }
    const auto z=cfd::rf::inverse(y);
    return cfd::rf::z_to_s(z,reference);
}

double limit_pn_junction_voltage(double proposed,double previous,double thermal_voltage,double saturation_current) {
    if(!(thermal_voltage>0.0)||!(saturation_current>0.0)||!std::isfinite(proposed)||!std::isfinite(previous))throw std::invalid_argument("invalid PN junction limiter input");
    const double vcrit=thermal_voltage*std::log(thermal_voltage/(std::sqrt(2.0)*saturation_current));
    if(proposed>vcrit&&std::abs(proposed-previous)>2.0*thermal_voltage){
        if(previous>0.0){const double argument=1.0+(proposed-previous)/thermal_voltage;return argument>0.0?previous+thermal_voltage*std::log(argument):vcrit;}
        return thermal_voltage*std::log(std::max(proposed/thermal_voltage,1.0));
    }
    return proposed;
}


std::vector<ParameterSweepPoint> parameter_sweep_dc(const Circuit& nominal,std::span<const double> values,
                                                     const CircuitParameterSetter& setter,const NewtonConfig& config) {
    if(values.empty()||!setter) throw std::invalid_argument("invalid parameter sweep controls");
    std::vector<ParameterSweepPoint> result;
    result.reserve(values.size());
    for(double value:values){
        if(!std::isfinite(value)) throw std::invalid_argument("non-finite sweep value");
        Circuit circuit=nominal;
        setter(circuit,value);
        auto op=circuit.dc_operating_point_homotopy(config);
        if(!op.converged) throw std::runtime_error("parameter sweep operating point failed");
        result.push_back({value,std::move(op)});
    }
    return result;
}

std::vector<ParameterSweepPoint> temperature_sweep_dc(const Circuit& nominal,
                                                       std::span<const double> temperatures,
                                                       const NewtonConfig& config) {
    return parameter_sweep_dc(nominal,temperatures,
        [](Circuit& circuit,double temperature){ circuit.set_device_temperature(temperature); },config);
}

NoiseSpectrum output_noise_spectrum(const Circuit& circuit,std::string_view output_node,
                                    double start,double stop,std::size_t count,double temperature,
                                    const NewtonConfig& config) {
    if(!(start>0.0)||!(stop>=start)||count==0U||!(temperature>0.0)
       ||!std::isfinite(start)||!std::isfinite(stop)||!std::isfinite(temperature))
        throw std::invalid_argument("invalid noise spectrum controls");
    NoiseSpectrum result;
    result.points.reserve(count);
    if(count==1U){result.points.push_back(circuit.output_noise(start,output_node,temperature,config));return result;}
    const double ratio=std::pow(stop/start,1.0/static_cast<double>(count-1U));
    for(std::size_t i=0;i<count;++i){
        const double frequency=i+1U==count?stop:start*std::pow(ratio,static_cast<double>(i));
        result.points.push_back(circuit.output_noise(frequency,output_node,temperature,config));
    }
    double variance=0.0;
    for(std::size_t i=1;i<result.points.size();++i){
        const auto& a=result.points[i-1U];const auto& b=result.points[i];
        variance+=0.5*(a.output_noise_density_v2_per_hz+b.output_noise_density_v2_per_hz)
            *(b.frequency_hz-a.frequency_hz);
    }
    result.integrated_rms_voltage=std::sqrt(std::max(0.0,variance));
    return result;
}


DominantPoleEstimate estimate_dominant_pole(const Circuit& circuit,std::string_view output_node,
                                            double start,double stop,std::size_t count,
                                            const NewtonConfig& config) {
    if(!(start>0.0)||!(stop>start)||count<3U) throw std::invalid_argument("invalid pole-estimation controls");
    const auto sweep=circuit.ac_log_sweep(start,stop,count,config);
    const double low=std::abs(circuit.voltage(sweep.front(),output_node));
    if(!(low>0.0)||!std::isfinite(low)) return {low,0.0,false};
    const double target=low/std::sqrt(2.0);
    for(std::size_t i=1;i<sweep.size();++i){
        const double previous=std::abs(circuit.voltage(sweep[i-1U],output_node));
        const double current=std::abs(circuit.voltage(sweep[i],output_node));
        if(previous>=target&&current<=target&&previous!=current){
            const double x0=std::log(sweep[i-1U].frequency_hz),x1=std::log(sweep[i].frequency_hz);
            const double fraction=std::clamp((target-previous)/(current-previous),0.0,1.0);
            return {low,std::exp(x0+fraction*(x1-x0)),true};
        }
    }
    return {low,0.0,false};
}

} // namespace cfd::circuit
