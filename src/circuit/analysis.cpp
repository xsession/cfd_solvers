#include "cfd/circuit/analysis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>

namespace cfd::circuit {
namespace {

template<class T>
std::vector<T> solve_dense_ls(std::vector<T> matrix,std::vector<T> rhs,std::size_t n) {
    for(std::size_t column=0;column<n;++column){
        std::size_t pivot=column;double best=std::abs(matrix[column*n+column]);
        for(std::size_t row=column+1U;row<n;++row){const double candidate=std::abs(matrix[row*n+column]);if(candidate>best){best=candidate;pivot=row;}}
        if(best<1.0e-24)throw std::runtime_error("singular rational-fit normal matrix");
        if(pivot!=column){for(std::size_t j=column;j<n;++j)std::swap(matrix[pivot*n+j],matrix[column*n+j]);std::swap(rhs[pivot],rhs[column]);}
        const T diagonal=matrix[column*n+column];
        for(std::size_t j=column;j<n;++j)matrix[column*n+j]/=diagonal;
        rhs[column]/=diagonal;
        for(std::size_t row=0;row<n;++row){if(row==column)continue;const T factor=matrix[row*n+column];if(std::abs(factor)==0.0)continue;for(std::size_t j=column;j<n;++j)matrix[row*n+j]-=factor*matrix[column*n+j];rhs[row]-=factor*rhs[column];}
    }
    return rhs;
}

Complex evaluate_polynomial(std::span<const Complex> coefficients,Complex x) {
    Complex value{};
    for(auto it=coefficients.rbegin();it!=coefficients.rend();++it)value=value*x+*it;
    return value;
}

std::vector<Complex> polynomial_roots(std::vector<Complex> coefficients) {
    double maximum=0.0;for(const auto c:coefficients)maximum=std::max(maximum,std::abs(c));
    if(maximum==0.0)return {};
    while(coefficients.size()>1U&&std::abs(coefficients.back())<maximum*1.0e-10)coefficients.pop_back();
    const std::size_t degree=coefficients.size()-1U;
    if(degree==0U)return {};
    if(degree==1U)return {-coefficients[0]/coefficients[1]};
    const Complex lead=coefficients.back();for(auto& c:coefficients)c/=lead;
    double radius=1.0;for(std::size_t i=0;i<degree;++i)radius=std::max(radius,1.0+std::abs(coefficients[i]));
    std::vector<Complex> roots(degree);
    for(std::size_t i=0;i<degree;++i){const double angle=2.0*std::numbers::pi*(static_cast<double>(i)+0.37)/static_cast<double>(degree);roots[i]=std::polar(radius,angle);}
    for(std::size_t iteration=0;iteration<300U;++iteration){
        double change=0.0;
        const auto old=roots;
        for(std::size_t i=0;i<degree;++i){
            Complex denominator{1.0,0.0};
            for(std::size_t j=0;j<degree;++j)if(j!=i)denominator*=old[i]-old[j];
            if(std::abs(denominator)<1.0e-24)denominator=Complex{1.0e-24,0.0};
            roots[i]=old[i]-evaluate_polynomial(coefficients,old[i])/denominator;
            change=std::max(change,std::abs(roots[i]-old[i]));
        }
        if(change<1.0e-12)break;
    }
    return roots;
}

} // namespace

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

SmallSignalDistortion small_signal_distortion(Circuit circuit,std::string_view source,
                                                     std::string_view output,double amplitude,
                                                     double h,const NewtonConfig& config) {
    if(!(amplitude>0.0)||!(h>0.0)||!std::isfinite(amplitude)||!std::isfinite(h))
        throw std::invalid_argument("invalid small-signal distortion controls");
    const double bias=circuit.voltage_source_dc(source);
    const auto evaluate=[&](double offset){
        circuit.set_voltage_source_dc(source,bias+offset);
        const auto op=circuit.dc_operating_point_homotopy(config);
        if(!op.converged)throw std::runtime_error("distortion operating point failed");
        return circuit.voltage(op,output);
    };
    const double fm2=evaluate(-2.0*h),fm1=evaluate(-h),f0=evaluate(0.0),fp1=evaluate(h),fp2=evaluate(2.0*h);
    const double d1=(fm2-8.0*fm1+8.0*fp1-fp2)/(12.0*h);
    const double d2=(-fp2+16.0*fp1-30.0*f0+16.0*fm1-fm2)/(12.0*h*h);
    const double d3=(fp2-2.0*fp1+2.0*fm1-fm2)/(2.0*h*h*h);
    const double fundamental=std::abs(d1*amplitude+d3*amplitude*amplitude*amplitude/8.0);
    const double second=std::abs(d2*amplitude*amplitude/4.0);
    const double third=std::abs(d3*amplitude*amplitude*amplitude/24.0);
    const double hd2=fundamental>0.0?second/fundamental:std::numeric_limits<double>::infinity();
    const double hd3=fundamental>0.0?third/fundamental:std::numeric_limits<double>::infinity();
    return {f0,d1,d2,d3,fundamental,second,third,hd2,hd3};
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

PoleZeroResult pole_zero_analysis(const Circuit& circuit,std::string_view output_node,
                                  const PoleZeroConfig& controls,const NewtonConfig& config) {
    if(!(controls.start_hz>0.0)||!(controls.stop_hz>controls.start_hz)
       ||controls.denominator_order==0U||controls.samples<3U
       ||controls.denominator_order>8U||controls.numerator_order>8U)
        throw std::invalid_argument("invalid pole-zero controls");
    const std::size_t unknowns=controls.denominator_order+controls.numerator_order+1U;
    if(controls.samples<unknowns+1U)throw std::invalid_argument("pole-zero fit needs more samples than coefficients");
    const auto sweep=circuit.ac_log_sweep(controls.start_hz,controls.stop_hz,controls.samples,config);
    const double scale=2.0*std::numbers::pi*std::sqrt(controls.start_hz*controls.stop_hz);
    std::vector<std::vector<Complex>> rows;rows.reserve(sweep.size());
    std::vector<Complex> targets;targets.reserve(sweep.size());
    for(const auto& sample:sweep){
        const Complex h=circuit.voltage(sample,output_node);
        const Complex x{0.0,2.0*std::numbers::pi*sample.frequency_hz/scale};
        std::vector<Complex> row(unknowns);Complex power=x;
        for(std::size_t j=0;j<controls.denominator_order;++j){row[j]=h*power;power*=x;}
        power={1.0,0.0};
        for(std::size_t j=0;j<=controls.numerator_order;++j){row[controls.denominator_order+j]=-power;power*=x;}
        rows.push_back(std::move(row));targets.push_back(-h);
    }
    std::vector<Complex> normal(unknowns*unknowns),rhs(unknowns);
    for(std::size_t r=0;r<rows.size();++r){
        for(std::size_t i=0;i<unknowns;++i){
            rhs[i]+=std::conj(rows[r][i])*targets[r];
            for(std::size_t j=0;j<unknowns;++j)normal[i*unknowns+j]+=std::conj(rows[r][i])*rows[r][j];
        }
    }
    double largest_diagonal=0.0;for(std::size_t i=0;i<unknowns;++i)largest_diagonal=std::max(largest_diagonal,std::abs(normal[i*unknowns+i]));
    const double ridge=std::max(1.0,largest_diagonal)*1.0e-14;
    for(std::size_t i=0;i<unknowns;++i)normal[i*unknowns+i]+=ridge;
    const auto coefficients=solve_dense_ls(std::move(normal),std::move(rhs),unknowns);
    std::vector<Complex> denominator(controls.denominator_order+1U,Complex{}),numerator(controls.numerator_order+1U,Complex{});
    denominator[0]={1.0,0.0};
    for(std::size_t j=0;j<controls.denominator_order;++j)denominator[j+1U]=coefficients[j];
    for(std::size_t j=0;j<=controls.numerator_order;++j)numerator[j]=coefficients[controls.denominator_order+j];
    double error2=0.0,signal2=0.0;
    for(std::size_t i=0;i<sweep.size();++i){
        const Complex x{0.0,2.0*std::numbers::pi*sweep[i].frequency_hz/scale};
        const Complex fit=evaluate_polynomial(numerator,x)/evaluate_polynomial(denominator,x);
        const Complex actual=circuit.voltage(sweep[i],output_node);
        error2+=std::norm(fit-actual);signal2+=std::norm(actual);
    }
    PoleZeroResult result;
    result.poles_rad_per_s=polynomial_roots(denominator);
    result.zeros_rad_per_s=polynomial_roots(numerator);
    for(auto& root:result.poles_rad_per_s)root*=scale;
    for(auto& root:result.zeros_rad_per_s)root*=scale;
    const auto order_roots=[](const Complex& a,const Complex& b){return a.real()==b.real()?a.imag()<b.imag():a.real()<b.real();};
    std::sort(result.poles_rad_per_s.begin(),result.poles_rad_per_s.end(),order_roots);
    std::sort(result.zeros_rad_per_s.begin(),result.zeros_rad_per_s.end(),order_roots);
    result.relative_rms_fit_error=std::sqrt(error2/std::max(signal2,std::numeric_limits<double>::min()));
    return result;
}

PeriodicSteadyStateResult periodic_steady_state(const Circuit& circuit,
                                                const PeriodicSteadyStateConfig& controls,
                                                const NewtonConfig& config) {
    if(!(controls.period_s>0.0)||!std::isfinite(controls.period_s)
       ||controls.samples_per_period<4U||controls.max_periods<2U
       ||!(controls.relative_tolerance>0.0)||!(controls.absolute_tolerance>0.0))
        throw std::invalid_argument("invalid periodic steady-state controls");
    if(controls.max_periods>std::numeric_limits<std::size_t>::max()/controls.samples_per_period)
        throw std::overflow_error("periodic steady-state step count overflow");
    const double dt=controls.period_s/static_cast<double>(controls.samples_per_period);
    const std::size_t total_steps=controls.max_periods*controls.samples_per_period;
    const auto transient=circuit.transient(dt,total_steps,controls.method,config);
    PeriodicSteadyStateResult result;
    std::size_t selected_period=controls.max_periods-1U;
    result.normalized_residual=std::numeric_limits<double>::infinity();
    for(std::size_t period=1U;period<controls.max_periods;++period){
        const std::size_t base=period*controls.samples_per_period;
        double residual=0.0;
        for(std::size_t phase=0;phase<=controls.samples_per_period;++phase){
            const auto& current=transient[base+phase];
            const auto& previous=transient[base-controls.samples_per_period+phase];
            if(current.node_voltage.size()!=previous.node_voltage.size())throw std::runtime_error("PSS node-vector size changed");
            for(std::size_t node=1U;node<current.node_voltage.size();++node){
                const double scale=controls.absolute_tolerance+controls.relative_tolerance*std::max(std::abs(current.node_voltage[node]),std::abs(previous.node_voltage[node]));
                residual=std::max(residual,std::abs(current.node_voltage[node]-previous.node_voltage[node])/scale);
            }
        }
        result.normalized_residual=residual;
        if(residual<=1.0){result.converged=true;result.periods=period+1U;selected_period=period;break;}
    }
    if(!result.converged)result.periods=controls.max_periods;
    const std::size_t base=selected_period*controls.samples_per_period;
    result.period.reserve(controls.samples_per_period+1U);
    const double origin=transient[base].time;
    for(std::size_t phase=0;phase<=controls.samples_per_period;++phase){
        auto point=transient[base+phase];point.time-=origin;result.period.push_back(std::move(point));
    }
    return result;
}

} // namespace cfd::circuit
