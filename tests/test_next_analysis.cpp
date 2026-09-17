#include "cfd/circuit/analysis.hpp"
#include "cfd/circuit/spice.hpp"
#include "cfd/rf/thin_wire_mom.hpp"
#include "cfd/rf/peec.hpp"
#include "cfd/rf/nport.hpp"
#include "cfd/rf/antenna.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
template<class F> void require_throws(F&& function,const char* message){
    bool threw=false;try{function();}catch(const std::exception&){threw=true;}require(threw,message);
}

void sparse_mna(){
    using namespace cfd::circuit;
    constexpr std::size_t sections=120U;
    Circuit ladder;
    std::vector<std::size_t> nodes;nodes.reserve(sections);
    for(std::size_t i=0;i<sections;++i)nodes.push_back(ladder.node("n"+std::to_string(i)));
    ladder.add_voltage_source("VIN",nodes.front(),0,1.0);
    for(std::size_t i=0;i+1U<sections;++i)
        ladder.add_resistor("R"+std::to_string(i),nodes[i],nodes[i+1U],1000.0);
    ladder.add_resistor("RL",nodes.back(),0,1000.0);

    NewtonConfig sparse;
    sparse.gmin=0.0;
    sparse.sparse_mna_threshold=1U;
    sparse.sparse_dense_fallback=false;
    sparse.sparse_relative_tolerance=1.0e-12;
    sparse.sparse_max_iterations=1000U;
    const auto op=ladder.dc_operating_point(sparse);
    require(op.converged,"forced sparse MNA DC solve");
    require(std::abs(op.node_voltage[nodes.front()]-1.0)<1.0e-10,"sparse MNA source voltage");
    require(std::abs(op.node_voltage[nodes.back()]-1.0/static_cast<double>(sections))<2.0e-10,
            "sparse MNA resistor-ladder analytical solution");

    Circuit transient;const auto in=transient.node("in"),out=transient.node("out");
    transient.add_voltage_source("V",in,0,1.0);
    transient.add_resistor("R",in,out,1000.0);
    transient.add_capacitor("C",out,0,1.0e-6);
    const auto points=transient.transient(1.0e-4,2U,TransientMethod::backward_euler,sparse);
    require(points.size()==3U&&points.back().node_voltage[out]>0.0&&points.back().node_voltage[out]<1.0,
            "forced sparse MNA transient solve");
}

void dynamic_dae_device(){
    using namespace cfd::circuit;
    constexpr double r=1000.0,c=1.0e-6;
    Circuit circuit;const auto in=circuit.node("in"),out=circuit.node("out");
    circuit.add_voltage_source("V",in,0,1.0,Complex{1.0,0.0});
    circuit.add_resistor("R",in,out,r);
    circuit.add_dynamic_device("QCAP",{out,0},[=](std::span<const double> voltage){
        const double v=voltage[0]-voltage[1],q=c*v;
        DynamicDeviceEvaluation result;
        result.terminal_current={0.0,0.0};
        result.current_jacobian={0.0,0.0,0.0,0.0};
        result.terminal_charge={q,-q};
        result.charge_jacobian={c,-c,-c,c};
        return result;
    });
    const auto op=circuit.dc_operating_point();
    require(op.converged&&std::abs(op.node_voltage[out]-1.0)<2e-8,"dynamic DAE DC conduction seam");
    const double corner=1.0/(2.0*std::numbers::pi*r*c);
    const auto ac=circuit.ac(corner,&op);
    const Complex exact_ac=Complex{1.0,0.0}/Complex{1.0,1.0};
    require(std::abs(ac.node_voltage[out]-exact_ac)<2e-8,"dynamic DAE G+jwC small-signal stamping");

    NewtonConfig forced_sparse;forced_sparse.sparse_mna_threshold=1U;forced_sparse.sparse_dense_fallback=false;
    const double dt=1.0e-4;const std::size_t steps=10U;
    const auto transient=circuit.transient(dt,steps,TransientMethod::backward_euler,forced_sparse);
    const double exact_be=1.0-std::pow(1.0/(1.0+dt/(r*c)),static_cast<double>(steps));
    require(std::abs(transient.back().node_voltage[out]-exact_be)<2e-8,"dynamic DAE backward-Euler charge history");
    require_throws([&]{(void)circuit.transient(dt,1U,TransientMethod::trapezoidal);},
                   "dynamic DAE trapezoidal must reject unsupported derivative history");
}

void electrothermal_device(){
    using namespace cfd::circuit;
    Circuit circuit;const auto electrical=circuit.node("electrical"),thermal=circuit.node("theta");
    circuit.add_voltage_source("V",electrical,0,10.0);
    ElectroThermalDeviceConfig thermal_config;thermal_config.ambient_temperature_k=300.0;
    thermal_config.thermal_resistance_k_per_w=10.0;thermal_config.thermal_capacitance_j_per_k=0.1;
    circuit.add_electrothermal_device("RHEAT",{electrical,0},thermal,thermal_config,
        [](std::span<const double> voltage,double){
            constexpr double resistance=100.0,conductance=1.0/resistance;
            const double v=voltage[0]-voltage[1],current=conductance*v;
            ElectroThermalDeviceEvaluation result;
            result.electrical_current={current,-current};
            result.electrical_jacobian={conductance,-conductance,-conductance,conductance};
            result.electrical_temperature_derivative={0.0,0.0};
            result.dissipated_power_w=v*v/resistance;
            result.power_voltage_derivative={2.0*v/resistance,-2.0*v/resistance};
            result.power_temperature_derivative_w_per_k=0.0;
            return result;
        });
    const auto op=circuit.dc_operating_point();
    require(op.converged&&std::abs(op.node_voltage[thermal]-10.0)<2e-7,
            "electrothermal DC self-heating balance");
    NewtonConfig sparse;sparse.sparse_mna_threshold=1U;sparse.sparse_dense_fallback=false;
    const double dt=0.1;const std::size_t steps=10U;
    const auto transient=circuit.transient(dt,steps,TransientMethod::backward_euler,sparse);
    const double exact=10.0*(1.0-std::pow(1.0/(1.0+dt/(thermal_config.thermal_resistance_k_per_w*thermal_config.thermal_capacitance_j_per_k)),static_cast<double>(steps)));
    require(std::abs(transient.back().node_voltage[thermal]-exact)<2e-6,
            "electrothermal thermal-RC transient");
}

void adaptive_transient(){
    using namespace cfd::circuit;
    Circuit circuit;const auto in=circuit.node("in"),out=circuit.node("out");
    circuit.add_voltage_source("VSTEP",in,0,1.0);
    circuit.add_resistor("R",in,out,1000.0);
    circuit.add_capacitor("C",out,0,1.0e-6);
    AdaptiveTransientConfig controls;controls.initial_step=1.0e-3;controls.minimum_step=1.0e-8;
    controls.maximum_step=1.0e-3;controls.relative_tolerance=2.0e-5;controls.absolute_tolerance=1.0e-8;
    const auto result=circuit.transient_adaptive(5.0e-3,controls);
    const double exact=1.0-std::exp(-5.0);
    require(result.accepted_steps>3U&&result.rejected_steps>0U,"adaptive transient must exercise LTE rejection");
    require(result.minimum_accepted_step<0.5*result.maximum_accepted_step,"adaptive transient must use variable time steps");
    require(std::abs(result.points.back().node_voltage[out]-exact)<3.0e-4,"adaptive RC step response accuracy");
    require(std::abs(result.points.back().time-5.0e-3)<1.0e-14,"adaptive transient must land on requested stop time");

    // Variable-step BDF2 uses unequal-step Gear coefficients after one BE
    // bootstrap interval. With the controller clamped to uniform accepted steps,
    // halving dt should recover the expected second-order trend.
    AdaptiveTransientConfig gear;gear.method=TransientMethod::bdf2;gear.relative_tolerance=1.0e6;gear.absolute_tolerance=1.0e6;
    gear.minimum_scale=1.0;gear.maximum_scale=1.0;gear.safety_factor=1.0;
    gear.initial_step=gear.minimum_step=gear.maximum_step=5.0e-4;
    const auto coarse=circuit.transient_adaptive(5.0e-3,gear);
    gear.initial_step=gear.minimum_step=gear.maximum_step=2.5e-4;
    const auto fine=circuit.transient_adaptive(5.0e-3,gear);
    const double coarse_error=std::abs(coarse.points.back().node_voltage[out]-exact);
    const double fine_error=std::abs(fine.points.back().node_voltage[out]-exact);
    require(fine_error<coarse_error*0.4,"adaptive variable-step BDF2 must show second-order RC refinement");

    AdaptiveTransientConfig variable=controls;variable.method=TransientMethod::bdf2;
    const auto nonuniform=circuit.transient_adaptive(5.0e-3,variable);
    const double nonuniform_error=std::abs(nonuniform.points.back().node_voltage[out]-exact);
    require(nonuniform.minimum_accepted_step<0.5*nonuniform.maximum_accepted_step,
            "adaptive BDF2 must exercise unequal accepted time steps");
    require(nonuniform_error<2.0e-4,"unequal-step adaptive BDF2 RC accuracy");
}

void pss_and_pole_zero(){
    using namespace cfd::circuit;
    constexpr double r=1000.0,c=1.0e-6,f=100.0;
    Circuit periodic;const auto in=periodic.node("in"),out=periodic.node("out");
    periodic.add_voltage_source("VS",in,0,0.0,Complex{1.0,0.0},SineWaveform{0.0,1.0,f,0.0,0.0,0.0});
    periodic.add_resistor("R",in,out,r);periodic.add_capacitor("C",out,0,c);
    PeriodicSteadyStateConfig pss;pss.period_s=1.0/f;pss.samples_per_period=128U;pss.max_periods=20U;
    pss.relative_tolerance=2.0e-4;pss.absolute_tolerance=1.0e-7;
    const auto steady=periodic_steady_state(periodic,pss);
    require(steady.converged&&steady.periods<pss.max_periods,"RC PSS must converge by cycle-to-cycle shooting residual");
    const auto measurement=fourier_measurement(steady.period,out,f);
    const double expected=1.0/std::sqrt(1.0+std::pow(2.0*std::numbers::pi*f*r*c,2));
    require(std::abs(measurement.magnitude-expected)<2.0e-3,"RC PSS amplitude must agree with small-signal solution");

    Circuit lowpass;const auto li=lowpass.node("in"),lo=lowpass.node("out");
    lowpass.add_voltage_source("VAC",li,0,0.0,Complex{1.0,0.0});lowpass.add_resistor("R",li,lo,r);lowpass.add_capacitor("C",lo,0,c);
    PoleZeroConfig pz;pz.start_hz=1.0;pz.stop_hz=1.0e5;pz.samples=40U;pz.denominator_order=1U;pz.numerator_order=0U;
    const auto lp=pole_zero_analysis(lowpass,"out",pz);
    require(lp.poles_rad_per_s.size()==1U&&lp.zeros_rad_per_s.empty(),"RC low-pass pole count");
    require(std::abs((lp.poles_rad_per_s[0].real()+1000.0)/1000.0)<2.0e-4&&std::abs(lp.poles_rad_per_s[0].imag())<1.0e-5,"RC low-pass pole location");
    require(lp.relative_rms_fit_error<1.0e-6,"RC low-pass rational fit accuracy");

    Circuit highpass;const auto hi=highpass.node("in"),ho=highpass.node("out");
    highpass.add_voltage_source("VAC",hi,0,0.0,Complex{1.0,0.0});highpass.add_capacitor("C",hi,ho,c);highpass.add_resistor("R",ho,0,r);
    pz.numerator_order=1U;const auto hp=pole_zero_analysis(highpass,"out",pz);
    require(hp.poles_rad_per_s.size()==1U&&hp.zeros_rad_per_s.size()==1U,"RC high-pass pole/zero count");
    require(std::abs((hp.poles_rad_per_s[0].real()+1000.0)/1000.0)<2.0e-4,"RC high-pass pole location");
    require(std::abs(hp.zeros_rad_per_s[0])<1.0e-3,"RC high-pass zero at origin");
}


void transformer_and_lossy_line(){
    using namespace cfd::circuit;
    Circuit transformer;const auto primary=transformer.node("p"),secondary=transformer.node("s");
    transformer.add_current_source("IIN",0,primary,1.0);
    transformer.add_resistor("RL",secondary,0,100.0);
    transformer.add_ideal_transformer("T1",primary,0,secondary,0,2.0);
    const auto op=transformer.dc_operating_point();
    require(op.converged,"ideal transformer DC operating point");
    require(std::abs(transformer.voltage(op,"p")-400.0)<5.0e-7,"ideal transformer reflected load voltage");
    require(std::abs(transformer.voltage(op,"s")-200.0)<5.0e-7,"ideal transformer secondary voltage");
    require(std::abs(transformer.voltage(op,"p")/transformer.voltage(op,"s")-2.0)<1.0e-12,"ideal transformer turns ratio");

    constexpr double frequency=1.0e9;const std::array<double,1> samples{{frequency}};
    Circuit line;const auto a=line.node("a"),b=line.node("b");
    cfd::rf::TemTransmissionLine lossy{50.0,299792458.0,0.1,1.0};
    line.add_tem_transmission_line("TL",{a,0},{b,0},lossy,samples,50.0);
    const auto sparams=two_port_s_parameters(line,{"a","0"},{"b","0"},frequency,50.0);
    require(std::abs(sparams.a11)<1.0e-10,"matched lossy transmission line reflection");
    require(std::abs(std::abs(sparams.a21)-std::exp(-0.1))<1.0e-9,"lossy transmission-line attenuation");
}



void compact_device_noise(){
    using namespace cfd::circuit;
    const auto has_positive=[](const NoiseResult& result,const std::string& prefix){
        for(const auto& item:result.contributions)if(item.source.rfind(prefix,0)==0U&&item.output_noise_density_v2_per_hz>0.0)return true;
        return false;
    };
    Circuit mos;auto vdd=mos.node("vdd"),gate=mos.node("g"),drain=mos.node("d");
    mos.add_voltage_source("VDD",vdd,0,5.0);mos.add_voltage_source("VG",gate,0,2.0);mos.add_resistor("RD",vdd,drain,2000.0);
    mos.add_mos_model("NM",{false,0.7,1.0e-3,0.02});mos.add_mosfet("M1",drain,gate,0,0,"NM");
    const auto mn=mos.output_noise(1.0e3,"d");
    require(has_positive(mn,"M1")&&mn.output_noise_density_v2_per_hz>0.0,"MOS channel thermal-noise contribution");

    Circuit jf;vdd=jf.node("vdd");drain=jf.node("d");
    jf.add_voltage_source("VDD",vdd,0,5.0);jf.add_resistor("RD",vdd,drain,1000.0);
    jf.add_jfet_model("NJ",{false,-2.0,1.0e-3,0.01});jf.add_jfet("J1",drain,0,0,"NJ");
    const auto jn=jf.output_noise(1.0e3,"d");
    require(has_positive(jn,"J1"),"JFET channel thermal-noise contribution");

    Circuit bjt;vdd=bjt.node("vcc");auto vbb=bjt.node("vbb"),base=bjt.node("b"),collector=bjt.node("c");
    bjt.add_voltage_source("VCC",vdd,0,5.0);bjt.add_voltage_source("VBB",vbb,0,5.0);
    bjt.add_resistor("RC",vdd,collector,1000.0);bjt.add_resistor("RB",vbb,base,100000.0);
    bjt.add_bjt_model("QN",{false,1e-15,100.0,1.0,300.0});bjt.add_bjt("Q1",collector,base,0,"QN");
    const auto qn=bjt.output_noise(1.0e3,"c");
    require(has_positive(qn,"Q1:collector")&&has_positive(qn,"Q1:base"),"BJT collector/base shot-noise contributions");
    const auto spectrum=output_noise_spectrum(bjt,"c",10.0,1.0e4,9U,300.0);
    require(spectrum.points.size()==9U&&spectrum.integrated_rms_voltage>0.0,
            "multi-source noise sweep and integration");
}

void distortion_and_rf_adaptive_tools(){
    using namespace cfd::circuit;
    const char* nonlinear_net=R"SPICE(
VIN in 0 0
BOUT out 0 V={V(in)+0.1*V(in)*V(in)+0.01*V(in)*V(in)*V(in)}
RLOAD out 0 1k
)SPICE";
    auto nonlinear=Circuit::parse_spice(nonlinear_net);
    const auto distortion=small_signal_distortion(nonlinear,"VIN","out",0.1,1.0e-3);
    require(std::abs(distortion.first_derivative-1.0)<1.0e-6,"distortion first derivative");
    require(std::abs(distortion.second_derivative-0.2)<2.0e-6,"distortion second derivative");
    require(std::abs(distortion.third_derivative-0.06)<2.0e-5,"distortion third derivative");
    require(std::abs(distortion.hd2-0.004999625)<2.0e-6&&std::abs(distortion.hd3-2.4998e-5)<2.0e-7,
            "small-signal harmonic distortion prediction");

    using namespace cfd::rf;
    const NoiseParameters noise{1.2,{0.2,0.1},10.0,50.0};
    const auto circle=noise_circle(noise,1.5);
    require(circle.valid&&circle.radius>0.0,"RF noise circle geometry");
    const Complex on_circle=circle.center+Complex{circle.radius,0.0};
    require(std::abs(noise_factor(noise,on_circle)-1.5)<1.0e-12,"RF noise circle constant-noise factor");
    require(std::abs(noise_factor(noise,noise.optimum_source_reflection)-noise.minimum_noise_factor)<1.0e-14,
            "RF optimum source reflection reaches Fmin");

    const std::array<double,1> reference{{50.0}};
    AdaptiveNetworkSweepConfig sweep_controls;sweep_controls.initial_points=3U;sweep_controls.maximum_points=129U;
    sweep_controls.relative_tolerance=5.0e-3;sweep_controls.absolute_tolerance=1.0e-6;
    const auto sweep=adaptive_network_sweep(1.0e3,1.0e7,reference,
        [](double frequency){
            ComplexMatrix matrix(1U);const Complex pole{1.0,frequency/1.0e5};matrix(0,0)=Complex{0.8,0.0}/pole;return matrix;
        },sweep_controls);
    require(sweep.converged&&sweep.points.size()>sweep_controls.initial_points,
            "adaptive RF network sweep refines curved response");
    require(sweep.points.front().frequency_hz==1.0e3&&sweep.points.back().frequency_hz==1.0e7,
            "adaptive RF network sweep preserves requested endpoints");
}

void rf_network_and_peec_extensions(){
    using namespace cfd::rf;
    ComplexMatrix z(2);z(0,0)={70.0,10.0};z(1,1)={45.0,-5.0};z(0,1)=z(1,0)={8.0,2.0};
    const std::array<Complex,2> complex_ref{{{50.0,10.0},{75.0,-15.0}}};
    const auto s=z_to_s_power_wave(z,complex_ref);
    const auto roundtrip=s_to_z_power_wave(s,complex_ref);
    double roundtrip_error=0.0;
    for(std::size_t r=0;r<2;++r)for(std::size_t c=0;c<2;++c)
        roundtrip_error=std::max(roundtrip_error,std::abs(roundtrip(r,c)-z(r,c)));
    require(roundtrip_error<1.0e-11,"complex power-wave Z/S roundtrip");
    const std::array<Complex,2> new_ref{{{60.0,-5.0},{40.0,7.0}}};
    const auto renorm=renormalize_s_power_wave(s,complex_ref,new_ref);
    const auto recovered=s_to_z_power_wave(renorm,new_ref);
    double renorm_error=0.0;
    for(std::size_t r=0;r<2;++r)for(std::size_t c=0;c<2;++c)
        renorm_error=std::max(renorm_error,std::abs(recovered(r,c)-z(r,c)));
    require(renorm_error<1.0e-11,"complex-reference power-wave renormalization preserves Z");

    ComplexMatrix active(2);active(0,0)={0.2,0.0};active(1,0)={1.4,0.0};active(0,1)={0.1,0.0};active(1,1)={0.1,0.0};
    require(!network_quality(active).passive,"active network identified before projection");
    const auto passive=enforce_passivity(active,0.999);
    require(network_quality(passive,1.0e-12).maximum_singular_value<=0.9990000001,
            "singular-value passivity projection");

    std::vector<NPortPoint> causal_points,advance_points;
    constexpr std::size_t causal_count=17U;constexpr double df=1.0e6;
    constexpr double sample_dt=1.0/(32.0*df);constexpr double delay=3.0*sample_dt;
    for(std::size_t k=0;k<causal_count;++k){
        const double frequency=df*static_cast<double>(k);
        NPortPoint causal;causal.frequency_hz=frequency;causal.s=ComplexMatrix(1U);causal.reference_impedance={50.0};
        causal.s(0,0)=std::exp(Complex{0.0,-2.0*std::numbers::pi*frequency*delay});
        NPortPoint advance=causal;advance.s(0,0)=std::exp(Complex{0.0,2.0*std::numbers::pi*frequency*delay});
        causal_points.push_back(std::move(causal));advance_points.push_back(std::move(advance));
    }
    const auto causal_gate=check_broadband_causality(causal_points,1.0e-10);
    const auto advance_gate=check_broadband_causality(advance_points,1.0e-3);
    require(causal_gate.causal&&causal_gate.negative_time_energy_ratio<1.0e-12,
            "broadband causality gate accepts a discrete pure delay");
    require(!advance_gate.causal&&advance_gate.negative_time_energy_ratio>0.99,
            "broadband causality gate rejects a phase advance");

    PeecFilamentSystem peec;
    peec.add_segment({{0,0,0},{0,0,0.1},5.0e-4,5.8e7});
    peec.add_segment({{0.02,0,0},{0.02,0,0.1},5.0e-4,5.8e7});
    const auto cap=peec.extract_capacitance();
    require(cap.size==2U&&cap.coefficient_of_potential_v_per_c[0]>0.0,
            "PEEC coefficient-of-potential extraction");
    require(std::abs(cap.coefficient_of_potential_v_per_c[1]-cap.coefficient_of_potential_v_per_c[2])<1.0e-6,
            "PEEC potential matrix reciprocity");
    require(cap.capacitance_f[0]>0.0&&cap.capacitance_f[3]>0.0&&cap.capacitance_f[1]<0.0,
            "PEEC Maxwell capacitance sign structure");
    const auto cap_er4=peec.extract_capacitance(4.0);
    require(std::abs(cap_er4.capacitance_f[0]/cap.capacitance_f[0]-4.0)<2.0e-12,
            "PEEC capacitance scales with permittivity");

    const auto match=antenna_match_metrics({50.0,0.0},50.0);
    require(std::abs(match.reflection_coefficient)<1.0e-15&&std::isinf(match.return_loss_db)&&std::abs(match.vswr-1.0)<1.0e-15,
            "antenna matched-feed metrics");
}

void antenna_optimization(){
    using namespace cfd::rf;
    const AntennaOptimizationVariable variable{0.30,0.20,0.80,0.12};
    AntennaOptimizationConfig controls;controls.parameter_tolerance=1.0e-7;controls.max_iterations=96U;
    const auto result=optimize_antenna(std::span<const AntennaOptimizationVariable>(&variable,1U),
        [](std::span<const double> parameters){
            const double detuning=200.0*(parameters[0]-0.5);
            const auto match=antenna_match_metrics({50.0,detuning},50.0);
            return std::norm(match.reflection_coefficient);
        },controls);
    require(result.converged&&result.evaluations>1U&&std::abs(result.parameters[0]-0.5)<2.0e-6,
            "bounded antenna matching optimization driver");
    require(result.objective<1.0e-10,"antenna optimizer reaches matched objective");
}

void multi_wire_mom(){
    using namespace cfd::rf;
    constexpr double frequency=3.0e8;
    const double wavelength=299792458.0/frequency;
    ParallelWireMomConfig config;config.frequency_hz=frequency;config.quadrature_order=8U;
    config.wires={{-0.05*wavelength,0.0,0.0,0.47*wavelength,0.001*wavelength,21U},
                  { 0.05*wavelength,0.0,0.0,0.47*wavelength,0.001*wavelength,21U}};
    config.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    const auto first=solve_parallel_thin_wires(config);
    double driven=0.0,induced=0.0;
    for(std::size_t i=0;i<21U;++i)driven=std::max(driven,std::abs(first.current_a[i]));
    for(std::size_t i=21U;i<42U;++i)induced=std::max(induced,std::abs(first.current_a[i]));
    require(driven>0.0&&induced>driven*1.0e-6,"parallel-wire MoM must produce mutual induced current");
    config.feeds={{1U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    const auto second=solve_parallel_thin_wires(config);
    double symmetry=0.0,reference=0.0;
    for(std::size_t i=0;i<21U;++i){
        symmetry=std::max(symmetry,std::abs(first.current_a[i]-second.current_a[21U+i]));
        symmetry=std::max(symmetry,std::abs(first.current_a[21U+i]-second.current_a[i]));
        reference=std::max(reference,std::abs(first.current_a[i]));
    }
    require(symmetry/std::max(reference,1.0e-30)<1.0e-9,"parallel-wire MoM reciprocal geometry symmetry");
    config.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    config.loads={{0U,10U,{1.0e6,0.0}}};
    const auto loaded=solve_parallel_thin_wires(config);
    require(loaded.feed_current_a[0]<first.feed_current_a[0],"series wire load must reduce driven feed current");

    // Rotate the original two-wire geometry onto an arbitrary 3-D axis. The
    // free-space dyadic kernel must be rotationally invariant.
    const double inv_d=1.0/std::sqrt(14.0),inv_n=1.0/std::sqrt(5.0);
    const WirePoint3 d{inv_d,2.0*inv_d,3.0*inv_d};
    const WirePoint3 n{2.0*inv_n,-inv_n,0.0};
    OrientedWireMomConfig rotated;rotated.frequency_hz=frequency;rotated.quadrature_order=8U;
    for(double sign:{-1.0,1.0}){
        const WirePoint3 center{sign*0.05*wavelength*n.x,sign*0.05*wavelength*n.y,0.0};
        rotated.wires.push_back({{center.x-0.5*0.47*wavelength*d.x,center.y-0.5*0.47*wavelength*d.y,center.z-0.5*0.47*wavelength*d.z},
                                 {center.x+0.5*0.47*wavelength*d.x,center.y+0.5*0.47*wavelength*d.y,center.z+0.5*0.47*wavelength*d.z},
                                 0.001*wavelength,21U});
    }
    rotated.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    const auto arbitrary=solve_oriented_thin_wires(rotated);
    require(arbitrary.current_a.size()==first.current_a.size(),"oriented-wire segment count parity");
    double rotation_error=0.0,rotation_reference=0.0;
    for(std::size_t i=0;i<first.current_a.size();++i){
        rotation_error=std::max(rotation_error,std::abs(arbitrary.current_a[i]-first.current_a[i]));
        rotation_reference=std::max(rotation_reference,std::abs(first.current_a[i]));
    }
    require(rotation_error/std::max(rotation_reference,1.0e-30)<2.0e-10,"oriented-wire MoM rotational invariance");

    rotated.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}},
                   {1U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    const auto zports=oriented_wire_port_impedance_matrix(rotated);
    require(zports.size()==2U&&zports(0,0).real()>0.0&&zports(1,1).real()>0.0,
            "two-antenna MoM port-impedance extraction");
    require(std::abs(zports(0,1)-zports(1,0))/std::max(std::abs(zports(0,1)),1.0e-30)<2.0e-8,
            "two-antenna MoM mutual impedance reciprocity");
    const std::array<double,2> z0{{50.0,50.0}};
    const auto sport=oriented_wire_port_s_parameters(rotated,z0);
    require(std::isfinite(std::abs(sport(1,0)))&&std::abs(sport(1,0))>1.0e-8,
            "two-antenna MoM coupling S21");

    OrientedWireMomConfig line_loaded=rotated;line_loaded.feeds.resize(1U);
    line_loaded.transmission_line_loads.push_back({0U,10U,1U,10U,{75.0,0.0},
                                                   {0.2,2.0*std::numbers::pi/wavelength},0.13*wavelength,1,-1});
    const auto line_loaded_result=solve_oriented_thin_wires(line_loaded);
    require(line_loaded_result.transmission_line_loss_w>0.0&&line_loaded_result.accepted_power_w>0.0,
            "thin-wire lossy transmission-line section load");
    require(std::abs(line_loaded_result.feed_impedance_ohm.front()-arbitrary.feed_impedance_ohm.front())>1.0e-4,
            "wire transmission-line section must alter feed impedance");

    OrientedWireMomConfig ground;ground.frequency_hz=frequency;ground.quadrature_order=8U;
    ground.ground_model=ThinWireGroundModel::pec_plane;ground.ground_plane_z_m=0.0;
    ground.wires={{{-0.235*wavelength,0.0,0.15*wavelength},{0.235*wavelength,0.0,0.15*wavelength},
                   0.001*wavelength,21U}};
    ground.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    const auto grounded=solve_oriented_thin_wires(ground);
    OrientedWireMomConfig shifted_ground=ground;shifted_ground.ground_plane_z_m=0.7*wavelength;
    shifted_ground.wires[0].start_m.z+=0.7*wavelength;shifted_ground.wires[0].end_m.z+=0.7*wavelength;
    const auto shifted_grounded=solve_oriented_thin_wires(shifted_ground);
    require(std::abs(grounded.feed_impedance_ohm.front()-shifted_grounded.feed_impedance_ohm.front())
            /std::max(std::abs(grounded.feed_impedance_ohm.front()),1.0e-30)<2.0e-10,
            "PEC-ground image kernel translation invariance");
    auto free_horizontal=ground;free_horizontal.ground_model=ThinWireGroundModel::free_space;
    const auto free_horizontal_result=solve_oriented_thin_wires(free_horizontal);
    require(std::abs(grounded.feed_impedance_ohm.front()-free_horizontal_result.feed_impedance_ohm.front())>1.0e-3,
            "PEC-ground image interaction changes antenna feed impedance");

    OrientedWireMomConfig lossy_wire;lossy_wire.frequency_hz=frequency;lossy_wire.quadrature_order=8U;
    OrientedThinWire copper{{0.0,0.0,-0.235*wavelength},{0.0,0.0,0.235*wavelength},0.001*wavelength,21U};
    copper.conductivity_s_per_m=5.8e7;copper.relative_permeability=1.0;copper.dielectric_loss_ohm_per_m=0.25;
    lossy_wire.wires={copper};lossy_wire.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    const auto lossy_result=solve_oriented_thin_wires(lossy_wire);
    require(lossy_result.accepted_power_w>0.0&&lossy_result.conductor_loss_w>0.0&&lossy_result.dielectric_loss_w>0.0,
            "thin-wire distributed conductor/dielectric loss accounting");
    require(lossy_result.radiated_power_w>0.0&&lossy_result.radiation_efficiency>0.0&&lossy_result.radiation_efficiency<1.0,
            "thin-wire radiation efficiency includes distributed losses");

    OrientedWireMomConfig connected;connected.frequency_hz=frequency;connected.quadrature_order=8U;
    connected.junction_tolerance_m=1.0e-9*wavelength;
    connected.wires={{{-0.48*wavelength,0.0,0.0},{0.0,0.0,0.0},0.001*wavelength,11U},
                     {{0.0,0.0,0.0},{0.48*wavelength,0.0,0.0},0.001*wavelength,11U}};
    connected.feeds={{0U,5U,{1.0,0.0}}};
    const auto joined=solve_oriented_thin_wires(connected);
    require(joined.current_a.size()==22U&&joined.feed_impedance_ohm.front().real()>0.0,
            "endpoint-connected wire junction solve");
    const Complex junction_residual=joined.current_a[10U]-joined.current_a[11U];
    const double junction_scale=std::max({std::abs(joined.current_a[10U]),std::abs(joined.current_a[11U]),1.0e-30});
    require(std::abs(junction_residual)/junction_scale<1.0e-10,
            "connected-wire junction enforces current continuity/KCL");

    OrientedWireMomConfig crossing=rotated;crossing.wires={
        {{-0.25*wavelength,0.0,0.0},{0.25*wavelength,0.0,0.0},0.001*wavelength,21U},
        {{0.0,-0.25*wavelength,0.0},{0.0,0.25*wavelength,0.0},0.001*wavelength,21U}};
    crossing.feeds={{0U,std::numeric_limits<std::size_t>::max(),{1.0,0.0}}};
    require_throws([&]{(void)solve_oriented_thin_wires(crossing);},
                   "crossing wires must be rejected until junction basis functions exist");
}
}

int main(){
    try{sparse_mna();dynamic_dae_device();electrothermal_device();adaptive_transient();pss_and_pole_zero();transformer_and_lossy_line();compact_device_noise();distortion_and_rf_adaptive_tools();rf_network_and_peec_extensions();antenna_optimization();multi_wire_mom();std::cout<<"next analysis tests passed\n";return 0;}
    catch(const std::exception& error){std::cerr<<"next analysis test failed: "<<error.what()<<'\n';return 1;}
}
