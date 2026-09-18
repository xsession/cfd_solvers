#include "cfd/battery/pack.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace cfd::battery;
namespace {
void check(bool ok, const char* msg) { if (!ok) { std::cerr << msg << '\n'; std::exit(1); } }
bool near(double a, double b) { return std::abs(a-b) < 1e-10 * (1.0 + std::abs(a)); }
template<class F> void rejects(F f) { bool threw=false; try { f(); } catch (const std::exception&) { threw=true; } check(threw,"expected rejection"); }
}
int main() {
    auto c = default_graphite_nmc_config();
    std::vector<LithiumIonCellConfig> configs(3,c);
    SeriesBatteryPack pack(configs);
    SingleParticleElectrolyteModel ref(c);
    for (int n=0;n<10;++n) {
        const auto r=ref.step(0.2,0.5);
        const auto p=pack.step(0.2,0.5);
        check(near(p.terminal_voltage_v,3*r.terminal_voltage_v),"series voltage sum and independent cell parity");
        check(near(p.net_terminal_charge_ah,0.2*0.5*(n+1)/3600),"terminal charge integral");
        check(p.balancing_heat_w==0,"balancing disabled");
    }
    // State must retain electrolyte polarization and concentration diagnostics.
    const auto query=ref.state(0.2);
    check(query.electrolyte_overpotential_v!=0,"SPMe state electrolyte polarization");
    check(query.electrolyte_min_concentration_mol_per_m3<query.electrolyte_max_concentration_mol_per_m3,"SPMe state concentration range");
    pack.reset(); ref.reset();
    check(pack.time_s()==0 && pack.cells()[0].time_s()==0,"reset clocks");
    const auto charged=pack.step(-0.1,1.0);
    check(charged.terminal_energy_wh<0 && charged.net_terminal_charge_ah<0,"charging sign");

    configs[0].positive.particle.initial_stoichiometry=0.30;
    PackControl control; control.balancing_current_a=0.05;
    SeriesBatteryPack balanced(configs,control);
    std::vector<SingleParticleElectrolyteModel> independent;
    for (auto cfg: configs) independent.emplace_back(cfg);
    auto p=balanced.step(-0.1,1.0);
    check(p.balancing_currents_a[0]==0.05 && p.balancing_currents_a[1]==0,"high voltage cell selectively bleeds");
    double cell_power=0;
    for (std::size_t i=0;i<configs.size();++i) {
        const double current=-0.1+p.balancing_currents_a[i];
        auto r=independent[i].step(current,1.0);
        check(near(p.cells[i].terminal_voltage_v,r.terminal_voltage_v),"heterogeneous cell reference parity");
        check(near(p.cells[i].current_a,current),"shunt Kirchhoff current balance");
        cell_power+=r.terminal_voltage_v*current;
    }
    check(near(cell_power,p.terminal_power_w+p.balancing_heat_w),"pack/shunt electrical power conservation");
    check(near(p.balancing_energy_wh,p.balancing_heat_w/3600),"shunt energy accounting");
    SeriesBatteryPack idle(configs,control);
    const auto rest=idle.step(0,1);
    check(rest.terminal_energy_wh==0 && rest.balancing_energy_wh>0,"rest balancing dissipates only cell energy");

    // Bleeding at rest must remove exactly I*dt/F moles from the negative electrode.
    const auto& neg = idle.cells()[0].negative_particle();
    const double removed = (configs[0].negative.particle.initial_stoichiometry
        * configs[0].negative.particle.maximum_concentration_mol_per_m3
        - neg.average_concentration()) * configs[0].area_m2
        * configs[0].negative.thickness_m * configs[0].negative.active_volume_fraction;
    check(std::abs(removed-0.05/faraday_constant)<1e-13,"bleed lithium inventory conservation");
    control.maximum_cell_voltage_v=3.0;
    SeriesBatteryPack overvoltage(configs,control);
    rejects([&]{ (void)overvoltage.step(-0.1,1); });
    check(overvoltage.time_s()==0,"overvoltage rollback");

    // A failure in the second cell occurs after the first trial cell advanced.
    configs=std::vector<LithiumIonCellConfig>(2,c);
    configs[1].contact_resistance_ohm=100;
    SeriesBatteryPack cutoff(configs);
    const auto before=cutoff.cells()[0].negative_particle().concentrations();
    rejects([&]{ (void)cutoff.step(0.2,1); });
    check(cutoff.time_s()==0 && cutoff.cells()[0].time_s()==0 && cutoff.cells()[1].time_s()==0,"cutoff atomic clocks");
    check(cutoff.cells()[0].negative_particle().concentrations()==before,"cutoff atomic particle state");
    const auto recovery=cutoff.step(0,1);
    check(recovery.net_terminal_charge_ah==0 && recovery.terminal_energy_wh==0,"rejected step did not alter integrals");

    configs=std::vector<LithiumIonCellConfig>(1,c);
    control=PackControl{}; control.maximum_temperature_k=c.initial_temperature_k+1e-8;
    SeriesBatteryPack hot(configs,control);
    rejects([&]{ (void)hot.step(0.2,1); });
    check(hot.time_s()==0,"trial-end thermal cutoff rollback");
    rejects([&]{ (void)hot.step(0,std::numeric_limits<double>::infinity()); });
    rejects([&]{ (void)hot.step(std::numeric_limits<double>::quiet_NaN(),1); });
    rejects([&]{ SeriesBatteryPack empty(std::span<const LithiumIonCellConfig>{}); });
    control.balancing_current_a=-1;
    rejects([&]{ SeriesBatteryPack invalid(configs,control); });
    c.electrolyte.cells=0;
    rejects([&]{ SingleParticleElectrolyteModel invalid(c); });
    std::cout << "v0.19.1 pack regression passed\n";
}
