#include "cfd/battery/thermal3d.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace cfd::battery;
namespace {
void check(bool ok,const char* msg) { if(!ok) { std::cerr<<msg<<'\n'; std::exit(1); } }
bool near(double a,double b,double eps=1e-9) { return std::abs(a-b)<eps*(1+std::abs(a)); }
template<class F> void rejects(F f) { bool threw=false; try {f();} catch(const std::exception&) {threw=true;} check(threw,"expected rejection"); }
}
int main() {
    // Analytical discrete 3D Neumann eigenmode exercises every spatial axis.
    ThermalGridConfig g; g.cells={4,3,2}; g.spacing_m={0.1,0.2,0.3};
    ThermalVoxel material{1000,{1,2,3}};
    ThermalGrid3D grid(g,{material});
    std::vector<double> t(24),power(24,0);
    double eigenvalue=0;
    for(std::size_t a=0;a<3;++a) {
        const double s=std::sin(std::numbers::pi/(2*static_cast<double>(g.cells[a])));
        eigenvalue+=4*material.conductivity_w_per_m_k[a]*s*s/(1000*g.spacing_m[a]*g.spacing_m[a]);
    }
    for(std::size_t z=0;z<2;++z) for(std::size_t y=0;y<3;++y) for(std::size_t x=0;x<4;++x)
        t[x+4*(y+3*z)]=300+5*std::cos(std::numbers::pi*(static_cast<double>(x)+0.5)/4)
            *std::cos(std::numbers::pi*(static_cast<double>(y)+0.5)/3)
            *std::cos(std::numbers::pi*(static_cast<double>(z)+0.5)/2);
    grid.reset(t); const auto mode=grid.step(2,power);
    for(std::size_t i=0;i<t.size();++i)
        check(near(grid.temperature()[i],300+(t[i]-300)/(1+2*eigenvalue),1e-11),"3D anisotropic diffusion eigenmode");
    check(std::abs(mode.stored_energy_change_j)<1e-9,"insulated mode energy conservation");

    // Heterogeneous capacities/conductivities: internal face flux cancels.
    std::vector<ThermalVoxel> materials(24,material);
    for(std::size_t i=0;i<24;++i) {
        materials[i].volumetric_heat_capacity_j_per_m3_k+=100*static_cast<double>(i);
        materials[i].conductivity_w_per_m_k={1+static_cast<double>(i),0.2,4};
        power[i]=0.1*static_cast<double>(i);
    }
    ThermalGrid3D heterogeneous(g,materials); heterogeneous.reset(t);
    const auto energy=heterogeneous.step(3,power);
    check(near(energy.stored_energy_change_j,3*energy.supplied_power_w),"heterogeneous source energy balance");
    check(std::abs(energy.energy_balance_error_j)<1e-8,"reported energy residual");
    const auto saved=heterogeneous.temperature();
    power[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{(void)heterogeneous.step(1,power);});
    check(heterogeneous.temperature()==saved,"invalid source atomic rollback");

    // Exact backward-Euler single-volume convection including half-cell resistance.
    ThermalGridConfig one; one.cells={1,1,1}; one.spacing_m={1,1,1}; one.initial_temperature_k=310;
    one.boundaries[0]={10,300};
    ThermalGrid3D cooled(one,{ThermalVoxel{10,{2,2,2}}});
    const std::vector<double> zero{0}; const double conductance=1/(0.5/2+1.0/10);
    const auto cool=cooled.step(2,zero);
    check(near(cooled.temperature()[0],(10*310+2*conductance*300)/(10+2*conductance)),"analytical convection decay");
    check(near(cool.stored_energy_change_j,-2*cool.outward_cooling_power_w),"cooling energy balance");

    // Thermal capacitance equals configured cell mass*cp; zero lumped cooling.
    auto cell=default_graphite_nmc_config(); cell.thermal.heat_transfer_coefficient_w_per_m2_k=0;
    std::vector<LithiumIonCellConfig> cells(2,cell);
    ThermalGridConfig pg; pg.cells={2,1,1}; pg.spacing_m={0.01,0.01,0.01};
    const double capacity=cell.thermal.mass_kg*cell.thermal.heat_capacity_j_per_kg_k;
    ThermalVoxel pm{capacity/1e-6,{1,1,1}};
    ElectrothermalBatteryPack coupled(cells,{},ThermalGrid3D(pg,{pm}),{{0},{1}});
    SingleParticleElectrolyteModel reference(cell);
    const auto ref=reference.step(0.2,1);
    const auto result=coupled.step(0.2,1);
    check(near(result.electrical.cells[0].temperature_k,ref.temperature_k,1e-12),"no double-counted cell heating");
    check(near(result.thermal.supplied_power_w,2*ref.heat_generation_w),"conservative cell heat deposition");
    check(near(coupled.electrical().cells()[0].temperature_k(),coupled.thermal().temperature()[0]),"thermal feedback owns cell temperature");
    const auto next=coupled.step(0.2,1);
    check(next.electrical.time_s==2,"coupled continuation clock");
    check(std::abs(next.thermal.energy_balance_error_j)<1e-10,"coupled energy conservation");

    PackControl tight; tight.maximum_temperature_k=pg.initial_temperature_k+1e-8;
    ElectrothermalBatteryPack cutoff(cells,tight,ThermalGrid3D(pg,{pm}),{{0},{1}});
    const auto initial=cutoff.thermal().temperature();
    rejects([&]{(void)cutoff.step(0.2,1);});
    check(cutoff.electrical().time_s()==0 && cutoff.thermal().temperature()==initial,"thermal cutoff rolls back both systems");
    check(cutoff.electrical().cells()[0].negative_particle().average_stoichiometry()==cell.negative.particle.initial_stoichiometry
          || near(cutoff.electrical().cells()[0].negative_particle().average_stoichiometry(),cell.negative.particle.initial_stoichiometry),"cutoff preserves lithium");
    rejects([&]{ElectrothermalBatteryPack bad(cells,{},ThermalGrid3D(pg,{pm}),{{0},{0}});});
    rejects([&]{ElectrothermalBatteryPack bad(cells,{},ThermalGrid3D(pg),{{0},{1}});});
    auto badg=pg; badg.cells[0]=0;
    rejects([&]{ThermalGrid3D bad(badg);});
    badg=g; badg.maximum_iterations=1;
    ThermalGrid3D failed(badg,materials); failed.reset(t);
    power[0]=3;
    rejects([&]{(void)failed.step(1,power);});
    check(failed.temperature()==t,"unconverged solve rollback");

    SeriesBatteryPack external(cells);
    const std::vector<double> external_t{300,301}; external.set_external_temperatures(external_t);
    const auto isothermal=external.step(0.2,1);
    check(isothermal.cells[0].temperature_k==300,"external temperature suppresses lumped solve");
    external.reset(); const auto lumped=external.step(0.2,1);
    check(lumped.cells[0].temperature_k>cell.initial_temperature_k,"reset restores lumped ownership");
    std::cout<<"v0.19.2 thermal regression passed\n";
}
