#include "cfd/battery/thermal3d.hpp"
#include <iostream>
int main() {
    using namespace cfd::battery;
    std::vector<LithiumIonCellConfig> cells(3,default_graphite_nmc_config());
    cells[0].contact_resistance_ohm=0.05;
    ThermalGridConfig grid; grid.cells={6,2,2}; grid.spacing_m={0.01,0.01,0.01};
    for(auto& b:grid.boundaries) b={8,298.15};
    const double capacity=cells[0].thermal.mass_kg*cells[0].thermal.heat_capacity_j_per_kg_k;
    ThermalVoxel material{capacity/(8e-6),{2,0.5,0.5}};
    std::vector<std::vector<std::size_t>> regions(3);
    for(std::size_t i=0;i<24;++i) regions[(i%6)/2].push_back(i);
    ElectrothermalBatteryPack pack(cells,{},ThermalGrid3D(grid,{material}),regions);
    std::cout<<"time_s,voltage_v,cell_heat_w,cooling_w,peak_temperature_k,energy_residual_j\n";
    for(int n=0;n<60;++n) {
        const auto r=pack.step(0.5,1);
        std::cout<<r.electrical.time_s<<','<<r.electrical.terminal_voltage_v<<','<<r.thermal.supplied_power_w<<','
            <<r.thermal.outward_cooling_power_w<<','<<r.thermal.maximum_temperature_k<<','<<r.thermal.energy_balance_error_j<<'\n';
    }
}
