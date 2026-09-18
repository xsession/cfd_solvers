# cfd_solvers v0.19.2

## 3-D thermal field

`cfd/battery/thermal3d.hpp` adds `ThermalGrid3D` and `ElectrothermalBatteryPack`.
The grid uses regular, cell-centred finite volumes, x-fastest indexing, voxel heat capacities and diagonal directional conductivities. Heat sources are **watts per voxel**, not W/m³. A single material is broadcast; alternatively provide one material per voxel.

For each internal face, the conductance is the face area divided by the sum of the two half-cell conduction resistances. The face heat flux enters adjacent cells with opposite signs. Boundary conductance includes both the half-cell conduction resistance and the convective resistance. Zero convection insulates a boundary; zero normal conductivity blocks transfer through that face.

Backward Euler solves:

`C_i * (T_new_i - T_old_i) = dt * (Q_i - sum_j G_ij*(T_new_i-T_new_j) - G_boundary_i*T_new_i + ambient_source_i)`.

The symmetric positive-definite system uses the shared matrix-free conjugate-gradient implementation. Solving for the temperature increment avoids letting a roughly 300 K baseline dominate the residual target. Failed convergence or invalid output leaves the old field intact. Controls expose iteration limit and tolerance; diagnostics include supplied power, outward boundary cooling, stored-energy change, energy-balance error, maximum temperature and iteration count. Cooling power is negative when the environment heats the grid.

## Pack coupling

Each SPMe cell maps to a nonempty, disjoint list of voxel indices. Unassigned voxels can represent passive pack material. The summed mapped heat capacity must agree with the corresponding cell's configured `mass_kg * heat_capacity_j_per_kg_k` within relative 1e-8. This prevents silently changing the cell's thermal inertia during coupling.

Each timestep:

1. Advance the series SPMe string at the previous mapped temperatures with its internal lumped thermal update disabled.
2. Distribute each cell's irreversible heat uniformly by mapped voxel volume, conserving total watts.
3. Solve implicit 3-D conduction and boundary cooling.
4. Check the maximum temperature across the entire grid, including passive material, against the pack limit.
5. Feed heat-capacity-weighted region temperatures back into the cells and commit both systems together.

A cutoff, invalid input or failed thermal solve rolls back the electrical state, particle inventories, clocks, cumulative terminal accounting and temperature field. Returned electrical voltages/powers belong to the electrical substep; returned cell temperatures are the thermal endpoint. The next electrical step uses those new temperatures.

This is first-order, explicitly partitioned coupling, not a monolithic electrothermal solve. Only existing SPMe temperature dependencies receive feedback; no new temperature-dependent transport/material laws are introduced. The grid supplies initial temperatures and cooling boundaries; cell lumped ambient/cooling parameters are unused under external ownership. Shunt dissipation remains external and separately reported; this release does not place balancing resistors in the thermal grid. The underlying model still omits reversible/entropic heating and retains its previously documented reduced-model limitations.

The standalone `set_external_temperature()` and series `set_external_temperatures()` APIs enable external thermal ownership. Calling the existing cell/series `reset()` restores internal lumped heating. Reconstruct an `ElectrothermalBatteryPack` to restart its coupled initial condition.

## Example

Build and run:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/cfd-example-battery-thermal > battery_thermal.csv
```

The example uses three series cells, each mapped to eight voxels in a 6x2x2 grid, directional conductivity, six convective faces and increased contact resistance in one cell. It reports voltage, heat, cooling, peak temperature and energy residual over a 60-second discharge.

Minimal two-cell setup:

```cpp
#include "cfd/battery/thermal3d.hpp"
using namespace cfd::battery;
std::vector<LithiumIonCellConfig> cells(2, default_graphite_nmc_config());
ThermalGridConfig grid;
grid.cells = {2, 1, 1};
grid.spacing_m = {0.01, 0.01, 0.01};
const double capacity = cells[0].thermal.mass_kg
                      * cells[0].thermal.heat_capacity_j_per_kg_k;
ThermalVoxel material{capacity / 1e-6, {1, 1, 1}};
ElectrothermalBatteryPack pack(cells, {}, ThermalGrid3D(grid, {material}), {{0}, {1}});
auto result = pack.step(0.2, 1.0);
```

## Validation and limits

`cfd-v0192-thermal-tests` validates a three-axis anisotropic discrete Neumann eigenmode, heterogeneous heat conservation, an analytical one-volume convective cooling solution including half-cell resistance, independent lumped-cell heating equivalence, conservative heat deposition, temperature feedback, thermal ownership reset, hotspot rollback, invalid mappings/materials and nonconvergent-solve rollback.

The CPU field supports orthogonal voxels and diagonal anisotropy only. No adaptive/curved mesh, radiation, fluid coolant dynamics, contact resistance layer, full tensor conductivity or GPU kernel is claimed. Pack limits are checked at sampled endpoints rather than by within-step event localization. Full DFN/P2D remains the sole open Phase 16A item.

Tracker: **717/825 (86.9%)**, Phase 16A **17/18 (94.4%)**.

### Build-host results

- Release/OpenMP build, Python bindings and examples: passed.
- CTest: all **165 targets passed** across the initial run (162 passes, three launch failures) and focused rerun (3/3 passes). Three cached executables had missing permissions or empty contents; relinking those targets resolved all launch failures.
- Focused AddressSanitizer/UBSan thermal regression: passed. Leak detection was disabled because this runner does not support LeakSanitizer under ptrace.
- Three-cell 60-second thermal example: passed; final peak 298.238 K, final energy residual approximately -6.99e-14 J.
- Native HDF5 was unavailable; inline XDMF fallback configured. HDF5, MPI and SYCL/GPU paths were not validated here.

Logs and example output are retained in `docs/validation_v0192/`. Existing historical benchmark references were not changed.
