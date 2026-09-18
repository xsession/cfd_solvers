# cfd_solvers v0.19.1

## Series battery pack workflow

`cfd/battery/pack.hpp` provides `SeriesBatteryPack`, a heterogeneous series string of existing SPMe cells. Each cell retains its own electrode, electrolyte, thermal and degradation configuration.

```cpp
#include "cfd/battery/pack.hpp"
using namespace cfd::battery;
std::vector<LithiumIonCellConfig> cells(3, default_graphite_nmc_config());
cells[0].positive.particle.initial_stoichiometry = 0.30;
PackControl control;
control.balancing_current_a = 0.05;
SeriesBatteryPack pack(cells, control);
auto result = pack.step(-0.1, 1.0); // 100 mA charging for one second
```

Positive terminal current means discharge. An enabled shunt adds its positive bleed current to its cell's terminal current. During charge, it reduces net charging of that cell; at rest it discharges the selected cell. At the start of each step, cells above both the minimum balancing voltage and the lowest-cell voltage plus the threshold receive the configured bleed current. Ties do not bleed. Zero balancing current disables control.

The shunts are ideal regulated dissipative current sinks, held constant during a step. They do not model a fixed resistor, PWM switching, active charge transfer or controller hysteresis. Their heat is reported separately and is not injected into the cell's lumped thermal model.

Pack voltage is the sum of loaded cell voltages. The implementation obeys:

- `I_cell[i] = I_pack + I_bleed[i]`;
- `sum(V_cell[i] * I_cell[i]) = P_terminal + P_bleed`;
- charge is the signed terminal current integral (not per-cell capacity);
- terminal and balancing energies use endpoint/backward-rectangle quadrature.

Build and run `./build/cfd-example-battery-pack` for a 60-second, three-cell charging example with CSV voltage, power and energy output.

## Limits and reset

Voltage and temperature limits are checked at the unloaded control sample, loaded start and trial endpoint. A violation throws `std::domain_error`; invalid step/control input throws `std::invalid_argument`. Every step uses trial cell copies, so a later-cell failure cannot leave earlier cells advanced. State and cumulative accounting are committed together. `reset()` resets all cells, clocks and integrals.

This is endpoint protection, not within-step event localization. Shorten timesteps near limits. The inherited SPM/SPMe reduced models, concentration clipping, illustrative OCV functions and degradation diagnostics retain their existing limitations. This is not a qualified battery-management controller. Parallel strings, temperature-dependent property feedback, 3-D pack heat flow and DFN/P2D are not implemented here.

## SPMe query correction

`SingleParticleElectrolyteModel::state()` now reports electrolyte polarization and actual concentration extrema. Previously a direct SPMe query inherited the SPM implementation and omitted those terms. The legacy base-class API remains nonvirtual; query SPMe through its concrete type. Grids smaller than three electrolyte cells are rejected during construction.

## Validation

New `cfd-v0191-pack-tests` covers series voltage parity against independently stepped cells, heterogeneous shunt-current parity, Kirchhoff power conservation, signed charging and rest balancing, charge and energy integrals, reset, second-cell cutoff rollback, endpoint thermal cutoff, nonfinite inputs and invalid grids. The original v0.19.0 battery regression is retained.

Phase 16A: 16/18. Overall tracker: 716/825 (86.8%). Full DFN/P2D and 3-D battery thermal coupling remain open.

### Results on this build host

- GCC Release/OpenMP build with Python bindings and examples: passed.
- CTest: **164/164 passed** (`OMP_NUM_THREADS=2`, four parallel tests).
- Focused pack test under AddressSanitizer + UBSan: passed with `ASAN_OPTIONS=detect_leaks=0`; LeakSanitizer is unsupported under this runner's ptrace environment.
- Three-cell charging CSV example: 60 steps completed.
- HDF5 was requested but not installed; CMake used its inline-XDMF fallback. Native HDF5 paths were not validated. MPI/SYCL/GPU paths were not exercised.

Logs are included under `docs/validation_v0191/`. The retained v0.19.0 benchmark reference was not regenerated or represented as new performance evidence.
