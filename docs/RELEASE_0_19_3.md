# cfd_solvers v0.19.3

## DFN/P2D distributed porous-electrode model

`cfd/battery/dfn.hpp` adds `DoyleFullerNewmanModel`, the distributed limit of the
single-particle battery models. Instead of collapsing each electrode to one
spherical particle and the electrolyte to averaged states, the model resolves the
through-cell axis x in [0, L] (L = L_neg + L_sep + L_pos) into one electrolyte
node per electrode node plus the separator. Each electrode node owns one
`spherical_diffusion_particle` whose surface reaction flux is set from a
local Butler-Volmer-consistent partition of the imposed terminal current, so the
solid concentration, the electrolyte concentration and the solid/electrolyte
potentials are all fields along the cell.

Per time step:

1. **Reaction partition.** The full terminal current is distributed across the
   active nodes of each electrode, weighted by local exchange-current activity
   `i0 * 2*sqrt(x(1-x))` evaluated at the previous-step surface stoichiometry.
   Each electrode carries the full current (SPMe convention), which keeps the
   global lithium balance closed.
2. **Solid diffusion.** Each per-node spherical particle advances with the same
   conservative implicit spherical finite-volume diffusion used by SPM/SPMe; its
   outward surface molar flux is the local partitioned reaction current divided
   by the active surface density and Faraday's constant.
3. **Electrolyte concentration.** An implicit finite-volume diffusion solve
   (backward Euler, harmonic-mean face conductances) advances the electrolyte
   concentration with Neumann (zero-flux) boundaries and the local Butler-Volmer
   source `(1-t+)*i_k/F`. The scheme is conservative for the electrolyte lithium
   mole count; the full-cell balance is closed by the electrode partition.
4. **Potentials.** The electrolyte potential is reconstructed from the running
   sum of the local areal reaction current (zero at both collectors, peaking in
   the separator) integrated through the Bruggeman conductivity. The solid
   potential integrates the local reaction-driven solid current from the negative
   collector and is flat in the separator.
5. **Outputs.** `DfnStepResult` reports the area-averaged OCV and reaction
   overpotential, the distributed electrolyte ohmic plus Nernst concentration
   polarization, the distributed solid ohmic drop, terminal voltage, the
   electrode-end surface stoichiometries and the electrolyte concentration
   extrema. `total_lithium_mol()` sums the particle and electrolyte inventories
   for conservation checks.

The model is isothermal (the cell temperature tracks the initial condition) and
CPU-only for this release; no GPU kernel or 3-D thermal coupling is claimed for
the DFN state. Reversible/entropic heating is omitted, consistent with the rest of
the battery family.

## Example

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/cfd-solve battery-dfn
```

`battery-dfn` runs a 0.5 A discharge and reports the open-circuit voltage, loaded
terminal voltage, reaction overpotential, electrode-end surface stoichiometries
and the conserved total lithium inventory.

Minimal usage:

```cpp
#include "cfd/battery/dfn.hpp"
using namespace cfd::battery;
auto cfg = default_graphite_nmc_dfn_config();
DoyleFullerNewmanModel cell(cfg);
auto open = cell.state();            // OCV at rest
for (int i = 0; i < 120; ++i)
    (void)cell.step(0.5, 1.0);       // 0.5 A discharge, 1 s steps
double lithium = cell.total_lithium_mol();
```

## Validation and limits

`cfd-v0193-dfn-tests` validates: open-circuit voltage in a physical range and
terminal voltage equal to OCV at zero current; discharge depleting the negative
and lithiating the positive electrode with loaded voltage below OCV and a positive
Butler-Volmer overpotential; global lithium conservation when the electrolyte
source is removed (transference number 1); positive and distributed electrolyte
concentration; distributed solid and electrolyte potential drops under load;
charge reversing the stoichiometry and raising the terminal voltage above OCV; and
rejection of an invalid electrode node count. A `battery-dfn` cfd-solve smoke case
is wired into CTest.

Phase 16A is now **18/18** and the tracker is **718/825 = 87.0%**.

## Docker deployment system

A multi-stage `Dockerfile` builds the C++20 solvers (Release, OpenMP, optional
MPI/HDF5) in a builder stage and ships the binaries plus the source tree in a
slim non-root runtime stage. `docker-compose.yml` exposes:

- `cfd` — run any `cfd-solve` case (`docker compose run cfd lbm-d3q19-cpu`),
- `cfd-test` — the full CTest suite (`--profile test`),
- `cfd-mpi` — the distributed d3q19 case under `mpirun` (`--profile mpi`, MPI build).

`cfd-docker.bat` + `cfd-docker.ps1` provide a no-bash Windows helper:
`cfd-docker.bat build | run [case] [--threads N] | bench | test | mpi [--ranks N] |
list | ps | clean`. The build args (`CMAKE_BUILD_TYPE`, `CFD_ENABLE_OPENMP`,
`CFD_ENABLE_MPI`, `CFD_ENABLE_NATIVE_ARCH`, `CFD_ENABLE_HDF5`) are configurable on
both the image build and the compose file.
