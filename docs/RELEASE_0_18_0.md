# cfd_solvers v0.18.0

## Phase 15 semiconductor TCAD baseline

v0.18.0 starts the new semiconductor-device solver family with an independent 1-D finite-volume Poisson + electron/hole drift-diffusion implementation.

### Added

- `cfd::tcad::SemiconductorDevice1D`.
- stable Bernoulli/Scharfetter-Gummel flux evaluation;
- net-doping fields and PN/PIN helpers;
- nondegenerate Boltzmann equilibrium initialization;
- ohmic and Schottky carrier boundaries;
- SRH/Auger/radiative recombination helpers;
- doping/high-field mobility helpers;
- DC Gummel operating point and contact-voltage sweep;
- electron/hole/total current diagnostics;
- electric-field and Joule-heating output for electrothermal coupling.

### Scope boundary

This release does not yet claim transient or small-signal drift-diffusion, C-V, MOS/BJT/IGBT device solvers, Fermi-Dirac statistics, insulating/gate boundaries, or TCAD-to-SPICE extraction.

### Validation

- default Python/HDF5 matrix: **143/143 CTest targets passed**;
- Python C ABI reports **0.18.0**;
- strict TCAD source/test compilation passes `-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror`;
- symmetric abrupt-PN built-in voltage matches `2 V_T ln(N/ni)`;
- uniform n-type Scharfetter-Gummel current matches `q mu n E`;
- Joule-heating volume integral matches terminal electrical power.

### Tracker

Phase 15: **20/30 = 66.7%**. Overall: **690/807 = 85.5%**.
