# cfd_solvers v0.18.2

v0.18.2 advances Phase 15 semiconductor TCAD from diode/transient analysis into MOS electrostatics, gate boundaries and circuit extraction.

## Added

- `ContactType::insulating` and `ContactType::gate` with zero carrier flux.
- Oxide-capacitance Robin electrostatic gate boundary and gate-charge diagnostic.
- Fermi-Dirac `F_{1/2}` / inverse helpers and incomplete donor/acceptor ionization functions.
- `MosCapacitor1D` quasi-static surface-potential, gate-charge and C-V workflow.
- `LongChannelMosfet` gradual-channel charge-sheet drift-diffusion current model with body effect, channel-length modulation, reversal symmetry and I-V sweeps.
- Four-terminal `[D,G,S,B]` TCAD-derived SPICE `StaticDeviceEvaluator`.
- v0.18.2 regression coverage for gate blocking, MOS C-V, MOSFET current and SPICE integration.

## Scope

The MOSFET implementation is a long-channel charge-sheet drift-diffusion baseline, not a 2-D field-resolved MOSFET. Fermi/incomplete-ionization helpers are available as a framework, while the base 1-D transport solver still uses the Boltzmann closure by default. BJT and IGBT/power-device structures remain Phase-15 gaps.

## Tracker

- Phase 15: **28/30 = 93.3%**.
- Overall: **698/807 = 86.5%**.

## Validation

- Complete default matrix: **145/145 CTest targets passed**.
- Python ABI: **0.18.2**.
- HDF5-enabled default configuration: passed.
- Strict TCAD/MOS compile: `-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror`.
- v0.18.0 and v0.18.1 TCAD regressions remain passing.

