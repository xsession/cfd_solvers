# cfd_solvers v0.18.3

v0.18.3 completes Phase 15 semiconductor TCAD with bipolar and MOS-gated bipolar power-device foundations.

## Added

- `BipolarJunctionTransistor1D` abrupt NPN/PNP emitter-base-collector structure.
- Diffusion-length/base-transport/injection-efficiency parameter derivation and reciprocal Ebers-Moll terminal transport.
- BJT built-in-junction, beta, saturation-current, power and doping-profile diagnostics.
- Three-terminal `[C,B,E]` TCAD-derived SPICE static evaluator.
- `IgbtPowerDevice` MOS-gated bipolar foundation with conductivity-modulated drift resistance.
- Temperature-scaled channel mobility and drift resistance.
- Turn-off bipolar tail-current storage and lumped thermal-RC self-heating.
- Three-terminal `[C,G,E]` IGBT SPICE evaluator.
- v0.18.3 regression coverage plus v0.18.0-v0.18.2 compatibility checks.

## Scope

The BJT is a 1-D charge-control/Ebers-Moll baseline whose parameters are derived from semiconductor geometry/materials; it is not a 2-D/3-D field-resolved three-contact drift-diffusion solve. The IGBT is a compositional power-device foundation, not a full Hefner/high-injection ambipolar or latch-up model.

## Tracker

- Phase 15: **30/30 = 100%**.
- Overall: **700/807 = 86.7%**.

## Validation

- Complete default matrix: **146/146 CTest targets passed**.
- Python ABI: **0.18.3**.
- HDF5-enabled default configuration: passed.
- Strict BJT/IGBT compile: `-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror`.
- v0.18.0, v0.18.1 and v0.18.2 TCAD regressions remain passing.
