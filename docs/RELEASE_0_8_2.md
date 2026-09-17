# cfd_solvers v0.8.2

v0.8.2 deepens the v0.8.1 circuit/RF analysis work without changing the clean-room scope. The tracker is **439/580 (75.7%)** overall, with RF Phase 10 at **32/64 (50.0%)** and SPICE Phase 11 at **69/98 (70.4%)**.

## Added

- adaptive variable-step BDF2/Gear integration with exact unequal-step coefficients;
- Richardson full-step/two-half-step LTE control appropriate for the second-order method;
- arbitrary-orientation mutually coupled straight-wire MoM in 3-D free space;
- a dyadic Green-kernel implementation projected onto source/observation wire tangents;
- ideal lossless transformer convenience stamping with a reflected-load regression;
- lossy sampled TEM-line attenuation validation and reconciled microstrip/stripline/CPW/TE10 circuit-element tracking;
- rotational-invariance and second-order RC refinement regressions.

## Compatibility

The original `solve_parallel_thin_wires()` and `solve_center_fed_thin_wire()` APIs remain available. The parallel implementation delegates to the generalized orientation-independent solver. Adaptive transient defaults remain backward Euler; callers opt into BDF2 with `AdaptiveTransientConfig::method`.

## Deliberate limits

Adaptive order selection above BDF2, exact descriptor-system pole-zero extraction, harmonic balance, connected/bent wire junction basis functions, finite/conductive ground and NEC-grade canonical convergence remain open.
