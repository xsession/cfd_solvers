# cfd_solvers v0.8.1

v0.8.1 is a focused analysis-depth checkpoint on top of the broad v0.8.0 RF/SPICE integration.

The authoritative tracker reports **433/577 capabilities (75.0%)** validated: Phase 10 RF/antenna/microwave is **31/63 (49.2%)** and Phase 11 SPICE/circuit/compact-model simulation is **64/96 (66.7%)**.
The release GCC/OpenMP CPU matrix passes **47/47 CTest targets**.

## Added

- adaptive variable-step circuit transient integration using backward-Euler step doubling and local-truncation-error control;
- complex rational-fit pole/zero extraction for low-order small-signal SISO transfer functions;
- periodic steady-state shooting-by-settling with phase-aligned cycle convergence;
- coupled parallel-wire thin-wire MoM with independent segmentation, delta-gap feeds and lumped series loads;
- a dedicated focused regression target covering analytical RC and reciprocal two-wire cases.

## Deliberate limitations

The pole-zero path is a fitted baseline rather than an exact generalized-eigenvalue MNA/DAE solve. PSS is transient settling rather than harmonic balance/Newton shooting. The wire solver is restricted to parallel z-directed conductors and does not yet provide NEC-grade arbitrary junction geometry, finite/conductive ground or canonical benchmark convergence claims.
