# v0.8.1 adaptive transient, PSS, pole-zero and coupled-wire MoM

This checkpoint deliberately closes several narrow analysis baselines without claiming production SPICE or NEC equivalence.

## Adaptive transient

`Circuit::transient_adaptive()` uses backward Euler with Richardson step doubling. Every candidate interval is solved once at the full step and twice at half step. The difference at the common endpoint forms a normalized local-truncation-error estimate over node voltages. Accepted states use the two-half-step solution, so accepted samples remain more accurate than the error-estimator path. Step growth/shrinkage uses a safety factor and configurable clamps, and the final interval is shortened so the solution lands exactly on the requested stop time.

The implementation supports the same current transient device set as the fixed-step circuit solver: R/C/L, mutual inductance, independent and controlled sources, diode, MOS Level-1, BJT, JFET, switches and memoryless static-device callbacks. Dynamic compact-model state callbacks and higher-order variable-step BDF/Gear formulas remain open.

## Pole-zero baseline

`pole_zero_analysis()` samples the configured small-signal AC transfer function on a logarithmic frequency grid and fits

`H(s) = B(s) / A(s)`

with a complex least-squares rational model in a normalized frequency variable. Polynomial roots are then mapped back to the continuous-time s plane and reported in rad/s. The caller selects numerator and denominator orders; the returned relative RMS fit error is a mandatory quality indicator.

This baseline is useful for low-order linearized circuits and regression work. It is not yet the exact descriptor-system/generalized-eigenvalue pole-zero algorithm required for large sparse MNA/DAE systems.

## Periodic steady state

`periodic_steady_state()` runs a phase-aligned fixed-step transient over successive periods and compares corresponding phase samples in consecutive cycles. Convergence is declared only when every non-ground node satisfies the configured normalized cycle-to-cycle tolerance. The selected converged period is returned with time shifted to `[0, period]` for Fourier/THD post-processing.

This is a robust settling/shooting baseline, not harmonic balance. Strongly nonlinear or very high-Q circuits can require many periods and will ultimately need Newton shooting with a state-transition Jacobian and/or harmonic-balance infrastructure.

## Coupled parallel-wire MoM

`solve_parallel_thin_wires()` lifts the previous single straight dipole reference into a coupled system of multiple parallel z-directed wires. Each wire may have its own length, radius and segmentation. The solver assembles all self and mutual Pocklington-kernel terms into one dense complex matrix and supports multiple delta-gap feeds plus lumped series loads.

The regression checks nonzero induced current, reciprocal symmetry under swapping two identical wires and lower feed current after adding a resistive series load. The original `solve_center_fed_thin_wire()` API is now a one-wire wrapper around the generalized implementation.

Bent wires, junction continuity basis functions, finite/conductive ground, image theory, transmission-line wire sections and NEC-grade canonical convergence remain open and are still unchecked in the tracker.

## v0.8.2 continuation: unequal-step Gear/BDF2 and oriented wires

`AdaptiveTransientConfig::method = TransientMethod::bdf2` now upgrades the adaptive path after one accepted backward-Euler bootstrap interval. The BDF2 derivative uses the exact coefficients for the current and previous unequal time steps. LTE is estimated by full-step versus two-half-step Richardson comparison with the second-order divisor and cubic controller exponent. This closes adaptive BDF2/Gear, but does not claim automatic BDF3+ order selection.

The RF reference now also exposes `solve_oriented_thin_wires()`. Each straight PEC wire is defined by arbitrary 3-D endpoints. Mutual/self interactions use the free-space dyadic Green function projected onto source and observation tangents, while the existing parallel-wire API delegates to the new solver for compatibility. Disjoint wires are enforced deliberately: touching/crossing conductors still require junction-continuity basis functions and remain open.

### RF-circuit completion in the same checkpoint

`Circuit::add_ideal_transformer()` composes an ideal voltage-ratio constraint with the corresponding current relation, preserving instantaneous power in the controlled-source MNA baseline. The regression drives the primary from a current source and verifies both the turns ratio and the reflected resistive load.

The existing sampled `TemTransmissionLine` circuit bridge is now explicitly validated with non-zero attenuation: a matched line recovers `|S21| = exp(-alpha L)` while maintaining negligible reflection. Existing tested microstrip, stripline, coplanar-waveguide and rectangular-waveguide TE10 sampled elements are now reflected accurately in the tracker.
