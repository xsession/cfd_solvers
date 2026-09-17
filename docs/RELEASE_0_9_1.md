# cfd_solvers 0.9.1

## Scope

v0.9.1 continues the CST-class expansion without changing the clean-room architecture. The release concentrates on the dependency chain needed for practical harmonic electromagnetics and adjacent multiphysics: complex sparse algebra, 3-D curl-conforming FEM, ports/Q/adaptivity, cable propagation, relativistic particles/wakes, and nonlinear low-frequency magnetics.

## 3-D frequency-domain Maxwell

The new `cfd::em::solve_driven_edge_maxwell_3d` path uses lowest-order first-kind Nedelec edge elements on the existing Tet4 mesh. Global edges are oriented by node index and local element signs map to that deterministic orientation. Boundary faces are detected topologically; every edge belonging to an exterior face is constrained to zero tangential electric field for the PEC baseline.

The complex curl-curl system is assembled directly into `ComplexCsrBuilder`. The shared complex sparse solve maps the system to the equivalent real block form and uses CSR + ILU(0)-GMRES. Indefinite systems that are a poor fit for ILU(0) retry with sparse unpreconditioned GMRES rather than silently falling back to dense algebra.

Supported excitation in this baseline:

- volume impressed current density `J(x)`;
- directed current excitation on an interior mesh edge, used as a minimal lumped/edge port seam.

The result reconstructs element-centroid electric and magnetic phasors.

### Validation

A manufactured PEC field

`E = y_hat sin(pi x) sin(pi z)`

is driven by the corresponding analytical impressed current. The regression checks sparse convergence and error reduction between Tet4 mesh resolutions. Edge-current excitation is independently verified to produce a non-zero field response.

## Wave ports and resonator Q

`rectangular_waveguide_modes` provides analytical rectangular PEC TE/TM transverse modes sorted by cutoff frequency. Propagating modes are normalized so the time-average forward Poynting power is exactly 1 W; evanescent modes use unit transverse-electric L2 normalization. The WR-90 TE10 cutoff and 1 W normalization are regression tested.

`resonator_quality_3d` computes:

- electric stored energy;
- magnetic stored energy;
- volume conductive/dielectric loss;
- good-conductor wall loss through surface resistance;
- `Q = omega W / P_loss`.

The regression verifies the expected `Q proportional sqrt(sigma_wall)` scaling when conductor loss dominates.

## Adaptive RF FEM

The first 3-D RF adaptivity loop is intentionally simple and inspectable:

1. reconstruct element fields;
2. compute interior-face tangential E/H jumps;
3. accumulate element indicators;
4. apply Dorfler bulk marking;
5. choose the longest edge of marked tetrahedra;
6. bisect the entire incident edge star so neighboring tetrahedra remain conforming.

This is a local refinement baseline, not a production anisotropic RF mesher.

## Rational RF reduced-order model

The existing adaptive log-frequency N-port sweep can now be compressed into a stable common-pole rational model

`H(s) = D + sum_k R_k (-p_k)/(s-p_k)`

with negative real poles logarithmically spanning the sampled band. Each matrix entry is fitted by complex least squares. A one-pole analytical response is recovered at frequencies not present in the training sweep.

## Cable/harness baseline

`solve_multiconductor_cable` solves the frequency-domain multiconductor telegrapher equations for full per-unit-length R/L/G/C matrices. A second-order `[1/1]` Pade propagation step is cascaded along the cable and the input currents are solved against a full complex matrix load termination.

Validation includes:

- a matched lossless line with analytical characteristic impedance and phase;
- a coupled two-conductor line producing crosstalk;
- frequency-dependent shield transfer impedance;
- shield-current and effective-height field-to-cable voltage coupling helpers.

## Relativistic particles, boundaries and wakes

The relativistic Boris update advances proper velocity `u=gamma v` and preserves gamma in the magnetic-only regression. Axis-aligned particle walls support absorption and specular reflection; absorbing walls can report expected secondary macro-particle weight through a smooth configurable yield curve.

The wake baseline provides:

- a causal resonator point-charge longitudinal wake;
- bunch line-charge convolution;
- wake-to-longitudinal-impedance transform.

The impedance regression peaks at the configured resonator frequency.

## Nonlinear magnetics, coils and force/torque

`Magnetostatics2D::solve_nonlinear` adds elementwise Picard integration of a B-dependent reluctivity law. `PiecewiseLinearBHCurve` provides a readable saturation baseline. `RectangularStrandedCoil2D` supplies current density and exports flux linkage/inductance to a circuit-facing scalar seam.

`maxwell_stress_boundary_force_2d` integrates the Maxwell stress tensor on a selected boundary patch. Uniform-field force and torque are checked analytically.

## Validation status

- normal CPU CTest matrix: **55/55 pass**;
- focused v0.9.1 ASan+UBSan executable: **pass** with leak detection;
- integration tracker: **488/638 (76.5%)**;
- Phase 12: **22/50 (44.0%)**.

## Explicit limitations

This release does not claim production CST parity. In particular, the following remain open:

- arbitrary cross-section FEM wave-port eigenmodes and port de-embedding;
- impedance/radiation/open boundaries for 3-D frequency-domain FEM;
- 3-D cavity eigenvalue extraction from the edge-FEM matrices;
- nonlinear eddy-current material integration and solid-conductor external-circuit coils;
- moving-band electrical-machine interfaces;
- RWG surface MoM, MLFMM and SBR/PO;
- eye/BER/PDN compliance workflows;
- Vay pushing, EM-PIC, collisions/plasma ionization and breakdown;
- heterogeneous anatomical SAR averaging and closed-loop EM/thermal tissue feedback.
