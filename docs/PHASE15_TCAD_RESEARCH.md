# Phase 15 semiconductor TCAD clean-room research

## Reference boundary

Phase 15 uses DEVSIM public documentation and standard semiconductor-device equations as a capability and validation reference. DEVSIM's public README describes a finite-volume TCAD simulator with user-defined PDEs, DC, small-signal AC, transient analysis, extended precision and 1-D/2-D/3-D simulation. The pinned reference for this study is `43b41ca845184c47e22b72d144db7e7db8509377` under Apache-2.0.

The cfd_solvers implementation is independent C++20 code. It does not translate DEVSIM classes or source layout.

## v0.18.0 numerical baseline

The first TCAD release deliberately starts with the classical non-degenerate 1-D drift-diffusion system:

- Poisson: `d2 psi / dx2 = -q(N_D-N_A+p-n)/epsilon`.
- Electron continuity using conventional electron current.
- Hole continuity using conventional hole current.
- Einstein relation `D = mu V_T`.
- Scharfetter-Gummel exponential fitting on each finite-volume edge.
- Ohmic contacts using local charge-neutral equilibrium carrier densities.
- A Schottky carrier-density boundary baseline based on a thermionic Boltzmann barrier.

The Bernoulli operator is evaluated with a series near zero and asymptotic branches for large positive/negative arguments. This avoids cancellation in low-field cells and overflow in strongly biased cells.

## Model helpers included

- arbitrary net-doping callback/arrays;
- PN and PIN profile generators;
- Caughey-Thomas doping mobility;
- high-field velocity-saturation mobility;
- SRH, Auger and radiative recombination;
- DC Gummel iteration;
- contact-voltage I-V sweeps;
- electric-field/current-density diagnostics;
- Joule-heating output for electrothermal coupling.

## Explicitly not claimed in v0.18.0

After v0.18.1, the following remain tracker-open: Fermi-Dirac/incomplete-ionization statistics, insulating/gate boundaries, BJT/MOS structures, MOSFET/IGBT device solvers, and compact-model extraction. The v0.18.1 small-signal path is a terminal/quasi-static linearization baseline, not a full distributed frequency-domain carrier Jacobian. The Schottky boundary is a carrier-density baseline and is not yet a full thermionic-emission/contact-current law.

## Validation strategy

v0.18.0 checks:

1. the Bernoulli detailed-balance identity;
2. zero recombination at mass-action equilibrium;
3. expected mobility degradation with doping and field;
4. PN built-in potential `2 V_T ln(N/ni)` for symmetric abrupt doping;
5. near-zero equilibrium terminal current;
6. exact uniform-bar Scharfetter-Gummel drift current `q mu n E`;
7. face-to-face current continuity;
8. Joule heating integrated against volume equals terminal electrical power;
9. Schottky Boltzmann boundary density;
10. PN/PIN majority-carrier/profile behavior.

## v0.18.1 transient and small-signal analysis

v0.18.1 adds three analysis capabilities while preserving the v0.18.0 finite-volume/Scharfetter-Gummel DC path.

### Backward-Euler transient

The electron and hole continuity equations add the implicit accumulation term `dx^2/(D dt)` to the diagonal while using the previous accepted carrier state on the right-hand side. Poisson remains quasi-static and is re-solved inside each transient Gummel iteration. The reported left-terminal transient current is the conductive current plus dielectric displacement current `epsilon A dE/dt`.

### Small-signal terminal linearization

The current AC baseline linearizes the converged nonlinear device around a DC operating point with symmetric finite differences. It reports

`Y(omega) = dI/dV + j omega dQ/dV`,

where `Q` is the electrostatic charge on the left terminal. This is a useful terminal impedance/C-V baseline, but it is intentionally **not** claimed to be a full distributed frequency-domain linearization of the coupled Poisson/electron/hole Jacobian.

### C-V workflow

C-V sweeps bias the left contact, solve DC, and evaluate `dQ/dV` using the same symmetric perturbation. A uniform charge-neutral bar reproduces the geometric capacitance `epsilon A/L`.

### Temperature dependence

Optional temperature scaling adds a silicon Varshni bandgap law, intrinsic-density scaling, power-law electron/hole mobility scaling, and density-of-states scaling relative to the construction-time material reference state. Existing 300 K behavior is unchanged.

### v0.18.1 validation

The new regression verifies analytic conductance `q mu n A/L`, geometric capacitance `epsilon A/L`, admittance/impedance reciprocity, backward-Euler voltage-step displacement current `C DeltaV/Delta t`, conductive step current `G DeltaV`, C-V consistency, and expected temperature trends for bandgap, intrinsic density, mobility and density of states.

## v0.18.2 MOS/gate physics and compact-model bridge

v0.18.2 adds blocking semiconductor boundaries and device structures while keeping the original Scharfetter-Gummel transport path intact.

### Insulating and gate boundaries

`ContactType::insulating` enforces zero electrostatic normal derivative and zero electron/hole Scharfetter-Gummel flux. `ContactType::gate` uses the same zero-carrier-flux condition but applies an oxide-capacitance Robin condition,

`(epsilon_s/dx) (psi_surface-psi_adjacent) = C_ox (V_gate - V_fb - psi_surface)`.

The electron/hole boundary unknown is eliminated from the first/last continuity row using the exact zero-Scharfetter-Gummel-flux ratio, rather than approximating zero flux as equal carrier densities.

### Fermi statistics and incomplete ionization

The common TCAD helpers now include the normalized complete Fermi-Dirac integral `F_{1/2}(eta)`, its numerical inverse, carrier-density helpers `N_c F_{1/2}(eta_n)` / `N_v F_{1/2}(eta_p)`, and donor/acceptor ionized-density functions with binding energy and degeneracy. This is a reusable statistics/ionization framework; the v0.18.2 1-D drift-diffusion state equations still default to the existing Boltzmann carrier closure unless a higher-level model explicitly uses the Fermi helpers.

### MOS capacitor

`MosCapacitor1D` solves the standard quasi-static p-type MOS surface-potential equation using the exact Boltzmann semiconductor charge integral and an oxide capacitance. It provides gate charge, surface potential, finite-difference C-V sweeps, accumulation, depletion, and quasi-static inversion behavior. This standalone formulation is used for robust MOS electrostatics validation while the lower-level `SemiconductorDevice1D` gate boundary remains available for general transport problems.

### Long-channel MOSFET

`LongChannelMosfet` is a gradual-channel charge-sheet drift-diffusion model. The drain current is the analytic integral of

`I_D = W mu_n Q_inv(V) dV/dx`, with `Q_inv=-C_ox(V_GS-V_T-V)`,

including pinch-off, a simple channel-length-modulation factor, optional body effect, drain/source reversal, operating-point derivatives, and I-V table generation. It is deliberately labeled a long-channel compact drift-diffusion baseline and does not claim 2-D Poisson/channel-field parity.

### TCAD to SPICE

`make_long_channel_mosfet_spice_evaluator()` exports the model through the circuit solver's generic four-terminal `StaticDeviceEvaluator` interface `[D,G,S,B]`. Terminal currents satisfy KCL and the centered finite-difference Jacobian is stamped directly by `Circuit::add_static_device()`. This is the first explicit TCAD-to-SPICE compact-model extraction/interface seam.

### v0.18.2 validation

The new regression checks the Boltzmann limit and inverse round-trip of `F_{1/2}`, incomplete-ionization trends, exact zero carrier current at gate/insulating faces, MOS accumulation/depletion/inversion C-V behavior, analytical long-channel linear/saturation current, drain/source reversal, compact-device KCL/Jacobian conservation, and convergence of the extracted evaluator inside the actual SPICE Newton solver.

After v0.18.2, only the BJT baseline and IGBT/power-device foundation remain open in Phase 15.

## v0.18.3 bipolar and power-device completion

The final machine-counted Phase-15 wave adds two intentionally scoped device foundations. Public capability references include DEVSIM's BJT publication/example repository and the numerical BJT device family documented by CIDER/ngspice. The implementation remains independent C++20 and does not translate either source tree.

### BJT structure/model baseline

`BipolarJunctionTransistor1D` represents an abrupt emitter/base/collector stack with explicit region widths and majority doping. Minority-carrier diffusion coefficients follow the Einstein relation, diffusion lengths follow `L=sqrt(D tau)`, and the finite-base transport factor is `alpha_T=1/cosh(W_B/L_n)`. Emitter/collector injection efficiencies combine base electron injection with minority-hole injection in the adjacent n/p regions. Forward/reverse common-base gains and saturation currents are then used in a reciprocal Ebers-Moll terminal relation satisfying `alpha_F I_ES = alpha_R I_CS`.

This is a physics-parameterized 1-D charge-control baseline, not a field-resolved three-contact BJT PDE solve. It nevertheless exposes the spatial doping profile, built-in junction voltages, transport/injection diagnostics, NPN/PNP polarity, temperature scaling and a three-terminal `[C,B,E]` SPICE evaluator.

### IGBT / power-device foundation

`IgbtPowerDevice` composes the existing long-channel nMOS drift-diffusion channel with a low-gain PNP transport path. The MOS channel provides the gate-controlled injection current; a bounded bipolar gain represents conductivity modulation; drift-region resistance decreases with injected channel current; electron mobility and drift resistance respond to junction temperature. The transient state adds an exponential bipolar tail-current store and a lumped thermal RC.

The model intentionally does not claim Hefner-model, high-injection ambipolar PDE, latch-up or lifetime-profile parity. It is a stable power-device foundation with explicit hooks for those later refinements. The static `[C,G,E]` evaluator can be stamped directly into the existing SPICE Newton solver.

### v0.18.3 validation

The regression checks abrupt NPN doping polarity, reciprocity, diffusion-derived forward/reverse gain, NPN/PNP current symmetry, BJT KCL/Jacobian conservation and SPICE convergence. IGBT checks cover gate cutoff, MOS-driven bipolar current enhancement, conductivity-modulated drift resistance, temperature-dependent channel mobility, self-heating, turn-off tail decay, evaluator KCL and SPICE convergence. The v0.18.0/v0.18.1/v0.18.2 TCAD regressions remain passing.

After v0.18.3, Phase 15 is 30/30 machine-counted capabilities complete.
