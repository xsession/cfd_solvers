# Phase 16A battery cell physics research and clean-room implementation notes

## Scope

Phase 16A introduces lithium-ion cell models on top of the existing electrochemistry, nonlinear-solver, thermal and impedance infrastructure. PyBaMM is used as a capability and validation reference only; the C++20 implementation is independently written from standard battery equations and public literature.

Pinned public reference:

- PyBaMM: `71d9ff424f876debe67132b86f900d77e511fe8e`
- License tracked in `THIRD_PARTY.md` / `docs/upstreams.json`.

## v0.19.0 model hierarchy

### Spherical solid diffusion

`SphericalDiffusionParticle` solves conservative radial diffusion in a sphere with shell volumes and surface-area fluxes. A prescribed outward surface molar flux changes the average concentration by the exact control-volume balance. The surface concentration is reconstructed with a radial-gradient correction for reaction kinetics.

### Single Particle Model (SPM)

The SPM uses one representative spherical particle per electrode. Applied current maps to opposite-sign surface lithium fluxes in the negative and positive particles. Terminal voltage combines electrode open-circuit potentials, symmetric Butler-Volmer reaction polarization, solid/contact ohmic loss and the cell thermal state.

### Single Particle Model with Electrolyte (SPMe)

The SPMe extends the SPM with a 1-D through-cell electrolyte concentration field, region-dependent porosity/effective diffusivity, electrolyte ohmic/concentration polarization and reconstructed electrolyte/solid potential profiles.

This is a practical reduced-order electrolyte baseline; it is not the full distributed porous-electrode DFN/P2D system.

## Thermal models

Two levels are present:

1. lumped cell temperature driven by irreversible electrochemical loss and convective cooling;
2. `ThermalSlab1D`, a through-cell finite-volume conduction model with volumetric heat and surface convection.

A future 3-D pack thermal field should reuse the existing general FEM/FVM thermal infrastructure rather than duplicate it here.

## Degradation states

The v0.19.0 baseline tracks independent reduced-order state variables for:

- SEI thickness growth;
- plated lithium inventory during charging;
- loss of active material;
- concentration-gradient-driven particle cracking/damage.

These are intentionally compact engineering state models, not full coupled degradation PDEs.

## Analysis workflows

- Battery EIS uses a Randles-style charge-transfer/double-layer network plus finite-length diffusion contribution.
- Drive-cycle simulation accepts piecewise-constant current segments, limits the maximum integration step, and accumulates discharged capacity and electrical energy.

## Explicit remaining Phase 16A gaps

- full Doyle-Fuller-Newman / P2D porous-electrode model;
- 3-D cell/pack thermal coupling;
- multi-cell pack-control / balancing interface.

Those stay open in the machine-counted tracker.

## Validation strategy

v0.19.0 regressions check:

- exact spherical-particle lithium conservation against prescribed surface flux;
- discharge stoichiometry direction in both electrodes;
- loaded voltage below OCV during discharge;
- SPMe electrolyte concentration/potential gradients;
- analytic 1-D uniform thermal energy balance;
- monotonic degradation-state evolution under an accelerated test configuration;
- low/high-frequency EIS ordering;
- drive-cycle sample count, capacity and finite energy accounting.

## Performance baseline

Battery usage is represented in `examples/phase16a_battery_drive_cycle.cpp`. It runs an SPMe discharge/rest/charge profile and emits the same `CFD_BENCH` timing record used by the other major solver-family examples.
