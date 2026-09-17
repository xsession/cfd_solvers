# Phase 5 - FDTD Baseline Expansion

The openEMS-class branch now contains a 3-D Cartesian Yee baseline plus a more complete 1-D validation harness.

New validated primitives:

- per-cell dielectric permittivity and conductivity preprocessing;
- Crank-Nicolson-style conductive electric-field coefficient update;
- soft and hard electric-field sources;
- first-order Mur absorbing boundaries;
- reusable time probes and trapezoidal DFT monitors;
- 3-D PEC enclosure stability/CFL smoke test.

The Mur validation launches the same Gaussian pulse in PEC and absorbing domains and requires the residual Mur energy after exit to be below 1e-4 of the PEC reflected-energy level.

UPML/CPML, dispersive Drude/Lorentz/Debye materials, ports/S-parameters, NF2FF, SAR and accelerator kernels remain tracked work.
