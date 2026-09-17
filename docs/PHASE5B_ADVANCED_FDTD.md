# v0.7.1 - advanced 1-D FDTD boundary/source/circuit baseline

This checkpoint advances three openEMS-class capability targets using independent C++20 implementations rather than translated upstream source. openEMS documents UPML/Mur boundaries, TFSF excitation and lumped RLC elements as core features; those capabilities are used as the feature map only.

## CPML

`Maxwell1D` adds a convolutional-PML boundary option alongside the existing matched-loss PML. Electric and magnetic derivatives use graded conductivity and kappa profiles plus auxiliary convolution memories. `cpml_alpha_fraction` is exposed for complex-frequency shifting and defaults to zero because the broadband normal-incidence regression is more accurate without an alpha shift; applications targeting evanescent/low-frequency absorption can tune it explicitly.

The outermost cells retain a PEC termination behind the absorbing layer. The CPML test measures attenuation relative to both initial and PEC residual energy rather than claiming it always outperforms Mur for every 1-D normal-incidence waveform.

## TFSF

`set_tfsf_source()` defines the left boundary of a +x total-field region. The correction uses this solver's `Ez/Hy` orientation and evaluates the magnetic incident field at both the Yee half-time and half-cell spatial offset. This staggering is required to suppress artificial leakage into the scattered-field side. The homogeneous regression reaches a scattered/total peak ratio of about `5.2e-6`.

The current API is intentionally one-dimensional and +x only; arbitrary 3-D plane-wave boxes remain future work.

## Lumped parallel RLC

`set_parallel_lumped_rlc()` maps a field-aligned cell voltage `V=E*l` and terminal current through a specified cross-sectional area. A parallel capacitor contributes an equivalent local permittivity, a parallel resistor contributes local conductivity, and a parallel inductor advances current density as an auxiliary state. Dispersive ADE and lumped models are rejected on the same cell in this baseline to avoid an unvalidated compound constitutive update.

Validation separately checks analytical inductor current ramp, resistor dissipation and capacitor numerical stability.

## Remaining openEMS-class gaps

The tracker still leaves 3-D CPML/TFSF, cylindrical grids, anisotropic media, NF2FF, SAR, HDF5, accelerator FDTD, distributed FDTD and shared geometry/material integration for later phases.
