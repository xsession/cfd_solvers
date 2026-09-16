# Phase 6 - Optics Baseline Expansion

The Optiland-class branch now validates both geometric and first-order optics rather than only vector Snell/reflection helpers.

Implemented and tested:

- sequential spherical/planar real-ray tracing with circular vignetting and CPU batch execution;
- RMS spot-radius analysis;
- paraxial surface tracing and back-focal-distance calculation;
- Sellmeier material model and a small built-in N-BK7/FusedSilica catalog;
- Jones-to-Stokes polarization conversion;
- Fresnel dielectric interface coefficients including total internal reflection;
- normal-incidence multilayer characteristic-matrix reflectance;
- quarter-wave anti-reflection coating validation.

The next optics blocks are coordinate breaks, conic/aspheric/freeform surfaces, general apertures, wavelength/field containers integrated into tracing, birefringence/coating polarization, non-sequential/ghost paths, diffraction and optimization/tolerancing.
