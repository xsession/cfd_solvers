# v0.6.0 release checkpoint

v0.6.0 broadens the unified C++20 solver framework beyond the v0.5.0 FVM/electrochemistry milestone while preserving the same clean-room and validation rules.

## Capability status

At this checkpoint `scripts/integration_status.py` reports 132/397 tracked capabilities complete (33.2%). Phase 0 project/release traceability is complete. The principal breadth advances are Elmer-class FEM at 22/45, openEMS-class FDTD at 8/26, and Optiland-class optics at 12/43.

## FEM

- Line2, Tri3, Quad4, Tet4, Hex8, Prism6 and Pyramid5 reference elements.
- Shape functions, reference gradients, Gaussian quadrature and isoparametric Jacobians.
- Tri3 assembled CSR and matrix-free Laplace actions.
- Mixed Dirichlet/Neumann/Robin scalar diffusion.
- 2-D Poisson/heat/linear elasticity, electrostatics and DC conduction.
- Axisymmetric elasticity with the `2*pi*r` measure and hoop strain.
- 3-D Tet4 Poisson with mesh-refinement convergence.

## FDTD

- Existing 1-D Yee baseline extended with heterogeneous epsilon/conductivity coefficients.
- Hard and soft sources.
- First-order Mur absorbing boundary.
- Reusable time and DFT monitors.
- 3-D Cartesian Maxwell baseline with PEC enclosure.

## Optics

- Sequential real-ray spherical/planar surface tracing and apertures.
- Paraxial first-order tracing and back focal distance.
- Sellmeier dispersion with small built-in material catalogue.
- Jones/Stokes primitives and dielectric Fresnel coefficients.
- Normal-incidence multilayer characteristic-matrix coating calculation.

## Reproducibility

The release gate requires GCC/OpenMP, Clang serial fallback, GCC serial, targeted ASan+UBSan smoke, `git diff --check`, a clean source-only rebuild, and rebuild/CTest after applying the exact v0.5.0->v0.6.0 patch to the v0.5.0 tree.

Hardware-only MPI/SYCL execution remains a separate validation item; no local MPI or AdaptiveCpp/GPU runtime is available in the development container.
