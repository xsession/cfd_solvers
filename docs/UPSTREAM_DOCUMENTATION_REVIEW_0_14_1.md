# Upstream documentation review - v0.14.1 LBM completion wave

This is a clean-room capability/documentation review. The implementation in `cfd_solvers` is independent C++20/SYCL code; upstream source is not vendored or mechanically translated.

## FluidX3D capability surface

Reviewed:

- Project documentation/README: https://github.com/ProjectPhysX/FluidX3D

The current public documentation calls out several useful capability boundaries: Smagorinsky-Lilly LES, VOF/free-surface flow with PLIC curvature, local body-force fields, passive or two-way immersed-boundary particles, GPU voxelization, and 16-bit population/DDF storage while arithmetic remains higher precision.

### Decisions applied here

1. **LES is a collision-time correction, not a new global solver.** `EsotericPullSolver` can compute a local Smagorinsky relaxation time from the non-equilibrium stress tensor while retaining the existing streaming/layout path.
2. **Free surface is kept as a separate transport model.** `free_surface.*` owns fill fraction, fluid/interface/gas classification, normals/curvature and continuum-surface-force acceleration. The first implementation is a conservative periodic VOF baseline; it does **not** claim analytic PLIC parity.
3. **Particle coupling is explicit about one-way versus two-way behavior.** A regularized trilinear delta interpolates velocity to the particle and, in two-way mode, spreads the equal-and-opposite drag back into the new local LBM acceleration field.
4. **Memory precision is decoupled from arithmetic precision.** The compressed solver stores 16-bit encoded populations but reconstructs float values for collision/math. The codecs are independent IEEE half, BF16, shifted-half, and fixed-deviation experiments; they are not copies or claimed bit-compatible implementations of FluidX3D FP16S/FP16C.
5. **GPU geometry work stays behind the common geometry contract.** The optional SYCL voxelizer emits the same solid mask expected by the CPU path instead of introducing a GPU-only mesh representation.

## Curved-wall interpolation

Reviewed OpenLB public documentation/reference material:

- Bouzidi boundary API: https://www.openlb.net/DoxyGen/html/d3/d53/setBouzidiBoundary_8h.html
- Boundary module overview: https://www.openlb.net/DoxyGen/html/df/ddc/boundary3D_8h.html

The transferable design point is a link-wise wall distance and an interpolated bounce-back update with a half-way bounce-back fallback. `D2Q9Solver` now stores an optional per-link wall fraction and applies an independently written two-branch Bouzidi/Firdaouss/Lallemand interpolation. A wall fraction of 0.5 is regression-tested against the existing half-way result.

## Q-criterion

Reviewed NASA description of the velocity-gradient invariants and Q-criterion:

- NASA NTRS 20130001604: https://ntrs.nasa.gov/citations/20130001604

For incompressible flow, the implemented diagnostic uses

`Q = 0.5 * (||Omega||^2 - ||S||^2)`

with the symmetric strain tensor `S` and antisymmetric spin tensor `Omega` from central finite differences of the 3-D velocity field. Regression cases cover solid-body rotation and simple shear.

## Immersed-boundary interpolation/spreading

Reviewed the standard immersed-boundary architecture and recent regularized-delta literature, including:

- Gruninger & Griffith, Journal of Computational Physics 546 (2026), regularized delta interpolation/spreading.

The implementation deliberately starts with a compact trilinear tensor-product kernel and point-particle drag. It verifies discrete action/reaction conservation. It is **not** yet a resolved finite-size surface IB method, lubrication/contact model, or multi-particle collision engine.

## Cumulant / central-moment collision

Reviewed descriptions of central-moment and cumulant collision spaces, including:

- Coreixas et al., cross-platform many-core LBM collision model comparison: https://pmc.ncbi.nlm.nih.gov/articles/PMC8084255/
- cumulant-space multiphase solver description: https://pmc.ncbi.nlm.nih.gov/articles/PMC7911600/
- D2Q9 central-moment formulation background: https://pmc.ncbi.nlm.nih.gov/articles/PMC7516464/

The v0.14.1 `CumulantD2Q9Solver` is a low-Mach baseline. It transforms pulled populations to central moments/cumulants, relaxes second-order cumulants with the viscosity-controlled rate, relaxes higher cumulants toward factorized equilibrium, reconstructs raw moments, and maps them back to populations. Taylor-Green regression checks mass conservation and kinetic-energy decay. This is intentionally narrower than a production D3Q27 cumulant model with optimized forcing/boundaries.

## Validation boundary and remaining Phase-2 work

Executed on this host:

- `cfd-tests`: passed.
- `cfd-completion-tests`: passed.
- `cfd-v0140-hpc-runtime-tests`: passed.
- `cfd-v0141-lbm-les-qcriterion-tests`: passed.

The SYCL voxelization code was API/syntax-checked against the optional SYCL seam, but no SYCL accelerator runtime is installed on this packaging host. Real NVIDIA/AMD/Intel hardware CI therefore remains unchecked rather than being inferred from compilation.

After this wave, the deliberate Phase-2 gaps are hardware CI on all three GPU vendors, device-specific kernel fusion/autotuning, multiple GPUs per MPI rank, and adaptive domain repartitioning.
