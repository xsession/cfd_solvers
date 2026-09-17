# Release 0.15.2 - transported RANS, Reynolds stress and SST-DES

## Highlights

- Adds a shared implicit finite-volume turbulence-scalar transport path with variable diffusivity, first-order upwind convection, fixed-value/zero-gradient patch handling, ILU0-preconditioned GMRES and Euler/BDF2/Crank-Nicolson time integration.
- Adds a transported Spalart-Allmaras model with the standard `nuTilda` working variable, `fv1/fv2`, clipped modified strain, wall destruction and the nonlinear gradient source.
- Adds standard k-epsilon transport with turbulent production, dissipation, configurable production limiting and bounded eddy viscosity.
- Adds k-omega SST transport with F1/F2 blending, blended diffusion/destruction coefficients, production limiting and cross diffusion.
- Integrates an SST-DES hybrid dissipation switch using cell/grid scale, `C_DES`, and optional F1/F2 zonal shielding.
- Adds a six-component Reynolds-stress transport framework plus epsilon, with LRR-style slow/rapid pressure-strain redistribution, exact tensor production from the supplied velocity gradient and a realizability projection.
- Closes all five remaining Phase-3 turbulence transport/framework tracker items.

## Shared turbulence transport

`include/cfd/solvers/fvm/rans_transport.hpp` is the public entry point. All transported turbulence variables use the same cell-centered polyhedral finite-volume equation:

`d(phi)/dt + div(U phi) - div(Gamma grad(phi)) = source - sink*phi`.

The spatial matrix uses implicit first-order upwind convection and orthogonal face diffusion with a harmonic face diffusivity. Destruction terms are linearized into a non-negative diagonal sink when the model naturally permits it. The same equation supports Euler, BDF2 and off-centered Crank-Nicolson through the existing FVM temporal primitives. Model-specific lower bounds are applied only after a converged linear solve.

## Spalart-Allmaras

`SpalartAllmarasTransport` evolves `nu_tilde`. The model includes molecular-plus-working-variable diffusion, the `Cb2 |grad(nu_tilde)|^2 / sigma` nonlinear source, `fv1/fv2`, modified-strain clipping and wall destruction through `fw`. The public eddy-viscosity output reuses the existing clean-room constitutive helper.

The implementation intentionally follows the no-trip baseline: transition trip terms are not introduced in this release.

## k-epsilon

`KEpsilonTransport` advances `k` and `epsilon` from one shared old-time state. Production is formed from kinematic eddy viscosity and the supplied strain-rate magnitude. The `epsilon` destruction and the `k` dissipation are treated semi-implicitly. Diffusion uses the configurable `sigma_k` and `sigma_epsilon` coefficients.

The model exposes both kinematic and dynamic eddy viscosity and caps the viscosity ratio with a configurable robustness limit.

## k-omega SST and DES

`KOmegaSSTTransport` computes F1 and F2 from wall distance, local `k`, `omega`, molecular viscosity and the `grad(k) . grad(omega)` cross-diffusion indicator. Model coefficients are blended cell-by-cell. The SST viscosity limiter reuses the existing constitutive helper, while the `omega` equation carries signed cross diffusion using source/sink splitting for stability.

When `des_enabled` is set, the k-equation destruction rate is multiplied by an SST-DES factor based on the modeled turbulent length scale and `C_DES * Delta`. `Delta` defaults to the cube root of cell volume but can be supplied explicitly. `none`, `f1` and `f2` zonal-filter modes are available. DDES/IDDES-specific shielding is deliberately left as a later enhancement rather than being mislabeled as this baseline DES formulation.

## Reynolds-stress transport framework

`ReynoldsStressTransport` evolves `Rxx`, `Ryy`, `Rzz`, `Rxy`, `Rxz`, `Ryz` and `epsilon`. It computes the exact Reynolds-stress production tensor

`P_ij = -(R_ik dU_j/dx_k + R_jk dU_i/dx_k)`

from either a velocity field or a directly supplied velocity-gradient tensor. The default closure uses LRR-style slow return-to-isotropy and rapid pressure-strain redistribution. An additional source callback provides a clean extension seam for buoyancy, rotation, wall-reflection or application-specific closures.

After each segregated transport step the symmetric stress tensor is projected back into the realizable covariance set: normal stresses are bounded non-negative, pairwise Cauchy-Schwarz limits are enforced, and off-diagonal stresses are uniformly reduced if needed to restore a non-negative determinant.

## Validation

Focused target:

- `cfd-v0152-rans-transport-tests`

The regression covers:

- exact strain-rate recovery for a linear simple-shear velocity field;
- SA production and positive `nu_tilde`;
- k-epsilon decay without production and growth with homogeneous shear;
- SST F1/F2 bounds and transported production;
- RANS-vs-DES k-dissipation switching on a fine grid;
- Reynolds-stress shear anisotropy and realizability preservation;
- fixed-value turbulent-scalar diffusion bounded by patch extrema.

The complete default CPU regression matrix passes **122/122** tests.

See `UPSTREAM_DOCUMENTATION_REVIEW_0_15_2.md` for the clean-room documentation mapping and model-boundary notes.
