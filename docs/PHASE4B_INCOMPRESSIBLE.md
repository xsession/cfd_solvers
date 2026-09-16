# Phase 4B - transient incompressible FVM baseline

Phase 4B starts building time-dependent incompressible flow on top of the Phase-4A mesh, operator and linear-solver layer. This checkpoint is a baseline, not OpenFOAM feature parity.

## Reusable face schemes on `PolyMesh`

`cfd/fvm/schemes.hpp` adds data-oriented face operations that use the existing owner/neighbour sign convention:

- distance-weighted linear scalar/vector interpolation;
- first-order upwind scalar/vector interpolation using oriented face flux;
- face volumetric flux from an interpolated velocity field;
- conservative scalar and vector convective divergence.

A global face flux is positive in the stored `Face::area` direction, i.e. owner to neighbour for an internal face. Cell loops flip that sign automatically for the neighbour cell. This keeps conservation explicit and avoids atomics.

The regression suite verifies that upwind transport preserves a constant scalar, linear interpolation gives the exact divergence of a linearly varying scalar in uniform flow on the Cartesian-hexa mesh, positive flux selects the owner value, and uniform vector momentum has zero convective divergence.

## Periodic staggered Navier-Stokes baseline

`cfd::fvm::Incompressible2D` is a periodic MAC/staggered finite-volume baseline with:

- x-velocity on x-normal faces and y-velocity on y-normal faces;
- conservative flux-form momentum convection;
- centered viscous diffusion;
- explicit momentum prediction;
- matrix-free CG pressure Poisson solve;
- conservative pressure correction of face velocities;
- periodic pressure nullspace removal;
- Taylor-Green analytical velocity-error measurement.

The staggered baseline deliberately reuses the same discrete face-flux/divergence principle as the Phase-4A pressure projection. It is useful as a correctness oracle for the later collocated `PolyMesh` SIMPLE/PISO implementation.

## Taylor-Green regression

The normal CLI case uses a 64x64 periodic domain, `nu=0.01`, `dt=1e-4`, amplitude `0.05`, and 100 steps. In the development container it reports approximately:

- simulated time: `0.01`;
- divergence RMS: `8.3e-15`;
- velocity RMS error against analytical viscous Taylor-Green decay: `1.5e-7`;
- final pressure correction: converged.

A separate regression advances 24x24 and 48x48 cases to the same physical time with the fine time step reduced by four. The fine-grid velocity error must be less than 45% of the coarse-grid error. In the current development run the measured errors were `5.56e-7` and `1.39e-7`, a ratio of about `0.250`. This is an accuracy/refinement gate, not a claim about the asymptotic order on arbitrary meshes.

## Phase-4C continuation

The collocated `PolyMesh` momentum path, Rhie-Chow-style flux, non-orthogonal pressure correction, SIMPLE/PISO/PIMPLE-style loops and wall/pressure boundary objects are now implemented in Phase 4C. See `PHASE4C_COLLOCATED.md`. Remaining work includes bounded higher-order convection/gradient limiters, a general sparse fvMatrix layer, advanced cyclic/coupled patches and turbulence models.
