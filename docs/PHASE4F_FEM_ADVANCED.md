# Phase 4F - Nonlinear, adaptive and additional FEM physics

This checkpoint extends the reusable Elmer-class FEM branch beyond linear elliptic and small-strain proof solvers. The implementation remains clean-room C++20 and reuses the common CSR/Krylov infrastructure rather than introducing per-physics solver stacks.

## Shared nonlinear solve

`cfd::core::newton_solve` provides a small generic Newton driver with:

- residual and CSR Jacobian callbacks;
- ILU(0)-preconditioned restarted GMRES for each Newton linearization;
- residual-based convergence controls;
- backtracking line search so a full Newton step is accepted only when the nonlinear residual decreases.

The first validation problem is

`-div(k grad u) + beta u^3 = f`

on a Tri3 mesh with a manufactured sine solution. The nonlinear FEM path converges in a few Newton iterations and reuses the same sparse infrastructure as the collocated FVM solver.

## Adaptive Tri3 mesh refinement

The adaptive baseline contains three explicit stages:

1. a residual/normal-flux-jump estimator for P1 Poisson problems;
2. Dorfler bulk marking based on squared estimator contribution;
3. conforming longest-edge refinement.

Selected edges are globally shared, and triangles touched by one, two or three split edges use conforming subdivision templates. Boundary edges are split while retaining their patch identifiers. The regression requires that a marked subset produces fewer elements than full red refinement and reduces manufactured-solution error.

This is a baseline h-adaptivity implementation. It does not yet include coarsening, anisotropic metric refinement, high-order error estimators, distributed refinement or repartitioning.

## Darcy porous flow

`Darcy2D` solves saturated isotropic pressure flow using

`v = -(K/mu) (grad(p) - rho g)`.

Pressure uses the shared scalar Tri3 diffusion operator, while element velocity is reconstructed from the pressure gradient. A pressure-driven rectangular channel reproduces the analytical constant Darcy velocity and near-zero transverse flow.

Brinkman/Stokes-Darcy coupling remains a separate unchecked tracker item.

## Magnetostatics

`Magnetostatics2D` solves the scalar out-of-plane vector-potential form

`-div(nu grad(A_z)) = J_z`,

then reconstructs

`B = (dA_z/dy, -dA_z/dx)`.

The manufactured sine-potential regression validates the scalar magnetostatic baseline. Frequency-domain eddy currents and 3-D H(curl)/edge-element Maxwell remain future work.

## Modal/eigenfrequency baseline

The common core now includes a sparse generalized symmetric eigensolver for

`K phi = lambda M phi`.

It uses inverse iteration with M-orthogonal deflation for the lowest modes and shared PCG inner solves. `BarModal1D` assembles Line2 axial stiffness plus consistent mass and validates the first fixed-free longitudinal modes against the analytical frequencies.

The current eigensolver is intended as a readable baseline. Lanczos/LOBPCG, spectral transforms, distributed eigenproblems and large-mode-count workflows remain future optimizations.
