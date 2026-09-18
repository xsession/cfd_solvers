# cfd_solvers v0.16.8 - implicit, articulated and flexible multibody dynamics

v0.16.8 closes four of the five remaining Phase-13 algorithmic gaps. The goal is not to reproduce Project Chrono internals, but to provide readable independent baselines for implicit second-order time integration, reduced-coordinate articulated dynamics and conservative FEM/multibody interface coupling.

## Implicit generalized-coordinate dynamics

`ImplicitGeneralizedIntegrator` advances a dense generalized second-order system

`M a = f(q, v, t)`

with nonlinear Newton iterations. The acceleration at `n+1` is the Newton unknown and Newmark kinematics reconstruct displacement and velocity inside each residual evaluation.

Provided parameterizations:

- average-acceleration Newmark (`beta=1/4`, `gamma=1/2` by default);
- generalized-alpha from a requested high-frequency spectral radius `rho_inf`;
- HHT-alpha for `alpha in [-1/3,0]`.

The generalized-alpha implementation evaluates inertia at `n+1-alpha_m` and forces at `n+1-alpha_f`; Newmark is the `alpha_m=alpha_f=0` special case. The current baseline uses finite-difference Newton tangents and a dense pivoted linear solve, prioritizing correctness/readability for modest generalized-coordinate systems.

## Articulated reduced coordinates

`ArticulatedSystem` provides a fixed-base tree of revolute and prismatic links. Each link stores one generalized coordinate plus a rigid-body mass/inertia model. Forward kinematics reconstruct body poses and Jacobians reconstruct body velocities.

The solver assembles

`M(q) = sum_b (m_b Jv_b^T Jv_b + Jw_b^T I_b Jw_b)`

and projects gravity, body wrench, joint effort and damping into generalized forces before solving for `qdd`.

This is a dense O(n^3) reduced-coordinate baseline. It is intentionally not described as a Featherstone articulated-body algorithm: Coriolis/centrifugal terms and linear-time ABA remain worthwhile future performance/accuracy extensions.

## FEM flexible-body <-> multibody interface

`FlexibleBodyInterface` is a floating-frame/modal coupling seam for FEM interface nodes:

- rigid-frame translation/rotation plus modal deformation maps to world node positions;
- frame twist plus modal rates maps to world node velocities;
- FEM nodal reactions project to a net rigid-body force and torque;
- the same reactions project onto modal generalized forces.

The projection is work-conjugate. The regression explicitly checks that nodal virtual work equals reduced body/modal virtual work, which prevents a coupling adapter from silently creating or destroying mechanical work.

## Validation

Focused v0.16.8 coverage includes:

- average-acceleration Newmark energy conservation for an undamped oscillator;
- HHT high-frequency algorithmic damping;
- generalized-alpha parameter construction;
- exact single-prismatic-link generalized acceleration;
- symmetric coupled two-link articulated mass matrix;
- gravity projection into revolute generalized coordinates;
- FEM interface displacement/velocity transfer;
- resultant force and torque conservation;
- modal-force projection;
- FEM/multibody virtual-work conservation.

The complete default project matrix passes **140/140 CTest targets** and Python ABI **0.16.8**. The full result is recorded in the packaged `ctest_v0168_full.log`.

## Tracker impact

Phase 13 closes:

- implicit Newmark multibody integrator;
- generalized-alpha/HHT multibody integrator;
- articulated reduced-coordinate solver;
- FEM flexible-body <-> multibody coupling.

Phase 13 therefore advances from 32/37 to **36/37 = 97.3%** and overall machine-counted progress from 645/756 to **649/756 = 85.8%**. The only intentionally open Phase-13 capability is distributed-memory DEM qualification on a real multi-rank GPU-aware MPI system.
