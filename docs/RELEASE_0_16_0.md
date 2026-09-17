# cfd_solvers v0.16.0 - Phase 13 multibody and DEM foundation

v0.16.0 starts the solver-family expansion identified by the uncovered-domain audit: constrained multibody dynamics and granular/DEM mechanics. Project Chrono is used as a capability and validation reference under its permissive BSD-3-Clause license; the cfd_solvers implementation is independent and intentionally small/readable.

## Implemented baseline

- 6-DOF `RigidBodyState` with quaternion orientation, linear/angular velocity and world-space inertia application.
- `RigidBodySystem` with gravity, force/torque accumulation, attached sphere/plane collision geometry and deterministic timestepping.
- Semi-implicit Euler and velocity-Verlet rigid-body integration.
- Jacobian-row constraint representation with projected Gauss-Seidel impulses.
- Distance, spherical, revolute, prismatic, fixed and gear-ratio constraints.
- Sphere AABBs plus sweep-and-prune and BVH broad phases.
- Sphere-sphere and sphere-plane narrow-phase contact.
- Smooth penalty contact with Coulomb-limited tangential damping.
- Non-smooth unilateral impulse/contact-complementarity baseline with restitution and friction.
- Explicit spherical DEM with Hertz-type normal contact, friction, rolling resistance and cohesive-force baseline.
- Stokes-drag/reaction-force primitive for later unresolved CFD/DEM coupling.

## Validation

`tests/test_v0160_multibody_dem.cpp` covers:

- gravitational free-fall velocity and quaternion normalization,
- distance-constraint correction,
- revolute-joint locked/free angular DOFs,
- sweep-and-prune/BVH candidate parity,
- sphere-sphere penetration geometry,
- equal-and-opposite smooth contact forces,
- non-smooth restitution impulse response,
- finite DEM wall support,
- Stokes drag direction.

The new multibody source also compiles independently with `-Wall -Wextra -Wpedantic -Werror -Wconversion -Wshadow`. The complete default/Python-enabled project matrix passes **132/132 CTest targets**.

## Tracker impact

Phase 13 introduces 37 explicit capabilities, with 24 completed in this baseline. The machine-counted project total becomes **637/756 (84.3%)**. The percentage decreases slightly because the roadmap denominator now includes a substantial new solver family; previously implemented capability count is not removed.

Remaining high-value Phase-13 work includes convex/mesh contact, history-dependent Mindlin tangential springs, bonded-particle fracture, implicit Newmark/generalized-alpha/HHT stepping, actuators, flexible-body/FEM coupling, GPU DEM and distributed DEM.
