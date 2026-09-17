# cfd_solvers v0.16.1 - Phase 13 contact and collision maturity

v0.16.1 deepens the independent multibody/DEM family introduced in v0.16.0. Project Chrono remains a clean-room capability and validation reference; the implementation below is native cfd_solvers code built on the compact rigid-body/Jacobian architecture already present in the project.

## Added

### Motors and actuators
- bounded linear velocity motor represented as a projected Jacobian row;
- bounded angular velocity motor;
- force/torque limits converted to per-step impulse limits, so motors share the existing PGS solver and obey deterministic saturation.

### Convex contact
- body-local convex hull primitive;
- world-space support mapping;
- GJK intersection search;
- EPA penetration depth, normal and witness-point recovery;
- iteration diagnostics and separated/rotated-convex regression cases.

### Triangle mesh and persistent contact state
- sphere-to-triangle closest-point narrow phase;
- static world meshes or meshes attached to rigid bodies;
- stable triangle feature namespaces;
- persistent contact manifold cache keyed by body pair and feature;
- tangent-history projection when the contact normal changes;
- stale contact removal.

The current mesh path is deliberately a sphere-triangle baseline. General convex-vs-triangle/mesh BVH acceleration is future work even though the Phase-13 baseline mesh/contact-manifold capability is now present.

### History-dependent Mindlin DEM
- persistent tangential displacement integrated over contact lifetime;
- tangent-plane reprojection;
- tangential spring and damping;
- Coulomb return limiting;
- rolling-resistance and cohesion compatibility;
- ExplicitDemSystem uses persistent tangential history by default, with a compatibility switch for the earlier memoryless model.

## Validation

Focused v0.16.1 regressions cover:
- linear and angular motor target velocities;
- motor force saturation;
- overlapping and separated boxes through GJK/EPA;
- rotated convex overlap;
- sphere-triangle penetration and contact normal;
- persistent manifold ageing/removal;
- Mindlin tangential-history accumulation and restoring force.

The previous v0.16.0 regression remains passing, preserving rigid-body, constraint, sphere collision, impulse/penalty contact, DEM and CFD-drag behavior. The complete default project matrix passes **133/133 CTest targets**; the Python ABI smoke reports **0.16.1**. The modified multibody source and v0.16.1 regression also pass a standalone C++20 `-Wall -Wextra -Wpedantic -Werror` syntax check.

## Tracker

Closed Phase-13 items:
- motor/actuator constraint family;
- convex GJK/EPA narrow phase;
- mesh/triangle collision and persistent contact manifolds;
- history-dependent Mindlin tangential spring.

Phase 13 advances from **24/37 (64.9%)** to **28/37 (75.7%)**. Overall progress advances from **637/756 (84.3%)** to **641/756 (84.8%)**.

## Remaining high-value Phase-13 work

- implicit Newmark and generalized-alpha/HHT integration;
- articulated reduced-coordinate dynamics;
- bonded particles and fracture/bond failure;
- GPU/SYCL neighbor search and contact kernels;
- distributed-memory DEM;
- resolved and unresolved CFD/rigid/DEM coupling;
- flexible FEM-to-multibody coupling.
