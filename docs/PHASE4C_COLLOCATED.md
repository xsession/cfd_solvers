# Phase 4C - collocated pressure/velocity coupling baseline

Phase 4C moves the incompressible finite-volume path from the periodic staggered
correctness solver to a cell-centred/collocated solver on `PolyMesh`. The goal of
this checkpoint is a readable pressure-velocity coupling architecture that can
later grow into higher-order OpenFOAM-class algorithms without hiding the
stabilization or mesh corrections inside opaque kernels.

This is a clean-room implementation. OpenFOAM 14 is used as an architectural
reference for the separation between momentum prediction, pressure-free
`HbyA`, inverse momentum diagonal, pressure equation, non-orthogonal correction,
flux correction, and velocity reconstruction. No OpenFOAM source is copied or
mechanically translated.

## New modules

- `cfd/fvm/pressure_velocity.hpp`
  - velocity BCs: fixed value, zero gradient, slip;
  - pressure BCs: fixed value, zero gradient;
  - face orthogonal/non-orthogonal geometric decomposition;
  - Rhie-Chow-style collocated face flux;
  - conservative face-flux divergence diagnostics.
- `cfd/solvers/fvm/collocated_incompressible.hpp`
  - transient/pseudo-transient cell-centred momentum equation;
  - reusable momentum diagonal and pressure mobility `V/aP`;
  - non-orthogonal pressure equation solved by the shared matrix-free CG core;
  - SIMPLE, PISO and PIMPLE-style control loops;
  - conservative corrected face flux retained as the continuity field.
- `make_sheared_cartesian_hexa_mesh()`
  - deterministic affine xy-shear generator for non-orthogonal regression.

## Rhie-Chow face flux

For an internal face, the stored area vector points owner -> neighbour. The
owner-neighbour displacement `d` is used to split the face area into

```
Sf = E + T
E  = d * dot(Sf,d) / dot(d,d)
T  = Sf - E
```

With pressure-free momentum predictor `HbyA` and cell pressure mobility
`D = V/aP`, the face flux is built as

```
phi_f = interpolate(HbyA)_f . Sf
      - D_f * [ (p_N-p_P) * dot(Sf,d)/dot(d,d)
                + interpolate(grad(p))_f . T ]
```

The direct owner-neighbour pressure difference is important: a collocated
checkerboard pressure mode that can disappear from a cell-centred central
pressure gradient still creates a face-flux response. The regression suite
explicitly checks this mode.

## Non-orthogonal pressure correction

The orthogonal term is kept in the matrix-free SPD pressure operator. The
`grad(p).T` contribution is treated explicitly and refreshed through the
configured non-orthogonal correction loop. Fixed-pressure boundary faces add
their diagonal/source contribution; zero-gradient pressure faces contribute no
normal pressure flux.

Pure-Neumann pressure cases remain in the mean-zero compatible subspace. This
is used by the closed lid-driven cavity. Cases with at least one fixed-pressure
patch do not require nullspace removal.

## Momentum equation

The baseline momentum equation uses:

- implicit transient/pseudo-transient diagonal;
- implicit orthogonal viscous coupling;
- first-order upwind convective coupling when convection is enabled;
- explicit cell pressure gradient;
- configurable Jacobi-style momentum sweeps;
- `HbyA` and `V/aP` retained after the final sweep for pressure coupling.

The pseudo-time term is intentionally retained for SIMPLE robustness. It is not
yet the final sparse-matrix/under-relaxed steady momentum implementation.

## Coupling algorithms

### SIMPLE baseline

`iterate_simple()` performs one pseudo-time momentum prediction, one pressure
sequence, pressure relaxation and velocity relaxation. `solve_simple()` repeats
until the requested velocity/continuity gates or the iteration limit.

### PISO baseline

`step_piso()` freezes the physical time-source velocity, performs one momentum
prediction and the configured number of pressure corrections without SIMPLE
relaxation, then advances physical time.

### PIMPLE baseline

`step_pimple()` repeats the momentum + PISO pressure sequence for the configured
outer-corrector count before advancing time. Rebuilding `HbyA` in each outer
loop is what makes the current PIMPLE path materially stronger than simply
repeating the same pressure solve.

These are compact baseline control loops, not full OpenFOAM feature parity.
There is not yet a general fvMatrix class, turbulence coupling, dynamic mesh,
MRF, source framework, or advanced pressure-reference/constraint machinery.

## Validation cases

### Pressure-driven Poiseuille channel

The 32x16 CLI/test case uses fixed kinematic pressure at left/right, no-slip top
and bottom, and slip front/back. Convection is disabled so the exact laminar
solution is the pressure-driven parabola.

Development-container snapshot:

- profile relative L2 error: `~4.7e-3`;
- continuity L2: `~2.2e-8`;
- maximum x velocity: `~1.249e-2`;
- pressure CG converged.

### Lid-driven cavity

A 20x20 collocated SIMPLE case uses a unit moving top wall and stationary other
in-plane walls. The current low-Re regression checks topology/sign rather than
claiming a Ghia benchmark match:

- centre `u_x`: approximately `-0.191`, showing the primary clockwise vortex;
- continuity L2: order `1e-7`;
- pressure CG converged.

### Sheared mesh PISO/PIMPLE

A 12x10x1 mesh with xy shear `0.35` starts from a deliberately divergent field.
Typical development-container values:

- initial continuity L2: `~3.01e-1`;
- after one PISO step: `~1.20e-4`;
- after one two-outer-corrector PIMPLE step from the same initial field:
  `~2.00e-5`.

This case exercises non-orthogonal face decomposition and edge-independent
cell-centred pressure correction rather than only orthogonal Cartesian grids.

## CLI

```bash
./build/cfd-solve fvm-collocated-channel
./build/cfd-solve fvm-collocated-cavity
./build/cfd-solve fvm-collocated-skew
```

## Deliberate next work

- bounded higher-order convection schemes and gradient limiters;
- sparse reusable finite-volume matrix assembly/preconditioners;
- improved transient Rhie-Chow time correction;
- generalized periodic/cyclic/coupled patches;
- skewness-corrected face interpolation beyond the current non-orthogonal
  pressure correction;
- higher-Re cavity/channel benchmarks and mesh-convergence studies;
- turbulence models and wall functions.
