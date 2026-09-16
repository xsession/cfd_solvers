# Phase 4D / 7A Results: Shared Numerics and Electrochemistry Foundation

Version: 0.5.0

This checkpoint intentionally advances two cross-cutting tracks together:

1. the sparse/numerical infrastructure required by OpenFOAM-class FVM and Elmer-class FEM;
2. the first executable chemistry/electrochemistry/corrosion layer.

## Shared sparse linear algebra

Added:

- CSR matrix storage with validated sorted rows;
- duplicate-merging sparse assembly builder;
- parallel CSR matrix-vector products;
- Jacobi preconditioning;
- preconditioned conjugate gradient;
- BiCGStab;
- ILU(0) factorization/preconditioner;
- restarted GMRES.

The collocated incompressible FVM solver now uses ILU(0)-GMRES for the linearized momentum equation by default. The former fixed-point/Jacobi-like momentum sweep remains available as an A/B reference through `use_krylov_momentum=false`.

The pressure-correction split is preserved correctly: the Krylov path solves the full momentum equation including pressure-gradient forcing, then reconstructs the diagonal-split `HbyA` field used by Rhie-Chow/SIMPLE/PISO/PIMPLE pressure correction.

## Bounded higher-order finite-volume reconstruction

`FaceInterpolationScheme::bounded_linear` adds a MUSCL-style upwind reconstruction with a Barth-Jespersen-style local limiter on arbitrary `PolyMesh` geometry.

Regression gates verify:

- exact reconstruction of a linear scalar field on the Cartesian test mesh;
- no under/overshoot for a discontinuous bounded scalar field.

## Chemistry kinetics foundation

The clean C++ core now includes:

- species name, molar mass and charge metadata;
- Arrhenius rate constants;
- elementary mass-action reactions;
- stoichiometric source-term assembly.

This is intentionally much smaller than Cantera/Reaktoro/PHREEQC. Reversible reactions, thermodynamic phases, activity models, stiff integration and equilibrium remain tracked work.

## Electrochemistry foundation

Added reusable relations for:

- Nernst equilibrium potential;
- Butler-Volmer current density;
- Faradaic molar flux;
- corrosion penetration rate.

### Nonlinear 1-D corrosion cell

`solve_corrosion_cell_1d` closes the nonlinear feedback between electrolyte ohmic loss and interfacial Butler-Volmer kinetics:

- electrolyte surface potential changes with current;
- overpotential changes with surface potential;
- Butler-Volmer current is solved self-consistently;
- Faraday's law converts the anodic current to dissolution flux and recession rate.

The small-overpotential regression is compared against the analytical linearized Butler-Volmer + electrolyte-resistance result.

### Conservative Nernst-Planck transport

The initial 1-D finite-volume transport solver supports:

- multiple species;
- diffusion;
- advection;
- electromigration;
- periodic and no-flux boundaries;
- externally supplied electric potential;
- conservative species update and positivity diagnostics.

Validation includes Fourier-mode diffusion decay, exact periodic amount conservation, and charged-species electromigration conservation.

The same transport model now also runs directly on arbitrary `PolyMesh` geometry. The PolyMesh path adds:

- self-consistent electroneutral current-continuity potential;
- Poisson electrostatics for space-charge / Poisson-Nernst-Planck models;
- fixed-potential and insulating electrostatic boundaries;
- Scharfetter-Gummel exponential fitting for drift-dominated ionic transport;
- Butler-Volmer/Faradaic species-flux electrode boundaries;
- operator-split Poisson-Nernst-Planck stepping.

The electroneutral solver is validated against the analytical binary liquid-junction potential. The Poisson solver is validated against a uniformly charged 1-D slab with fixed-potential endpoints. Scharfetter-Gummel is checked against the direct exponential-fit face flux.

### Galvanic mixed-potential solver

A multi-reaction corrosion-potential solve balances arbitrary Butler-Volmer branches by electrode area. The symmetric two-reaction regression has an exact midpoint mixed potential and zero open-circuit net current.

## Tracker/provenance automation

`docs/INTEGRATION_TRACKER.md` is the authoritative checkbox list. `docs/upstreams.json` pins the upstream revisions used for the current feature research. Run:

```bash
./scripts/integration_status.py
```

to print completion percentages by phase and upstream family.

At the 0.5.0 checkpoint the tracker reports 92/387 completed items (23.8%). This percentage is intentionally strict: planned features are not counted as complete merely because a related proof-of-concept exists.

## Validation commands

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFD_ENABLE_OPENMP=ON -DCFD_ENABLE_NATIVE_ARCH=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/cfd-solve electrochem-corrosion1d
./build/cfd-solve electrochem-pnp1d
./build/cfd-solve electrochem-pnp-poly
./build/cfd-solve electrochem-galvanic
./scripts/integration_status.py
```

The release validation matrix additionally covers Clang, serial GCC and ASan+UBSan.
