# Deep-search: making the electromagnetic stack more useful

## Decision

The highest-value locally verifiable next increment was a small PEC surface-current
Method of Moments (MoM) path using Rao-Wilton-Glisson (RWG) basis functions. It connects
the existing thin-wire MoM, RF network, FDTD/NF2FF, FEM and geometry-import layers without
pretending that a dense reference assembler is already an electrically-large production
boundary-element solver.

The implementation is clean-room C++20 and uses only the repository's existing dense complex
matrix utilities. It accepts the native triangle surface type used by the OBJ/ASCII-STL readers,
builds both interior and boundary half-RWG
functions, assembles a regularized low-order electric-field integral equation, supports a
plane-wave or delta-gap excitation, reports a relative linear residual and feed impedance, and
evaluates a current-based far field.

## What the reference projects imply

The Bempp documentation identifies RWG as the lowest-order tangential H(div) surface space and
describes the Maxwell electric-field operator as the combination of a vector-potential term and
a surface-divergence term. Its operator interface also makes the dense-versus-FMM assembly
boundary explicit. Sources:

- https://bempp.com/handbook/api/function_spaces.html
- https://bempp.com/handbook/api/boundary_operators.html

SCUFF-EM describes a C++ surface-current boundary-element stack with scattering and RF workflows.
That confirms the practical value of making the current solution a first-class reusable result,
rather than exposing only a one-off antenna example:

- https://github.com/HomerReid/scuff-em/blob/master/README
- https://github.com/HomerReid/scuff-em/blob/master/doc/ImplementationSupplement.md

For the next optics/EM slice, MPB is the useful reference boundary: a photonic-band solver needs
periodic/Bloch wavevectors, a frequency-domain eigenproblem, iterative eigenvalue support and
explicit field output. That is a separate solver family and should not be hidden inside the RWG
baseline:

- https://mpb.readthedocs.io/

## Implemented now

- `include/cfd/rf/surface_mom.hpp` and `src/rf/surface_mom.cpp`
- boundary half-RWG and interior RWG topology construction
- small dense PEC EFIE reference assembly
- plane-wave and delta-gap excitation
- feed current/impedance and normalized residual diagnostics
- electric-current far-field post-processing
- didactic CSV example: `phase12_surface_mom`
- deterministic regression for topology, residual, scaling and far-field finiteness

## Explicit limits

The current increment does not close dielectric surface integral equations, singular self/edge
quadrature, Calderón preconditioning, MLFMM, fast GPU assembly, or NEC-grade convergence. The
regularized three-point rule is useful for small comparison cases and API development; it is not
a certification-quality replacement for a singular-quadrature BEM implementation.

## Recommended next order

1. Add canonical RWG convergence cases (PEC sphere/plate and current continuity/error metrics)
   with singular or Duffy-transformed quadrature.
2. Add dielectric PMCHWT/SIE blocks and surface equivalence coupling to the existing FDTD/FEM
   ports.
3. Add matrix-free near/far interaction separation and an optional FMM backend; keep the dense
   assembler as the deterministic reference.
4. Add physical-optics and shooting-and-bouncing-rays current/ray baselines for electrically large
   bodies, then define a tested MoM/FDTD/asymptotic handoff contract.
5. Implement a separate 2-D periodic Bloch eigenmode solver for the remaining photonic-band rows,
   with HDF5/VTK-independent field export first and optional native HDF5 second.
