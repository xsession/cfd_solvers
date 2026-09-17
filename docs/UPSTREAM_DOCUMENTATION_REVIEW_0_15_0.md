# Upstream documentation review - v0.15.0 FVM time integration and pressure multigrid seam

This is a clean-room documentation review. It records numerical/API lessons only; no upstream source was copied or mechanically translated.

## OpenFOAM time discretization

Reviewed documentation:

- OpenCFD time-scheme overview: https://doc.openfoam.com/2606/tools/processing/numerics/schemes/time/
- OpenCFD scalar transport time-scheme example: https://doc.openfoam.com/2606/tools/processing/numerics/schemes/time/example/
- OpenFOAM Foundation v13 `CrankNicolsonDdtScheme`: https://cpp.openfoam.org/v13/classFoam_1_1fv_1_1CrankNicolsonDdtScheme.html

Transferable lessons:

1. A production FVM time scheme belongs in equation assembly, not as a post-processing update. The old and old-old field levels participate directly in the matrix and RHS.
2. Backward/BDF2 needs an additional old time level and is formally second order after startup. `cfd_solvers` therefore uses Euler for the first BDF2 step, then switches to the constant-step BDF2 coefficients.
3. Crank-Nicolson blends the new and old spatial residuals. The Foundation documentation defines its implicit coefficient as `1/(1+ocCoeff)`, where zero off-centering recovers Euler and one is fully centered.
4. Scalar transport is a useful independent regression surface for time schemes before coupling them to pressure/velocity algorithms.

Implemented follow-through:

- `TemporalScheme::{euler,backward_bdf2,crank_nicolson}` is shared by production FVM solvers.
- `ScalarTransport` separates spatial assembly from transient assembly and uses BDF2/CN inside the CSR solve.
- `CollocatedIncompressible` applies the same history and CN spatial-residual blend to PISO/PIMPLE momentum equations, while SIMPLE deliberately retains Euler pseudo-time semantics.

## PETSc multigrid/AMG preconditioning

Reviewed documentation:

- PETSc `PCGAMG`: https://petsc.org/main/manualpages/PC/PCGAMG/
- PETSc `PCMG`: https://petsc.org/main/manualpages/PC/PCMG/
- PETSc KSP linear-solver manual, AMG section: https://petsc.org/main/manual/ksp/

Transferable lessons:

1. Multigrid is best represented as a preconditioner boundary around a Krylov solver, with the hierarchy/V-cycle implementation replaceable independently of the pressure equation.
2. AMG implementations need hierarchy-specific information and tuning; the pressure solver should not hard-code one third-party package or hierarchy representation.
3. A generic `r -> z` V-cycle seam allows native, PETSc, hypre, AmgX or future in-tree multigrid to share the same pressure PCG path.

Implemented follow-through:

- The collocated pressure correction can now use unpreconditioned CG, built-in Jacobi PCG, or an externally supplied AMG/multigrid cycle through `set_pressure_amg_cycle(...)`.
- The pressure operator remains matrix-free and the AMG seam is only a preconditioner callback; this keeps Rhie-Chow/non-orthogonal assembly independent of hierarchy ownership.
- Regression coverage proves the external preconditioner is invoked by the production pressure solve rather than existing only as an isolated adapter utility.

## Deliberate limits

This checkpoint does not claim a new in-tree algebraic coarsener. The production pressure equation now has the required multigrid/AMG integration path, while hierarchy construction remains a separate backend concern. Variable-step BDF2 is also not implemented yet; the new BDF2 path is the validated constant-step form used by the current fixed-`dt` solver controls.
