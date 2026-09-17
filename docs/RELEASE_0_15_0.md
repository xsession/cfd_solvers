# Release 0.15.0 - production FVM second-order time integration and pressure AMG seam

## Highlights

- Integrates Euler, BDF2/backward and off-centered Crank-Nicolson directly into the production scalar FVM equation solver.
- Integrates BDF2 and Crank-Nicolson into physical PISO/PIMPLE momentum prediction while preserving SIMPLE pseudo-time behavior.
- Adds pressure PCG with built-in Jacobi or external AMG/multigrid V-cycle preconditioning.
- Adds a focused transient-accuracy and pressure-preconditioner regression.
- Closes three Phase-3 integration-tracker items.

## Time integration

`TemporalScheme` is shared by scalar transport and collocated incompressible flow. BDF2 uses an Euler startup step before the old-old state exists. Crank-Nicolson uses the OpenFOAM-style off-centering control, with an implicit spatial weight of `1/(1+off_centering)`.

The scalar solver now stores a reusable spatial CSR operator and synthesizes the transient solve matrix for each scheme. The collocated momentum predictor uses the same equation-level split for PISO/PIMPLE. Time-level history advances once per physical timestep, not once per PIMPLE outer corrector.

## Pressure multigrid path

`PressurePreconditionerKind` selects:

- `none`: legacy matrix-free CG;
- `jacobi`: matrix-free PCG with a pressure diagonal preconditioner;
- `external_amg`: matrix-free PCG invoking the V-cycle callback installed by `set_pressure_amg_cycle(...)`.

The external seam intentionally accepts only residual and correction spans, so hierarchy storage and third-party solver objects do not leak into the FVM API.

## Validation

`cfd-v0150-fvm-temporal-pressure-tests` checks:

- transient scalar diffusion against a fine-time reference;
- BDF2 and centered Crank-Nicolson accuracy relative to backward Euler;
- BDF2 PISO and Crank-Nicolson PIMPLE history advancement;
- incompressibility after the second-order momentum solve;
- real invocation of an external pressure AMG/preconditioner callback.

On the focused scalar case with `dt=0.02`, RMS errors against the fine-time reference were approximately `1.93e-4` for Euler, `4.76e-5` for BDF2 and `3.20e-7` for centered Crank-Nicolson.
