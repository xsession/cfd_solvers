# Release 0.15.1 - solid heat, reactive-radiative thermal coupling and FVM HDF5 restart

## Highlights

- Adds a production-facing finite-volume solid heat-conduction wrapper with density, heat capacity, conductivity, fixed-temperature/adiabatic boundaries and volumetric heat generation.
- Adds temperature-dependent thermal material properties and a chemistry-plus-radiation energy-source bridge.
- Makes radiation-source ownership lifetime-safe by retaining a shared model handle inside the thermal source callback.
- Adds native HDF5 finite-volume checkpoints containing mesh topology/geometry, patch names, physical time/step and multiple cell scalar fields.
- Adds HDF5 checkpoint readback that reconstructs and revalidates `PolyMesh` for restart workflows.
- Preserves explicit fail-closed behavior when HDF5 is not compiled in.
- Closes three Phase-3 tracker items: solid heat conduction, radiation/chemistry energy coupling, and HDF5 checkpoint/field output.

## Solid heat conduction

`cfd::fvm::SolidHeatConduction` reuses the implicit scalar finite-volume stack instead of introducing a separate discretization. Thermal diffusivity is derived from `k/(rho*cp)`, while volumetric heat release is converted to the temperature-equation source by `1/(rho*cp)`. This keeps physical material units at the public API while retaining the tested scalar CSR/GMRES path internally.

The focused source-only regression verifies that integrated sensible-energy growth equals integrated volumetric power. A fixed hot/cold wall case verifies boundedness and the expected monotone thermal gradient.

## Reactive-radiative thermal transport

`ThermalTransport` now supports an owned `RadiationSourceModel` plus an optional chemistry heat-release callback. Both are evaluated against the current cell temperature and summed in the same volumetric source contract. Material `cp(T)` and `k(T)` callbacks remain separate so constitutive behavior is not mixed with source physics.

## FVM HDF5 checkpoint format

The checkpoint writer records:

- `/Mesh/Cells`: cell-center XYZ and volume;
- `/Mesh/FaceGeometry`: face-center XYZ and oriented area-vector XYZ;
- `/Mesh/FaceIndex`: owner, neighbour and patch indices;
- `/Mesh/PatchNames`: fixed-width UTF-8-compatible byte table for patch names;
- `/Meta/Time` and `/Meta/Step`;
- `/Meta/FieldNames`;
- `/Fields/<name>` for each cell scalar field.

Field names are restricted to `[A-Za-z0-9_.-]+`, duplicate names and non-finite field values are rejected, and every field must match the mesh cell count. The reader validates dataset shapes and rebuilds `PolyMesh`, which re-applies its owner/neighbour/patch consistency checks.

## Validation

Focused targets:

- `cfd-v0151-thermal-coupling-tests`
- `cfd-v0151-fvm-hdf5-checkpoint-tests`

The complete CPU regression matrix passes **121/121** tests. The HDF5 checkpoint test passes both with HDF5 disabled (expected explicit rejection) and with system HDF5 enabled (real round trip of a sheared mesh, six patches, time/step metadata, temperature and pressure fields).

See `UPSTREAM_DOCUMENTATION_REVIEW_0_15_1.md` for the clean-room documentation mapping.
