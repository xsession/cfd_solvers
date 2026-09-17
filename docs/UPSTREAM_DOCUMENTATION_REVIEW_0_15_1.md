# Upstream documentation review - v0.15.1 thermal coupling and FVM restart I/O

This is a clean-room documentation review. The implementation in `cfd_solvers` is independent; upstream source code was not copied or mechanically translated.

## OpenFOAM thermal and energy architecture

Reviewed public documentation:

- OpenFOAM 2.3.0 thermal modelling: https://openfoam.org/release/2-3-0/thermal/
- OpenFOAM 2.2.0 thermophysical/multiphase energy: https://openfoam.org/release/2-2-0/thermophysical-multiphase-energy/
- OpenFOAM 5.0 heat/thermophysical modelling: https://openfoam.org/release/5-0/
- OpenFOAM 9 model/source framework: https://openfoam.org/release/9/
- OpenFOAM 2.0 chemistry/pyrolysis: https://openfoam.org/release/2-0-0/chemistry/

Transferable design lessons:

1. **Solid conduction is a first-class region equation.** Solid temperature evolution should use the same finite-volume mesh/operator conventions as fluid energy where possible, but with material density, heat capacity and conductivity explicit in its contract.
2. **Energy sources should be composable.** OpenFOAM's model/source architecture separates physical source models from the core energy equation. The new `ThermalTransport` source callback follows that boundary, allowing chemistry and radiation contributions to be summed without hard-coding a particular chemistry mechanism into the heat solver.
3. **Thermophysical properties may depend on temperature.** The transport layer therefore accepts heat-capacity and conductivity callbacks evaluated against the evolving temperature field.
4. **Radiation belongs in the energy balance, not in a post-processing-only path.** The new reactive/radiative source bridge evaluates the existing `RadiationSourceModel` at each thermal step and combines it with chemistry heat release.
5. **Restart data needs equation state plus mesh identity.** A useful checkpoint must carry the cell/face topology, patch names, physical time/step and field arrays together; a scalar-only dump is not sufficient for a reliable restart boundary.

## v0.15.1 implementation mapping

- `SolidHeatConduction` wraps the production implicit `ScalarTransport` equation with `rho*cp`/conductivity scaling, fixed-temperature/adiabatic boundaries and volumetric heat sources.
- `ThermalTransport::set_reactive_radiative_source(...)` composes chemistry heat release and a lifetime-safe shared radiation model into one volumetric energy source.
- `ThermalTransport::set_material(...)` provides temperature-dependent `cp` and conductivity callbacks with explicit validation.
- `fvm_checkpoint.hpp` provides a stable HDF5 checkpoint boundary for a `PolyMesh`, cell scalar fields, physical time and step index.
- HDF5 files store cell geometry, face geometry/connectivity, boundary-patch names, field names and field arrays. `read_fvm_hdf5_checkpoint(...)` reconstructs a validated `PolyMesh` rather than handing callers raw HDF5 handles.
- Builds without HDF5 fail explicitly; no alternate binary format is silently substituted.

## Validation boundary

The focused thermal regression checks source-energy conservation, bounded fixed-wall conduction, chemistry/radiation source superposition and Boussinesq sign behavior. The checkpoint regression is run twice: once with HDF5 disabled to verify the fail-closed contract, and once with the system HDF5 library enabled to verify an actual mesh/metadata/two-field write-read round trip.

Remaining Phase-3 work after this checkpoint is concentrated in high-order compressible reconstruction, full RANS/DES transport, multiphase/dispersed-phase solvers and topology-changing AMR.
