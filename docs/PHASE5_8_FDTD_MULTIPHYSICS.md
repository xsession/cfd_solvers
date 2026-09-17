# v0.7.0 - FDTD dispersion and coupled-multiphysics foundation

## Scope

This checkpoint advances two previously under-covered families without changing their numerical identities:

- Phase 5 extends the openEMS-class FDTD baseline with dispersive media, PMC, wave-port/S-parameter handling and VTK output.
- Phase 8 introduces explicit data exchange and partitioned coupling infrastructure, then validates it on a one-way DC-conduction -> Joule-heating -> transient-thermal FEM chain.

## Dispersive 1-D FDTD

`Maxwell1D` now supports three independent auxiliary differential equation (ADE) material models:

- Debye first-order polarization relaxation;
- Drude free-carrier polarization;
- Lorentz damped resonant polarization.

The nondispersive dielectric permittivity stored by the Yee update is interpreted as the high-frequency/infinite-frequency permittivity for dispersive cells. Polarization state is reset when a region is replaced by an ordinary dielectric with `set_material()`.

The current validation is a stability/response baseline, not a complete frequency-domain material-fit certification. Each model is excited by a finite pulse, remains finite under its regression parameters, and develops nonzero polarization. Production use still needs broadband analytical reflection/transmission validation and multi-pole support.

## PMC, ports and VTK

The 3-D Cartesian Maxwell baseline can apply either PEC or PMC enclosure semantics. The PMC regression verifies the magnetic boundary state and finite energy.

`WavePort1D` provides a TEM-like decomposition of E/H samples into forward/backward waves and integrates them using the existing DFT monitor. A synthetic known reflection/transmission signal recovers the expected complex S11/S21 values.

`write_maxwell3d_vtk_ascii()` writes the common logical Maxwell3D lattice as legacy ASCII VTK `STRUCTURED_POINTS` with E and H vectors. It is intended as a portable debugging/interoperability baseline; HDF5/XDMF and stagger-aware production output remain future work.

## Multiphysics field registry and transfers

`FieldRegistry` stores named fields with:

- component/entity counts;
- location (cell/face/node/edge/grid point/surface/global);
- topology family;
- SI base-dimension exponents and scale;
- producer metadata.

The first conservative transfer operators are deliberately small and testable:

- exact interval-overlap remapping of 1-D cell averages, preserving the integrated scalar exactly;
- oriented integrated face flux -> cell volumetric rate on `PolyMesh`, with exact internal cancellation in the volume-weighted global sum.

General unstructured cell-to-cell and FEM<->FVM projection remain unchecked in the tracker.

## Partitioned coupling

`solve_partitioned_fixed_point()` provides reusable fixed-point coupling with optional Aitken relaxation and explicit convergence controls. The regression solves the cosine fixed point and converges in six iterations in the development build.

## Electro-thermal coupling

`JouleHeatingCoupler2D` couples existing FEM modules on a shared Tri3 mesh:

1. solve DC conduction;
2. reconstruct element current density;
3. compute element volumetric Joule heat `q = |J|^2 / sigma = sigma |E|^2`;
4. pass that element-wise source into transient `Heat2D`;
5. advance the thermal state.

The manufactured case uses a 2 m x 1 m rectangle, 4 V potential difference and conductivity 5 S/m. The exact uniform field magnitude is 2 V/m and the exact Joule source is 20 W/m^3. The regression recovers that source to numerical tolerance and verifies a finite interior temperature increase with fixed-temperature boundaries.

## Limitations kept explicit

- No monolithic block coupler yet.
- No FEM<->FVM or structured<->unstructured geometric projection yet.
- No temperature-dependent electrical conductivity feedback yet; the current electro-thermal chain is one-way.
- No CPML/TFSF/NF2FF/SAR or FDTD GPU production kernel yet.
- The 1-D port is a TEM-like validation baseline, not a general waveguide mode solver.
