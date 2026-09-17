# Upstream documentation review for v0.15.2 turbulence transport

This review records public model concepts used to guide an independent implementation. No upstream source code was copied. The project keeps its own C++20 data structures, finite-volume assembly, linear solvers and APIs.

## Public references reviewed

### Spalart-Allmaras

OpenFOAM Foundation Source Code Guide, `SpalartAllmaras` model:

- https://cpp.openfoam.org/v13/classFoam_1_1RASModels_1_1SpalartAllmaras.html

Useful public model facts: the model is a one-equation RANS closure; the documented baseline omits the trip term; the modified strain is limited for robustness; and the documented default family includes `Cb1`, `Cb2`, `Cw2`, `Cw3`, `Cv1`, `Cs`, `sigmaNut` and `kappa`.

Clean-room mapping: `SpalartAllmarasTransport` owns the transport equation and obtains eddy viscosity through the pre-existing independent constitutive helper. The FVM operator, temporal integration and linear solve are project-native.

### Standard k-epsilon

OpenFOAM Foundation Source Code Guide, `kEpsilon` model:

- https://cpp.openfoam.org/v12/classFoam_1_1RASModels_1_1kEpsilon.html

The public documentation identifies the standard two-equation k-epsilon model and the familiar coefficient family `Cmu`, `C1`, `C2`, `sigmak` and `sigmaEps`.

Clean-room mapping: the implementation advances k and epsilon as separate bounded finite-volume equations, with eddy-viscosity production derived from strain-rate magnitude and semi-implicit dissipation.

### k-omega SST

OpenFOAM Foundation Source Code Guide, `kOmegaSST` model:

- https://cpp.openfoam.org/v13/classFoam_1_1kOmegaSST.html

The public model boundary is the two-equation SST formulation with blending between near-wall and outer-region coefficient sets, an eddy-viscosity limiter and cross-diffusion behavior.

Clean-room mapping: v0.15.2 computes F1/F2 and cross diffusion independently from the cell fields and gradients. It does not depend on OpenFOAM field classes, matrix classes, run-time selection machinery or wall-function implementation.

### Reynolds-stress models

OpenFOAM Foundation Source Code Guide, LRR model:

- https://cpp.openfoam.org/v13/classFoam_1_1RASModels_1_1LRR.html

The public documentation describes the Launder-Reece-Rodi Reynolds-stress model, generalized gradient diffusion and optional wall reflection.

Clean-room mapping: this release implements the transport framework and an LRR-style slow/rapid pressure-strain baseline. It adds a project-native realizability projection and an explicit source extension seam. Gibson-Launder wall reflection is not claimed in v0.15.2.

### SST-DES

OpenFOAM Foundation Source Code Guide, `kOmegaSSTDES`:

- https://cpp.openfoam.org/v13/classFoam_1_1LESModels_1_1kOmegaSSTDES.html

The public documentation describes SST-DES as a hybrid RANS/LES model, exposes a `C_DES` coefficient, a turbulent length scale, a DES dissipation-rate multiplier and optional zonal filtering using SST blending functions.

Clean-room mapping: v0.15.2 integrates a DES dissipation multiplier directly into the project-native SST k equation. The grid scale is explicit and defaults to `cbrt(cell volume)`. F1/F2 zonal shielding can be selected. DDES and IDDES remain future variants.

## Architectural conclusions carried into cfd_solvers

1. Transported turbulence state must be separate from algebraic constitutive helpers.
2. The generic FVM equation should own advection, diffusion, temporal history, boundary treatment and linear solves; turbulence models should provide coefficients and sources.
3. Positivity/realizability safeguards are numerical model boundaries, not post-processing conveniences.
4. Hybrid RANS/LES must alter the transported model dynamics; a standalone DES length-scale helper is insufficient.
5. Reynolds-stress transport needs an explicit closure/source seam because pressure-strain, wall reflection, buoyancy and rotation corrections evolve independently.

## Intentional scope boundaries

- SA transition/trip terms are not implemented.
- k-epsilon compressibility/RDT extensions are not yet exposed.
- SST wall functions remain in the existing wall-model layer rather than being embedded in these transport classes.
- The Reynolds-stress baseline does not yet include Gibson-Launder wall reflection or SSG pressure-strain.
- SST-DES is implemented; DDES and IDDES are not yet claimed.
