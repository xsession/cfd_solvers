# Geant4 clean-room learning notes for cfd_solvers

Geant4 is used here as an architectural and validation-reference source only. No
Geant4 source code is copied or translated. The useful concepts for this project
are the separation between tracks, steps, geometry/material navigation, physics
processes and scoring.

## Source observations

- Geant4 describes itself as a simulation toolkit for the passage of particles
  through matter, with applications in high-energy, nuclear, accelerator,
  medical and space science.
- The public source tree is organized into large subsystems such as geometry,
  materials, particles, processes, tracking, physics lists and examples.
- Its process abstraction separates *proposed interaction length* from the
  subsequent process action. A process exposes `AlongStepDoIt`,
  `PostStepDoIt`, `AtRestDoIt` and GPIL methods that propose step lengths.
- The stepping manager selects a physical step, stores it on the step/track,
  invokes along-step processes, updates the track, then invokes post-step
  processes and sends the completed step to sensitive-detector/scoring logic.
- The process tree separates electromagnetic, hadronic, optical, decay,
  transportation, scoring, cuts and biasing families.

## Clean-room mapping implemented in v0.10.8

`cfd_solvers` now has a compact native transport layer rather than a Geant4 port:

| Geant4 idea | Native v0.10.8 counterpart |
| --- | --- |
| `G4Track` | `TransportTrack` |
| `G4Step` | `TransportStep` |
| Geometry/material navigation | `TransportRegion`, `TransportMaterial`, `locate_region`, `distance_to_box_boundary` |
| GPIL process proposal | `TransportProcess::physical_interaction_length_m` and geometry/user-step comparison |
| Along-step continuous physics | `continuous_energy_loss` stopping-power process |
| Post-step discrete physics | `discrete_interaction` secondary-production process |
| Cuts | `production_cut_energy_ev` and `TransportConfig::energy_cut_ev` |
| Scoring | `TransportScoring` total deposited energy, track length, step/secondary counts |

## Validation exercises

- Continuous loss deposits the expected `dE/dx * step` energy.
- A discrete interaction with shorter GPIL than the user and geometry limits
  creates a secondary above the production cut.
- A geometry boundary limits a track before a user step when the track reaches a
  slab interface.
- The run-level scorer accumulates steps, path length, energy deposition and
  secondaries.

## Future implementation path

1. Add sampled interaction lengths and deterministic RNG streams per track.
2. Add physics-list objects that own ordered process collections by particle.
3. Add hierarchical geometry and BVH acceleration instead of slab-only regions.
4. Add sensitive-detector/hit collections and dose scorers.
5. Add hadronic/optical process families as clean-room model hooks.
6. Couple particle energy deposition into the existing bioheat/FEM/FDTD
   workflows.


## v0.10.9 extension

v0.10.9 keeps the clean-room boundary and extends the local transport layer with stochastic interaction-length sampling, physics-list grouping, BVH lookup, sensitive-detector hit collection, dose-grid scoring and SAR projection to Pennes bioheat. These are intentionally compact reference abstractions; they are not Geant4 source translations and they do not attempt to reproduce Geant4's production physics data tables.
