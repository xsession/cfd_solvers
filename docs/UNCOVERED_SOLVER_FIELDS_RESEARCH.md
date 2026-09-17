I compared the current `cfd_solvers` capability tracker against other major open-source scientific solvers. Your project is already much broader than a CFD collection: chemistry/electrochemistry/corrosion are substantial, SPICE/RF circuitry is already a dedicated phase, and Phase 12 now contains relativistic EM-PIC, plasma chemistry, Geant4-style particle transport, SAR/bioheat, 3-D frequency-domain EM, and cable/EMC work.

That changes what is actually worth adding next.

### Best genuinely uncovered solver fields

| PriorityNew fieldStrong OSS referencesWhat it would add to `cfd_solvers` |                                              |                           |                                                                                                                           |
| ------------------------------------------------------------------------ | -------------------------------------------- | ------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| **1**                                                                    | **Multibody dynamics + DEM**                 | **Project Chrono**        | 6-DOF rigid bodies, joints, constrained DAEs, frictional contact, collision detection, granular/DEM, flexible bodies, FSI |
| **2**                                                                    | **Acoustics / ultrasound / photoacoustics**  | **k-Wave**                | k-space pseudospectral solver, nonlinear acoustics, heterogeneous media, power-law absorption, ultrasound transducers     |
| **3**                                                                    | **Semiconductor TCAD**                       | **DEVSIM**                | Poisson + carrier drift-diffusion, semiconductor contacts, DC/AC/transient device simulation                              |
| **4**                                                                    | **Battery cell physics**                     | **PyBaMM**                | SPM/SPMe/DFN models, thermal electrochemistry, degradation, particle cracking                                             |
| **5**                                                                    | **Nuclear neutronics**                       | **OpenMC**                | continuous-energy neutron/photon Monte Carlo, criticality, tallies, depletion/burnup                                      |
| **6**                                                                    | **Electrical power systems**                 | **pandapower / VeraGrid** | load-flow, OPF, short circuit, contingency, RMS/EMT grid simulation                                                       |
| **7**                                                                    | **Geoscience / THMC / geomechanics**         | **OpenGeoSys + PyLith**   | coupled thermo-hydro-mechanical-chemical porous/fractured media, faults, friction                                         |
| **8**                                                                    | **Crash / impact / blast dynamics**          | **OpenRadioss**           | explicit nonlinear structural dynamics, large deformation, high-speed contact                                             |
| **9**                                                                    | **Molecular dynamics / atomistics**          | **HOOMD-blue + LAMMPS**   | MD/MC, atomistic potentials, thermostats, neighbor lists, GPU particle kernels                                            |
| **10**                                                                   | **Electronic structure / quantum materials** | **Quantum ESPRESSO**      | DFT, plane waves, pseudopotentials, phonons/DFPT, ab-initio MD                                                            |

## 1. Project Chrono — probably the highest-value new family

**Project Chrono** is especially well matched to your architecture. It covers large rigid-body systems, constrained dynamics, flexible/deformable bodies, granular dynamics using smooth and nonsmooth contact methods, fluid-solid interaction, and sensor/robotics simulation. It is permissively licensed. ([GitHub](https://github.com/projectchrono/chrono?utm_source=chatgpt.com "GitHub - projectchrono/chrono: High-performance C++ library for multiphysics and multibody dynamics simulations · GitHub"))

This fills a conspicuous hole between your existing:

- structural FEM,
- particle tracking,
- collision/coalescence infrastructure,
- CFD,
- FSI/multiphysics transfer,
- and optimization framework.

I would create:

**Phase 13 — Multibody, rigid-body and granular dynamics**

Core targets:

```text
RigidBodyState
RigidBodySystem
Constraint/Jacobian framework

RevoluteJoint
PrismaticJoint
SphericalJoint
FixedJoint
DistanceConstraint
GearConstraint

BroadPhaseCollision
BVHCollision
NarrowPhaseContact

PenaltyContact
ComplementarityContact
CoulombFriction

ExplicitDEM
HertzMindlinContact
RollingResistance
CohesiveParticleContact

RigidBodyIntegrator
NewmarkIntegrator
GeneralizedAlphaIntegrator

CFD <-> rigid body
DEM <-> CFD
FEM <-> multibody
particle <-> contact geometry
```

This would also let `cfd_solvers` simulate mechanisms, robotics, suspensions, gears, granular hoppers, powders, rock particles and fluidized beds instead of treating particles mainly as transported point entities.

Chrono is also a better strategic reference than making old LIGGGHTS the main DEM source family; Chrono gives you DEM **and** the surrounding multibody/contact architecture in one system.

---

## 2. k-Wave — a numerical method you currently don't really have

Your EM FDTD and structural mechanics can represent some wave phenomena, but that is not equivalent to a dedicated high-performance acoustic solver.

**k-Wave** implements linear and nonlinear acoustic propagation in 1-D/2-D/3-D heterogeneous media using a **k-space pseudospectral method**. It includes frequency-dependent power-law absorption, nonlinear propagation and PML treatment. ([GitHub](https://github.com/ucl-bug/k-wave?utm_source=chatgpt.com "GitHub - ucl-bug/k-wave: A MATLAB toolbox for the time-domain simulation of acoustic wave fields · GitHub"))

This deserves its own phase:

**Phase 14 — Acoustics and ultrasound**

```text
linear acoustic pressure/velocity solver
heterogeneous rho/c
FFT spectral derivatives
k-space temporal correction
split-field PML

power-law acoustic absorption
fractional Laplacian operator
nonlinear Westervelt terms

point/plane/transducer sources
sensor arrays
time reversal
beamforming
photoacoustic initial-pressure propagation

acoustic intensity
radiation pressure
acoustic heating
```

Then couple it to:

```text
acoustics -> structural vibration
acoustics -> CFD
ultrasound -> Pennes bioheat
acoustics -> piezoelectric FEM
acoustics -> optimization/inverse problem
```

The particularly useful part here is that it introduces a **spectral numerical backend**, rather than simply another PDE implemented with FVM/FEM.

---

## 3. DEVSIM — semiconductor TCAD

This is one of the strongest opportunities because most of the mathematical building blocks already exist in your project.

**DEVSIM** is an Apache-licensed TCAD simulator supporting user-defined PDE systems and 1-D/2-D/3-D semiconductor simulations, including DC, transient and small-signal AC/impedance analyses. ([GitHub](https://github.com/devsim/devsim?utm_source=chatgpt.com "GitHub - devsim/devsim: TCAD Semiconductor Device Simulator · GitHub"))

A new phase could reuse:

- electrostatics,
- nonlinear Newton infrastructure,
- Scharfetter-Gummel fluxes already developed for electrochemistry,
- heat conduction,
- SPICE MNA,
- sparse solvers,
- adaptive meshes.

### Phase 15 — Semiconductor devices / TCAD

Start with:

```text
Poisson equation
electron continuity
hole continuity

electron/hole drift
diffusion
Scharfetter-Gummel discretization

doping profiles
Fermi statistics
mobility models

SRH recombination
Auger recombination
radiative recombination

ohmic contacts
Schottky contacts
insulating boundaries

PN diode
PIN diode
BJT
MOS capacitor
MOSFET
IGBT/power device foundation
```

Then:

```text
DC operating point
transient
small-signal AC
C-V
I-V sweeps
electrothermal self heating
TCAD <-> SPICE compact-model extraction
```

This could become an unusually strong differentiator for the unified project.

---

## 4. PyBaMM — battery-specific models

Your tracker already explicitly leaves **battery/electrolyzer porous electrochemistry + thermal + flow** unfinished.

**PyBaMM** is therefore almost a direct specification source for this missing vertical. It implements battery model families including:

- SPM,
- SPMe,
- Doyle-Fuller-Newman / P2D,
- thermal models,
- degradation,
- particle mechanics/cracking and related extensions. ([PyBaMM Dokumentáció](https://docs.pybamm.org/en/latest/source/api/models/lithium_ion/dfn.html?utm_source=chatgpt.com "Doyle-Fuller-Newman (DFN) — PyBaMM v25.10.3.dev50+g18461a186 Manual"))

I would make this:

### Phase 16A — Batteries

```text
Single Particle Model
SPMe
DFN/P2D

solid diffusion
electrolyte concentration
electrolyte potential
solid potential
Butler-Volmer electrodes

1D through-cell thermal
lumped thermal
3D thermal coupling

SEI growth
lithium plating
active-material loss
particle swelling/cracking

EIS
drive-cycle simulation
cell balancing interface
```

A lot of it can be built almost as a composition layer over Phase 7 rather than another independent PDE framework.

---

## 5. OpenMC — nuclear transport

Your Geant4-inspired work gives you general particle-through-material infrastructure, but **reactor neutronics is a different solver domain**.

OpenMC implements continuous-energy Monte Carlo neutron/photon transport under an MIT license, including nuclear-data-driven interactions, tallies and depletion/transmutation capabilities. ([GitHub](https://github.com/openmc-dev/openmc?utm_source=chatgpt.com "GitHub - openmc-dev/openmc: OpenMC Monte Carlo Code · GitHub"))

### Phase 16B — Neutronics

Important additions:

```text
continuous-energy nuclear cross sections
energy-dependent reaction sampling

neutron elastic scattering
inelastic scattering
capture
fission

secondary neutron sampling
delayed neutrons

constructive solid geometry
material mixtures

fixed-source transport
k-effective eigenvalue solver

tallies
flux spectra
reaction rates
heating

depletion / burnup
Bateman equations
nuclide transmutation
```

The existing generic Monte Carlo transport framework could be reused underneath this.

That makes this much more practical than starting nuclear simulation from zero.

---

## 6. Power-system solver

A complete engineering platform should arguably contain a network-scale electrical solver between SPICE and electromagnetic field simulation.

Two useful references emerged.

**pandapower** is BSD-licensed and provides electrical power-system analysis and optimization around established power-flow algorithms. ([GitHub](https://github.com/e2nIEE/pandapower?utm_source=chatgpt.com "GitHub - e2nIEE/pandapower: Convenient Power System Modelling and Analysis based on PYPOWER and pandas · GitHub"))

**VeraGrid** — formerly GridCal — covers an even broader capability map: AC/DC power flow, unbalanced 3-phase analysis, short circuits, AC/DC OPF, PTDF/LODF, stochastic/continuation power flow, contingency analysis, RMS and EMT studies. ([GitHub](https://github.com/SanPen/GridCal "GitHub - SanPen/VeraGrid: VeraGrid, a cross-platform power systems software written in Python with user interface, used in academia and industry. · GitHub"))

### Phase 17 — Power systems

```text
bus/branch network
transformers
generators
loads
shunts
HVDC

AC Newton-Raphson power flow
fast-decoupled load flow
DC load flow

3-phase unbalanced power flow
short-circuit analysis

optimal power flow
state estimation
continuation power flow

PTDF / LODF
N-1 contingencies

RMS dynamics
electromechanical generator models
EMT bridge
```

And an especially interesting unified chain becomes:

```text
Power grid
    ↓
Power electronics/SPICE
    ↓
Cable/network solver
    ↓
3D EM
    ↓
Thermal
```

That is a combination conventional solvers rarely put behind one API.

---

## 7. OpenGeoSys + PyLith — geological and porous THMC physics

You already have Darcy, Brinkman, porous transport, chemistry and mechanics. That is almost the entire foundation needed for a serious geological multiphysics layer.

**OpenGeoSys** focuses on coupled **thermo-hydro-mechanical-chemical (THMC)** processes in porous and fractured media, with applications such as geothermal systems, CO₂ storage and subsurface transport. ([GitHub](https://github.com/ufz/ogs/blob/master/README.md?utm_source=chatgpt.com "ogs/README.md at master · ufz/ogs · GitHub"))

**PyLith** adds earthquake/fault-oriented crustal deformation, fault interfaces and friction models, quasi-static and dynamic mechanics. ([GitHub](https://github.com/geodynamics/pylith?utm_source=chatgpt.com "GitHub - geodynamics/pylith: PyLith is a finite element code for the solution of dynamic and quasi-static tectonic deformation problems. · GitHub"))

### Phase 18 — Geomechanics / subsurface

```text
unsaturated Darcy/Richards flow
multiphase porous flow
poroelasticity
Biot equations

fracture/interface flow
matrix-fracture exchange

reactive porous transport
mineral precipitation/dissolution

THM
THMC

fault interfaces
slip weakening
rate-and-state friction

reservoir deformation
geothermal
CO2 sequestration
groundwater
```

This is another field where much of your core already exists.

---

## 8. OpenRadioss — explicit crash/impact mechanics

Your FEM side has structural mechanics, plasticity, contact and dynamics, but high-speed explicit nonlinear mechanics deserves its own execution model.

**OpenRadioss** is aimed at explicit dynamic event simulation and handles crash/impact-style nonlinear mechanics. It is AGPL-licensed, so for your MIT clean-room core I would use it strictly as a capability/validation reference rather than transplanting implementation. ([GitHub](https://github.com/OpenRadioss/OpenRadioss?utm_source=chatgpt.com "GitHub - OpenRadioss/OpenRadioss: OpenRadioss is a powerful, industry-proven finite element solver for dynamic event analysis · GitHub"))

Add:

```text
central-difference explicit FEM
mass lumping
stable timestep / wave-speed criterion

large strain
finite rotation

advanced plasticity
damage/failure

element deletion
hourglass control
shell elements

high-speed contact
self-contact

impact
crash
blast
penetration
```

It fits naturally on top of the existing structural core.

---

## 9. Molecular dynamics / atomistics

This introduces a completely new spatial scale.

**HOOMD-blue** is particularly interesting for your architecture because it is explicitly GPU-oriented and supports molecular dynamics and Monte Carlo with pair, bond and angular interactions; it uses a BSD license. ([GitHub](https://github.com/glotzerlab/hoomd-blue?utm_source=chatgpt.com "GitHub - glotzerlab/hoomd-blue: Molecular dynamics and Monte Carlo soft matter simulation on GPUs. · GitHub"))

**LAMMPS** should be the larger behavioral capability map. It has enormously broad classical MD support and portable accelerator work including Kokkos/SYCL, but its GPL source should remain outside the clean-room implementation. ([GitHub](https://github.com/lammps/lammps?utm_source=chatgpt.com "GitHub - lammps/lammps: Public development project of the LAMMPS MD software package · GitHub"))

### Phase 19 — Atomistics

```text
Verlet neighbor lists
cell-linked neighbor search

Velocity Verlet
Langevin
Brownian dynamics

NVE
NVT
NPT

Lennard-Jones
Morse
Coulomb
bond/angle/dihedral

Ewald
PPPM

energy minimization
molecular constraints

Monte Carlo
GPU particle sorting
domain decomposition
```

Eventually it could couple atomistic simulations to FEM or continuum thermal models.

---

## 10. Quantum ESPRESSO — electronic structure

This is the furthest jump in complexity, but it would close the physics stack all the way down to electrons.

Quantum ESPRESSO provides plane-wave/pseudopotential DFT, Car-Parrinello/ab-initio molecular dynamics, phonons via density-functional perturbation theory and associated electronic-structure tools. It remains GPL-based, so it should be a clean-room scientific reference only. ([GitHub](https://github.com/QEF/q-e?utm_source=chatgpt.com "GitHub - QEF/q-e: Mirror of the Quantum ESPRESSO repository. Please do not post Issues or pull requests here. Use gitlab.com/QEF/q-e instead. · GitHub"))

I would leave this until the lower levels are mature because the prerequisite infrastructure is substantial:

```text
FFT distributed backend
complex eigensolvers
plane-wave bases
pseudopotentials
SCF iteration
density mixing
Kohn-Sham solver
k-point integration
DFPT
```

But it would eventually create an impressive chain:

```text
DFT → atomistics → material properties
                    ↓
             FEM / CFD / EM
                    ↓
             device/system model
```

---

# Three horizontal solver projects also worth mining

These aren't new physics domains, but they may give more value per development hour than several new solvers.

**SUNDIALS** should be a major reference for your common time-integration layer. Its current solver family includes CVODE/CVODES, IDA/IDAS, ARKODE and KINSOL for stiff ODEs, DAEs, IMEX methods, nonlinear systems and sensitivities, with CPU/distributed/GPU support and a BSD license. ([GitHub](https://github.com/llnl/sundials?utm_source=chatgpt.com "GitHub - llnl/sundials: Official development repository for SUNDIALS - a SUite of Nonlinear and DIfferential/ALgebraic equation Solvers. Pull requests are welcome for bug fixes and minor changes. · GitHub"))

That suggests adding a common API roughly like:

```cpp
TimeIntegrator
 ├── ExplicitRK
 ├── DIRK
 ├── IMEX
 ├── BDF
 ├── VariableOrderBDF
 └── Multirate

NonlinearSolver
 ├── Newton
 ├── NewtonKrylov
 └── FixedPoint

SensitivitySolver
 ├── Forward
 └── Adjoint
```

**hypre** is similarly useful for distributed sparse linear algebra and multigrid/preconditioning; its licensing is friendly to this kind of architecture. ([GitHub](https://github.com/hypre-space/hypre?utm_source=chatgpt.com "GitHub - hypre-space/hypre: Parallel solvers for sparse linear systems featuring multigrid methods. · GitHub"))

And **AMReX** is worth studying for block-structured AMR, heterogeneous accelerator execution and scalable mesh infrastructure. That directly attacks still-open FVM and particle-AMR items in your tracker. ([GitHub](https://github.com/AMReX-Codes/amrex?utm_source=chatgpt.com "GitHub - AMReX-Codes/amrex: AMReX: Software Framework for Block Structured AMR · GitHub"))

## What I would actually implement next

Given the current **608/719-capability** repository state documented in the README, I would expand in this order rather than add still more overlapping CFD/EM implementations.

```text
Phase 13  Multibody + rigid bodies + DEM
          Reference: Project Chrono

Phase 14  Acoustics + ultrasound + photoacoustics
          Reference: k-Wave

Phase 15  Semiconductor TCAD
          Reference: DEVSIM

Phase 16  Energy/nuclear
          16A PyBaMM-class batteries
          16B OpenMC-class neutronics

Phase 17  Electrical power systems
          Reference: pandapower + VeraGrid

Phase 18  Geomechanics / porous THMC
          Reference: OpenGeoSys + PyLith

Phase 19  Atomistics / molecular dynamics
          Reference: HOOMD-blue + LAMMPS

Phase 20  Explicit impact/crash mechanics
          Reference: OpenRadioss

Phase 21  Electronic structure / DFT
          Reference: Quantum ESPRESSO
```

Meanwhile I would treat **WarpX as an advanced PIC reference rather than a new solver phase**. Your Phase 12 already has 1-D/2-D/3-D electrostatic and electromagnetic PIC, relativistic Boris/Vay pushers, MCC/plasma chemistry, particle boundaries and distributed-exchange foundations. WarpX is more useful for filling the remaining sophisticated PIC features—AMR, dynamic load balancing, GPU scaling and related algorithms—than for defining another field. Its modern implementation supports CPU/GPU execution, mesh refinement and advanced plasma-accelerator workflows. ([GitHub](https://github.com/BLAST-WarpX/warpx "GitHub - BLAST-WarpX/warpx: WarpX is an advanced Particle-In-Cell code. · GitHub"))

The **three additions with the strongest immediate architectural payoff are Chrono-class multibody/DEM, k-Wave-class acoustics, and DEVSIM-class semiconductor TCAD**. Together they add three fundamentally different solver families—constrained mechanical dynamics, pseudospectral wave propagation and semiconductor carrier transport—while reusing a large fraction of the runtime, sparse algebra, geometry, particle, thermal, electrostatic and multiphysics infrastructure you have already built.