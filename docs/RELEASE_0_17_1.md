# cfd_solvers v0.17.1 - acoustic inverse reconstruction and multiphysics coupling

v0.17.1 completes the machine-counted Phase 14 acoustics roadmap by connecting the v0.17.0 k-space backend to inverse reconstruction, thermal, structural, flow, piezoelectric and optimization workflows.

## Added

- time-reversal photoacoustic reconstruction by enforcing recorded sensor pressure in reversed time order as a Dirichlet condition on the sensor mask;
- acoustic pressure integration to 2-D FEM boundary nodal loads;
- pressure-trace driving of the existing dynamic 1-D structural bar solver;
- direct volumetric-heating input in the Pennes bioheat solver and acoustic intensity/absorption coupling;
- absorption/radiation-force acceleration helper plus direct body-acceleration source support in the collocated incompressible FVM solver;
- coupled 1-D piezoelectric finite-element receive/transmit solver with open-circuit receive mode;
- normalized sensor-trace L2 misfit and acoustic parameter fitting through the existing inverse coordinate-search workflow.

## Validation

The focused regression verifies time-reversal refocusing, integrated structural traction, dynamic acoustic loading, analytic uniform Pennes heating, acoustic-driven FVM acceleration, closed-form piezoelectric open-circuit voltage, transmit displacement, and recovery of a synthetic acoustic parameter.

The complete default project matrix passes **142/142 CTest targets** and Python ABI **0.17.1**. Phase 14 is **21/21 = 100%** and the overall tracker is **670/777 = 86.2%**.

## Scope limits

- time reversal currently uses the same 2-D model/grid as the forward model; regularized attenuation compensation is not yet implemented;
- structural coupling is one-way pressure loading and does not yet feed structural motion back into the acoustic grid;
- FVM coupling is a one-way mean radiation-force/body-acceleration baseline rather than a fully resolved acoustic-streaming model;
- piezoelectric FEM is currently a 1-D small-strain thickness-mode baseline;
- inverse fitting uses derivative-free coordinate search; acoustic adjoints and large-scale tomography remain future work.
