# cfd_solvers v0.16.7 - conservative CFD/DEM coupling

v0.16.7 connects the Phase-13 rigid-body/DEM family to the existing finite-volume flow stack with both unresolved many-particle coupling and resolved surface-traction coupling. The implementation is deliberately conservative: every fluid-to-particle force has an equal-and-opposite Eulerian reaction, and the resident SYCL path preserves that contract without staging complete fields through host memory.

## Unresolved CFD/DEM

- conservative spherical particle-volume projection to Eulerian cells;
- solid-volume and fluid-void-fraction fields;
- Schiller-Naumann drag with void-fraction hindrance;
- optional pressure-gradient force;
- optional Saffman lift force;
- per-particle force output;
- equal-and-opposite Eulerian force-density deposition;
- direct rigid-body force application;
- exact action/reaction conservation for the discrete nearest-cell projection before any physical void-fraction clipping.

## Resolved rigid-body loading

- arbitrary surface quadrature points with position, outward normal and area;
- pressure traction `-p n`;
- symmetric viscous-stress traction `tau n`;
- integrated resultant force and moment about an arbitrary center;
- direct rigid-body force/torque application;
- six-point octahedral sphere quadrature baseline that exactly cancels uniform pressure and exactly integrates the resultant of a linear pressure field.

## Resident SYCL coupling

`ResidentCfdDemSyclCoupler` reuses the existing resident `PolyMesh` and `ResidentDemSycl` state:

1. project particle volume to device-resident Eulerian solid fraction;
2. derive device-resident void fraction;
3. sample resident fluid velocity/pressure-gradient/vorticity at the nearest cell;
4. accumulate drag/pressure/lift into the resident DEM force arrays;
5. accumulate the opposite force in Eulerian cells;
6. convert that reaction to a fluid acceleration source;
7. feed the source directly to `ResidentIncompressibleSycl` momentum assembly.

`ResidentIncompressibleSycl` now exposes a persistent external-acceleration device field. The ordinary momentum RHS consumes it directly, so the CFD/DEM source does not require a CPU field reconstruction.

## Validation

- complete default matrix: **139/139 CTest targets passed**;
- Python ABI: **0.16.7**;
- resolved single-surface pressure force and torque;
- uniform pressure cancellation on a closed sphere;
- conservative particle-volume projection;
- CPU two-way particle/fluid momentum conservation;
- CPU-executing fake-SYCL resident coupling regression;
- resident particle acceleration from fluid drag;
- resident equal-and-opposite fluid acceleration source;
- zero FVM/DEM bulk host traffic inside the covered resident coupling step;
- strict `-Wall -Wextra -Wpedantic -Werror` compilation of the modified SYCL sources.

Physical GPU qualification is still separate from compile/emulation validation.

## Tracker impact

Phase 13 closes:

- resolved CFD <-> rigid-body surface traction/force/torque coupling;
- unresolved CFD <-> many-particle drag/void-fraction coupling.

This moves Phase 13 from 30/37 to **32/37 = 86.5%**, and overall machine-counted progress from 643/756 to **645/756 = 85.3%**. The distributed-memory DEM checkbox remains open until real multi-rank GPU-aware MPI execution is qualified.
