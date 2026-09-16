# LBM forcing and physical boundaries

Phase 2B adds a validated physical-boundary layer to the D2Q9 CPU reference solver and a generic body-force term to the optimized in-place CPU/SYCL kernels.

## Body force

`D2Q9Solver`, `EsotericPullSolver` and `OneStepPullSolver` use a Guo-style forcing source. The configuration stores acceleration rather than force density; internally the force density is `rho * acceleration`. Macroscopic velocity includes the half-step force correction.

This is useful for pressure-gradient-equivalent periodic channel flow without introducing inlet/outlet compressibility artifacts.

## Halfway bounce-back

A fluid node pulling from a solid neighbour reflects the population from the opposite direction at the fluid node. The wall therefore lies halfway between the fluid and solid cell centres.

The physical wall position matters when comparing against analytical solutions. For the channel regression with solid rows `y=0` and `y=ny-1`, the fluid height is `ny-2` lattice units and the first/last fluid-node coordinates are one half lattice unit from the walls.

## Moving walls

`set_wall_velocity()` marks a cell solid and stores a wall velocity. Incoming bounced populations receive the standard moving-wall momentum correction. The lid-driven cavity regression uses this path.

`set_solid()` explicitly means a stationary wall and therefore resets any wall velocity previously stored on that cell. This prevents stale moving-wall state when a case reuses or reclassifies boundary cells.

## Velocity inlet / pressure outlet

The D2Q9 CPU reference implements left velocity and right density/pressure boundaries using Zou-He reconstruction:

```cpp
solver.set_velocity_inlet_left(0.02F, 0.0F);
solver.set_pressure_outlet_right(1.0F, 0.0F);
```

Solid cells take precedence at corners. The current API intentionally keeps the first implementation narrow and explicit instead of introducing a generic boundary object hierarchy into the hot loop.

## Current backend coverage

| Feature | D2Q9 two-grid CPU | in-place CPU | in-place SYCL |
|---|---:|---:|---:|
| periodic streaming | yes | yes | yes |
| Guo body acceleration | yes | yes | implemented; accelerator validation pending |
| halfway stationary wall | yes | not yet | not yet |
| moving wall | yes | not yet | not yet |
| left velocity inlet | yes | not yet | not yet |
| right pressure outlet | yes | not yet | not yet |

The two-grid D2Q9 implementation is therefore the physical-boundary oracle for the next optimization pass. Boundary-aware in-place kernels should be added only after reproducing these regression results. Body-force parity is already checked on D2Q9, D3Q19 and D3Q27 single-grid CPU paths.
