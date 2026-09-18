# cfd_solvers v0.18.1

## Phase 15 transient / AC / C-V analysis

This release extends the clean-room 1-D semiconductor drift-diffusion solver with time-domain and terminal small-signal workflows.

### Added

- backward-Euler electron/hole continuity with quasi-static Poisson/Gummel coupling;
- terminal conductive plus displacement-current reporting;
- terminal small-signal conductance, capacitance, complex admittance and impedance;
- C-V bias sweeps from symmetric electrode-charge differentiation;
- left-terminal electrostatic charge diagnostic;
- silicon Varshni bandgap helper;
- temperature-scaled intrinsic density, mobility and density of states;
- integrated optional temperature dependence in `SemiconductorDevice1D::set_temperature`;
- v0.18.1 analytic regression coverage.

### Scope boundary

The small-signal implementation evaluates a terminal/quasi-static linearization `Y = dI/dV + j omega dQ/dV`. It does not yet assemble and solve the full complex distributed Poisson/electron/hole Jacobian, so high-frequency transit-time and distributed carrier dynamics are outside this release.

### Validation

- complete default matrix: **144/144**;
- Python ABI: **0.18.1**;
- uniform-bar conductance agrees with `q mu n A/L`;
- uniform-bar capacitance agrees with `epsilon A/L`;
- backward-Euler voltage step reproduces `C DeltaV/Delta t`;
- temperature scaling gives lower bandgap/mobility and higher intrinsic density at elevated temperature;
- modified TCAD source/test compile with warnings-as-errors.

### Tracker

Phase 15: **23/30 = 76.7%**. Overall: **693/807 = 85.9%**.
