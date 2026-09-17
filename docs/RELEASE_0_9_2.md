# cfd_solvers 0.9.2

## Scope

v0.9.2 continues the CST-class expansion above the v0.9.1 3-D RF-FEM/cable/particle foundation. This release focuses on practical system workflows that sit around the field solvers: signal-integrity eye analysis, power-integrity impedance/IR-drop/decoupling, standardized EMC waveform/probe helpers, Vay/MCC particle dynamics and heterogeneous voxel SAR transfer into the bioheat solver.

## Signal integrity

`cfd::em::analyze_nrz_eye` folds a uniformly sampled waveform by unit interval when the ideal symbol stream is known. It reports:

- zero/one decision means;
- zero/one standard deviations;
- decision threshold;
- eye height using a 3-sigma vertical margin;
- eye width in UI from worst-case phase separation;
- a Gaussian BER estimate.

The regression builds a deterministic NRZ waveform with small bounded jitter/noise and checks level separation, eye height, eye width and a very low BER estimate.

## Power integrity

The new PDN helpers cover two useful early-design tasks.

`pdn_parallel_impedance_ohm` combines ESR/ESL/C decoupling capacitor models and an optional plane impedance into a parallel impedance response. `optimize_decoupling_greedy` chooses capacitor library entries that reduce the worst target-impedance ratio across a frequency grid.

`solve_pdn_ir_drop_grid` solves a loaded 2-D resistive rail/grid with fixed voltage-source nodes. Positive load currents consume rail current, and the result reports the node voltages, minimum voltage, maximum drop and convergence metadata.

## EMC waveform/probe utilities

The EMC utilities add standardized source/probe building blocks:

- peak-normalized double-exponential pulses for ESD/lightning-style current or voltage waveforms;
- damped sine bursts for BCI/ringing-style stimuli;
- peak, RMS, impulse and energy metrics for sampled probe records.

These are intentionally waveform-level utilities; full certification-specific fixtures and probe libraries remain future work.

## Vay pusher and Monte-Carlo collisions

`vay_push` adds a proper-velocity Vay relativistic update for high-gamma crossed-field regimes. The focused regression verifies magnetic-only gamma conservation.

`apply_monte_carlo_collisions` adds a deterministic neutral-collision baseline:

- elastic events isotropically scatter velocity at fixed kinetic energy;
- ionization events remove a configured threshold energy and append a low-energy secondary macro-particle;
- statistics count elastic and ionization events and accumulated energy loss.

The tests force high-probability elastic and ionization events and check energy preservation/loss semantics.

## Heterogeneous voxel SAR and Pennes projection

`VoxelTissueGrid3D` stores per-voxel material indices and RMS electric field. `voxel_sar_w_per_kg` computes local SAR using each tissue's conductivity and density. `max_mass_averaged_sar_w_per_kg` computes local mass-averaged SAR by accumulating nearest voxels around every candidate centre, and `project_voxel_sar_to_pennes2d` mass-weights heterogeneous 3-D SAR into a 2-D Pennes mesh.

Validation checks material-dependent SAR, local 1 g averaging, stacked-tissue projection and subsequent implicit Pennes heating.

## Validation status

- normal CPU CTest matrix: **56/56 pass**;
- focused v0.9.2 ASan+UBSan executable: **pass** with leak detection for the new workflow slice;
- integration tracker: **495/638 (77.6%)**;
- Phase 12: **29/50 (58.0%)**.

## Explicit limitations

This release does not claim production CST or SI/PI/EMC certification parity. In particular, the following remain open:

- self-consistent electromagnetic Yee-grid PIC and charge-conserving current deposition;
- plasma chemistry/ionization source coupling and multipactor/corona threshold workflows;
- temperature-dependent dielectric/perfusion feedback to the EM solve;
- installed-antenna/co-site hybrid coupling;
- full compliance-mask/channel simulation and IBIS/AMI-style SI flows;
- detailed PDN package/PCB geometry extraction and optimizer constraints;
- RWG surface MoM, MLFMM and asymptotic SBR/physical-optics solvers.
