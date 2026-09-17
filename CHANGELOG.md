# Changelog

## 0.15.2 - transported RANS, Reynolds stress and SST-DES

- Added shared implicit finite-volume transport for turbulence variables with variable diffusion, upwind convection, semi-implicit sinks, fixed/zero-gradient boundaries and Euler/BDF2/Crank-Nicolson time integration.
- Added transported Spalart-Allmaras, standard k-epsilon and k-omega SST models with positivity and production safeguards.
- Added SST F1/F2 blending, cross diffusion and an integrated DES dissipation switch with optional F1/F2 zonal shielding.
- Added six-component Reynolds-stress transport plus epsilon, LRR-style pressure-strain redistribution and realizability projection.
- Added focused v0.15.2 turbulence regressions and advanced Phase 3 to 94/106 tracked capabilities.

## 0.15.0 - FVM second-order time integration and pressure AMG seam

- Integrated BDF2/backward and off-centered Crank-Nicolson into production scalar FVM assembly.
- Integrated second-order physical time stepping into PISO/PIMPLE momentum prediction while retaining SIMPLE pseudo-time semantics.
- Added matrix-free pressure PCG with Jacobi and external AMG/V-cycle preconditioner paths.
- Added focused temporal-accuracy, PISO/PIMPLE history and AMG-callback regressions.
- Closed three Phase-3 tracker items.

## 0.14.2

- Added device-keyed kernel fusion/autotuning database and benchmark selection API.
- Added multi-device accelerator groups/work slices per local rank plus SYCL queue groups.
- Added weighted adaptive contiguous repartition planning and migration-segment analysis.
- Added self-hosted NVIDIA/AMD/Intel GPU workflow definitions without claiming hardware qualification before successful runs.
- Added focused portability/autotuning regression and updated integration tracker to 597/719.


## 0.14.1 - LBM physics, geometry, precision and diagnostics wave

- Added a D2Q9 cumulant collision baseline and Smagorinsky-Lilly LES option.
- Added conservative VOF/free-surface transport, interface geometry and capillary-force coupling.
- Added interpolated curved-wall bounce-back and passive/two-way immersed-boundary particles.
- Added local LBM acceleration fields and optional SYCL triangle-mesh voxelization.
- Added 16-bit compressed population storage with FP16/BF16/shifted/custom codec experiments and explicit accuracy gates.
- Added 3-D Q-criterion output and a standalone interactive HTML scalar viewer.
- Advanced Phase 2 from 29/46 to 40/46 tracked capabilities.

## 0.14.0 - Phase 1 common HPC runtime completion

- Added NUMA topology discovery, compact/spread/explicit-node thread placement plans, per-thread CPU binding and placement-aware first touch.
- Added dynamic Linux affinity masks so systems above CPU 1023 do not silently fall back to fixed `cpu_set_t` limits.
- Added persistent-device SYCL CSR SpMV plus a device-resident conjugate-gradient Krylov baseline using USM and SYCL reductions.
- Added MPI ownership discovery and cached sparse-vector request exchange integrated directly into distributed CG iterations.
- Added focused v0.14.0 CPU/SYCL tests and extended the MPI regression with a distributed tridiagonal solve.
- Completed Phase 1 at 34/34 tracked capabilities.

## 0.13.0 - Phase 9 interoperability completion

- Added inline XDMF Tri3 mesh and nodal-scalar output with deterministic validation.
- Added optional native HDF5 datasets and XDMF heavy-data sidecars through `CFD_ENABLE_HDF5`.
- Added a stable shared C ABI and a pure-stdlib `ctypes` Python package.
- Added C++ boundary and end-to-end Python binding smoke tests.
- Completed the Phase 9 interoperability/workflow/UX tracker at 50/50.

## 0.12.0 - production controller deployment baseline

- Added TLS 1.2+ server wrapping with required operator-supplied certificate and key paths.
- Added SHA-256 bearer-token identities, constant-time verification and viewer/operator/admin roles.
- Added bounded rotating NDJSON audit history and authenticated live deployment probes.
- Added hardened systemd and Caddy reverse-proxy artifacts plus secret-free environment templates.
- Added focused production-controller regression coverage and generated Python syntax validation.

## 0.11.8 - server-sent controller status events

- Added `GET /api/events/stream` with SSE status events, monotonic per-connection IDs and configurable heartbeat comments.
- Updated the generated dashboard to use `EventSource` with automatic reconnect while retaining periodic polling as a compatibility fallback.
- Added controller stream configuration validation, OpenAPI/route-table exposure and focused v0.11.8 regression coverage.
- Added release, deployment, learning and clean-room notes for the live event seam.

## 0.11.7 - generated multi-server controller dashboard

- Added generated browser dashboard assets: `controller/index.html`, `controller/dashboard.js` and `controller/dashboard.css`.
- Added controller routes for `/`, `/dashboard.js`, `/dashboard.css` and `/api/events`.
- Added polling dashboard UI with summary cards, case table, stdout/stderr log preview and token-prompted control buttons.
- Added newline-delimited `controller/events.ndjson` event snapshot for downstream dashboards.
- Added configurable dashboard title and refresh cadence to `RemoteControllerConfig`.
- Added CLI smoke case `particle-multiserver-dashboard` and focused regression target `cfd-v0117-controller-dashboard-tests`.
- Updated deployment docs, SU2/csauto notes, study guide, learning resources, README, release notes and integration tracker.

## 0.11.6 - multi-server controller API scaffold

- Added generated dependency-free Python controller service scaffold for multi-server supervision artifacts.
- Added controller routes for health, status, case listing and stdout/stderr log tails.
- Added optional token-gated mutation route that writes `.cfd_control/<action>.directive` files for approved control actions.
- Added generated OpenAPI JSON, status JSON, route table, environment example, controller README and `multiserver_controller.sh` launcher.
- Added read-only controller planning mode for dashboard-only deployments.
- Added CLI smoke case `particle-multiserver-controller` and focused regression target `cfd-v0116-multiserver-controller-tests`.
- Updated deployment docs, SU2/csauto notes, study guide, learning resources, README, release notes and integration tracker.

## 0.11.5 - multi-server supervision hardening

- Added multi-server supervision plans on top of supervised remote job metadata.
- Added SSH/local access probes, Docker daemon probes and campaign-root writability checks.
- Added retry-wrapped remote launch command plans with configurable attempt count and delay.
- Added stdout/stderr tail command plans for dashboard/log-panel integration.
- Added sensitive command-display redaction for inline TOKEN/PASSWORD/SECRET/KEY-style values.
- Added persisted `multiserver_access_checks.sh`, `multiserver_retry_launch.sh`, `multiserver_tail_logs.sh`, `multiserver_dashboard.json` and `multiserver_supervision.tsv`.
- Added CLI smoke case `particle-multiserver-supervision` and focused regression target `cfd-v0115-multiserver-supervision-tests`.
- Updated deployment docs, SU2/csauto notes, study guide, learning resources, README, release notes and integration tracker.

## 0.11.4 - supervised multi-server execution and Docker healthchecks

- Added supervised multi-server execution planning on top of the v0.11.3 placement layer.
- Added remote host health-check command plans for campaign root, native/Docker runtime availability and disk-space inspection.
- Added per-case remote job plans with launch, status, cancel and fetch-log commands, plus `.cfd_run` stdout/stderr/pid/exit-code paths.
- Added persisted `multiserver_health.sh`, `multiserver_launch.sh`, `multiserver_status.sh`, `multiserver_fetch_logs.sh` and `multiserver_jobs.tsv` artifacts.
- Added remote job status table parsing and registry update helpers for running/done/failed/stopped distributed cases.
- Extended Docker deploy scaffolding with Compose healthchecks plus `deploy/healthcheck.sh` and `deploy/worker-entrypoint.sh`.
- Added CLI smoke case `particle-multiserver-execution` and focused regression target `cfd-v0114-multiserver-execution-tests`.
- Updated deployment docs, SU2/csauto notes, study guide, learning resources, README, release notes and integration tracker.

## 0.11.3 - multi-server campaign planning and Docker deployment scaffold

- Added `CampaignServerDescriptor` and multi-server campaign placement with server slots, online filtering, tag filtering and deterministic assignment.
- Added generated remote command plans for directory creation, rsync-style case synchronization and native/Docker solver launch commands.
- Added persisted `multiserver_assignments.tsv` and `multiserver_commands.sh` plan files for auditable distributed campaign launches.
- Added Docker/Compose deployment scaffold generation: Dockerfile, `compose.yaml`, `.dockerignore`, `deploy/deploy.sh`, `deploy/env.example`, `deploy/servers.example.tsv` and deployment README.
- Added CLI smoke case `particle-multiserver-deploy` and focused regression target `cfd-v0113-multiserver-deploy-tests`.
- Updated SU2/csauto clean-room notes, deployment docs, study guide, learning resources, README, release notes and integration tracker.

## 0.11.2 - campaign control directives and scheduler status refresh

- Added live-control directive files for solver-adapter-approved stop, extend, checkpoint and flush actions.
- Added directive discovery for generated case folders under `.cfd_control`.
- Added Slurm-style scheduler queue parsing and registry status application for queued/running/done/failed/cancelled jobs.
- Added campaign output refresh that scans discovered logs, residual histories and performance files to update registry status/objective/iteration fields without relaunching cases.
- Added CLI smoke case `particle-campaign-control` and focused regression target `cfd-v0112-campaign-control-tests`.
- Updated SU2/csauto clean-room notes, study guide, learning resources, README, release notes and integration tracker.

## 0.11.1 - local campaign runner and runtime doctor checks

- Added runtime doctor checks for solver adapter descriptors, case directories and launch executables across native/container/scheduler boundaries.
- Added POSIX native argv process launching for local solver runs with captured stdout/stderr logs and without shell command-string execution.
- Added campaign output discovery for residual/history and performance/timing files.
- Added a local campaign runner that executes pending cases, parses residual/performance outputs and updates `registry.tsv` objective/iteration/status fields.
- Added CLI smoke case `particle-campaign-local-runner` and focused regression target `cfd-v0111-local-campaign-runner-tests`.
- Updated SU2/csauto clean-room notes, study guide, learning resources, README, release notes and integration tracker.

## 0.11.0 - persistent campaign execution and external-solver workflow boundary

- Added persistent campaign folder generation with rendered per-case files, `doe_row.csv` snapshots, template manifests and registry TSV roundtrips.
- Added external-solver command planning for native, Docker, Apptainer/Singularity and Slurm-style wrappers while keeping actual process launch outside this safe reference layer.
- Added residual-history and performance-metric log parsers plus a finite-difference gradient-descent campaign optimization loop.
- Added CLI smoke case `particle-campaign-execution` and focused regression target `cfd-v0110-campaign-execution-tests`.
- Updated SU2/csauto clean-room notes, study guide, learning resources, README, release notes and integration tracker.

## 0.10.9 - stochastic transport, dose scoring and SU2/csauto campaign workflow

- Extended the Geant4-inspired particle-through-matter seam with stochastic sampled interaction lengths, physics-list bundles, production-cut overrides, region BVH lookup, sensitive-detector hit collection, dose-grid scoring and SAR projection into the Pennes bioheat solver.
- Added SU2/csauto-inspired workflow utilities: factorial and Latin-hypercube campaign generation, template placeholders and IF/ENDIF conditional rendering, solver-adapter descriptors/capabilities, campaign registry summaries, and finite-difference gradient/descent helpers for design studies.
- Added CLI smoke cases `particle-transport-dose-bvh` and `particle-campaign-doe` plus focused target `cfd-v0109-transport-campaign-tests`.
- Updated clean-room notes, learning resources, CST/SU2/csauto audit material, release notes and the integration tracker.

## 0.10.8 - Geant4-inspired particle-through-matter transport

- Added a clean-room Geant4-inspired track/step/process transport baseline for particle-through-matter workflows.
- Added `TransportTrack`, `TransportStep`, `TransportProcess`, `TransportRegion`, `TransportMaterial`, `TransportWorld` and scoring diagnostics.
- Added geometry/user/process step arbitration, continuous stopping-power losses, discrete secondary-producing interactions, production cuts and energy cuts.
- Added the `particle-geant4-transport` CLI smoke case and `cfd-v0108-geant4-transport-tests` focused regression target.
- Added Geant4 clean-room notes plus study-guide and learning-resource updates.

## 0.10.7 - distributed staggered PIC timestep bridge

- Added `run_serialized_distributed_staggered_pic_step_3d`, a rank-local distributed 3-D EM-PIC timestep bridge that pushes particles with domain-local electromagnetic fields, migrates ownership through the serialized exchange contract and returns six-component EM guard fields for the next local update.
- Added `DistributedStaggeredPicStep3DConfig` / `DistributedStaggeredPicStep3D` diagnostics for particle counts, displacement, serialized exchange accounting and guarded fields.
- Added the `particle-distributed-step3d` CLI smoke case and `cfd-v0107-distributed-step-tests` regression target.
- Updated the CST audit, study guide and learning resources with the distributed timestep ordering.

## 0.10.6 - distributed PIC exchange round and EM field guard vectors

- Added six-component electromagnetic field-vector guard exchange for Ex/Ey/Ez/Bx/By/Bz domain blocks.
- Added an end-to-end serialized rank-local distributed PIC exchange round that packs particle migration and scalar guard payloads, addresses logical ranks, serializes/deserializes envelopes and applies the reconstructed payloads.
- Added `particle-distributed-round3d` and `particle-field-guards3d` CLI smoke cases.
- Added `cfd-v0106-distributed-round-tests` focused regression/sanitizer coverage.
- Updated the CST/PIC roadmap, study guide and learning resources around field guard exchange and distributed exchange-round orchestration.

## 0.10.5 - serialized PIC transport and MPI bridge contract

- Added ABI-local serialized particle/guard transport envelopes for 3-D PIC domain exchange.
- Added deterministic serialization/deserialization, destination-rank exchange planning, payload diagnostics and malformed-envelope rejection tests.
- Added a guarded MPI allgather transport wrapper for the serialized envelope contract (`CFD_HAS_MPI` builds); real MPI runtime validation remains a hardware/CI item.
- Added `particle-serialized-transport3d` CLI smoke coverage and focused regression/sanitizer coverage.
- Added multi-species plasma mass-action chemistry utilities with charge-density diagnostics.
- Added Paschen gas-breakdown and parallel-plate multipactor threshold estimators with CLI/regression coverage.
- Updated CST/PIC roadmap documentation and study resources around transport contracts, MPI readiness and remaining production gaps.

## 0.10.4 - in-memory PIC transport contract for future MPI

- Added deterministic rank-topology mapping from PIC subdomains to logical communicator ranks.
- Added rank-addressed transport envelopes for particle-migration and scalar guard-cell messages.
- Added `InMemoryPicTransport3D` and a complete exchange helper that reproduces the serial migration/guard result through a send/receive-like contract.
- Added `particle-transport3d` CLI smoke case and `cfd-v0104-transport-tests` focused regression target.
- Expanded Phase 12 from 58/78 to 61/81 validated capabilities and the whole tracker to 527/669 (78.8%).
- Expanded the debug/CI-speed CPU validation matrix to 80 CTest targets; all pass. A focused v0.10.4 ASan+UBSan run with leak detection passes for the new in-memory transport slice.

## 0.10.3 - communicator-ready PIC migration and guard-message packets

- Added communicator-ready particle migration messages with retained-domain buckets, source-domain/source-index tracking, periodic wrapping and nonperiodic outside-particle records.
- Added communicator-ready scalar guard-cell messages that carry destination padded coordinates and reproduce the serial ghost-padded guard exchange when applied.
- Added `particle-comm-exchange3d` CLI smoke case and `cfd-v0103-comm-exchange-tests` focused regression target.
- Expanded Phase 12 from 58/78 to 61/81 validated capabilities and the whole tracker to 527/669 (78.8%).
- Expanded the debug/CI-speed CPU validation matrix to 78 CTest targets; all pass. A focused v0.10.3 ASan+UBSan run with leak detection passes for the new communicator-message slice.

## 0.10.2 - PIC domain decomposition and scalar guard exchange foundation

- Added Cartesian 3-D PIC domain descriptors with uneven global-cell block coverage and physical extents.
- Added serial particle migration bucketing between PIC subdomains, including periodic wrapping and nonperiodic outside-particle reporting.
- Added serial scalar guard-cell exchange that builds ghost-padded local blocks from neighbouring subdomains, with deterministic periodic and nonperiodic exterior handling.
- Added `particle-domain-exchange3d` CLI smoke case and `cfd-v0102-domain-exchange-tests` focused regression target.
- Expanded Phase 12 from 55/75 to 58/78 validated capabilities and the whole tracker to 524/666 (78.7%).
- Expanded the debug/CI-speed CPU validation matrix to 76 CTest targets; all pass. A focused v0.10.2 ASan+UBSan run with leak detection passes for the new domain-exchange/guard-cell slice.

## 0.10.1 - PIC particle sorting and guard-halo decomposition foundation

- Added stable 3-D PIC particle sorting by cell with cell offsets, original-index recovery and occupied-cell diagnostics.
- Added reusable 3-D guard-halo classification for future MPI/GPU particle-domain exchange; edge/corner particles can appear in multiple face lists.
- Added optional particle sorting inside `StaggeredElectromagneticPic3D` with configurable interval and diagnostics.
- Added `particle-sort-halo3d` CLI smoke case and `cfd-v0101-particle-sorting-tests` focused regression target.
- Expanded Phase 12 from 52/72 to 55/75 validated capabilities and the whole tracker to 521/663 (78.6%).
- Expanded the debug/CI-speed CPU validation matrix to 74 CTest targets; all pass. A focused v0.10.1 ASan+UBSan run with leak detection passes for the new particle sorting/guard-halo slice.

## 0.10.0 - local finite-volume 3-D EM-PIC current deposition

- Added `deposit_charge_conserving_current_3d_local`, a finite-volume local current reconstruction that satisfies periodic backward-difference continuity from old/new trilinear CIC charge densities.
- Added a `CurrentDeposition3DMode` selector to `StaggeredElectromagneticPic3DConfig`, allowing the staggered 3-D EM-PIC path to use either spectral continuity reconstruction or the new local finite-volume current path.
- Added configurable polynomial order for the 3-D absorbing sponge boundary profile.
- Added `particle-local-current3d` CLI smoke case and `cfd-v0100-local-current-tests`.
- Expanded Phase 12 from 49/69 to 52/72 validated capabilities and the whole tracker to 518/660 (78.5%).
- Expanded the debug/CI-speed CPU validation matrix to 72 CTest targets; all pass. A focused v0.10.0 ASan+UBSan run with leak detection passes for the new local-current slice.

## 0.9.9 - staggered 3-D Yee electromagnetic PIC foundation

- Added `StaggeredElectromagneticPic3D`, a true 3-D/3-V Yee-style EM-PIC reference with component-specific Ex/Ey/Ez and Bx/By/Bz staggered field locations.
- Added 3-D Yee CFL checking, staggered curl updates, trilinear staggered E/B gather and spectral Jx/Jy/Jz continuity-current coupling.
- Added reusable 3-D electric-wall and absorbing-sponge field-boundary primitives plus absorbing/specular particle-wall handling with secondary-yield reporting.
- Added `particle-staggered-em-pic3d` CLI smoke case and `cfd-v099-staggered-pic3d-tests` regression target.
- Expanded Phase 12 from 45/65 to 49/69 validated capabilities and the whole tracker to 515/657 (78.4%).
- Expanded the debug/CI-speed CPU validation matrix to 70 CTest targets; all pass. A focused v0.9.9 ASan+UBSan run with leak detection passes for the new staggered 3-D EM-PIC slice.

## 0.9.8 - periodic 3-D/3-V electromagnetic PIC baseline

- Added `ElectromagneticPic3D`, a compact periodic 3-D/3-V electromagnetic PIC baseline with full Ex/Ey/Ez and Bx/By/Bz nodal field storage.
- Added centered periodic Maxwell curl updates, a compact 3-D CFL guard, trilinear E/B gather and Vay particle pushing for 3-D macro-particles.
- Reused the v0.9.7 spectral 3-D charge-conserving current reconstruction as Jx/Jy/Jz source coupling into the electromagnetic field update.
- Coupled Monte-Carlo neutral collisions into the 3-D EM-PIC step and preserved secondary macro-particles through the 3-D particle representation.
- Added `particle-em-pic3d` CLI/CTest smoke coverage and focused regressions for bounded vacuum wave energy, current continuity, Poisson correction, field gather, collision coupling and CFL rejection.
- Expanded Phase 12 from 42/62 to 45/65 validated capabilities and the whole tracker to 511/653 (78.3%).
- Expanded the CPU test matrix to 68 CTest targets; all pass in the debug/CI-speed validation build. A focused v0.9.8 ASan+UBSan run with leak detection passes for the new 3-D EM-PIC slice.

## 0.9.7 - 3-D electrostatic PIC foundation

- Added `PicParticle3D` and `ElectrostaticPic3D`, a periodic 3-D electrostatic PIC baseline with trilinear CIC charge deposition, spectral Poisson Ex/Ey/Ez reconstruction, trilinear field gather and particle stepping.
- Added 3-D spectral charge-conserving current reconstruction from old/new particle charge density to validate continuity before full 3-D Yee EM-PIC.
- Added `particle-pic3d` CLI/CTest smoke coverage and focused regressions for 3-D charge conservation, manufactured periodic Poisson fields, current continuity and PIC stepping.
- Expanded Phase 12 from 39/59 to 42/62 validated capabilities and the whole tracker to 508/650 (78.2%).
- Expanded the CPU test matrix to 66 CTest targets; all pass in the debug/CI-speed validation build. A focused v0.9.7 ASan+UBSan run with leak detection passes for the new 3-D PIC slice.

## 0.9.6 - staggered Yee EM-PIC and self-study documentation

- Added `StaggeredElectromagneticPic2D`, a true Yee-style 2-D/3-V electromagnetic PIC baseline with component-specific Ex/Ey/Ez/Bx/By/Bz locations and a 2-D CFL guard.
- Integrated electric-wall and absorbing-sponge field boundaries directly into the staggered EM-PIC update.
- Integrated non-periodic absorbing/specular particle walls with secondary-yield macro-weight reporting into the staggered EM-PIC step.
- Reused the Vay pusher, charge-conserving in-plane current reconstruction, transverse `Jz` deposition and MCC collision coupling inside the staggered solver.
- Added `particle-staggered-em-pic2d` CLI/CTest smoke coverage and focused regressions for bounded vacuum TM energy, tangential-E wall enforcement, sponge damping, particle absorption/reflection, `Jz -> Ez` coupling and MCC coupling.
- Added `docs/IMPLEMENTED_SOLVER_STUDY_GUIDE.md` and `docs/LEARNING_RESOURCES.md` to teach the implemented solver families through local source files, tests, CLI cases and curated external resources.
- Expanded Phase 12 from 36/56 to 39/59 validated capabilities and the whole tracker to 505/647 (78.1%).
- Expanded the CPU test matrix to 64 CTest targets; all pass in the debug/CI-speed validation build. A focused v0.9.6 ASan+UBSan run with leak detection passes for the new staggered EM-PIC slice.

## 0.9.5 - periodic 2-D/3-V electromagnetic PIC baseline

- Added `ElectromagneticPic2D`, a periodic 2-D/3-V electromagnetic PIC solver with centered curl field updates, bilinear field gather, Vay particle push and optional periodic Poisson in-plane correction.
- Reused the v0.9.4 spectral charge-conserving 2-D current reconstruction for in-plane `Jx/Jy` and added CIC transverse `Jz` deposition coupled to the `Ez` update.
- Coupled Monte-Carlo neutral collisions into the 2-D EM-PIC step, including secondary macro-particles mapped back into the 2-D population.
- Added `particle-em-pic2d` CLI/CTest smoke coverage and focused regressions for vacuum TM-wave bounded energy, continuity residuals, `Jz -> Ez` coupling, field gather and collision coupling.
- Expanded Phase 12 from 34/54 to 36/56 validated capabilities and the whole tracker to 502/644 (78.0%).
- Expanded the CPU test matrix to 62 CTest targets; the focused v0.9.5 regression target and `particle-em-pic2d` smoke case pass.

## 0.9.4 - 2-D electrostatic PIC foundation and boundary primitives

- Added periodic 2-D electrostatic PIC with bilinear CIC charge deposition, spectral Poisson electric fields, bilinear field gather and particle stepping.
- Added 2-D spectral charge-conserving current reconstruction from old/new CIC charge density, providing a validated bridge toward full 2-D/3-V EM-PIC deposition.
- Added reusable 2-D particle wall handling plus electric-wall and absorbing-sponge field-boundary primitives for future metallic/absorbing PIC and FDTD boundaries.
- Added `particle-pic2d` CLI/CTest smoke coverage and a focused regression target for manufactured 2-D Poisson fields, continuity residuals, field gather and boundary behavior.
- Expanded Phase 12 from 31/51 to 34/54 validated capabilities and the whole tracker to 500/642 (77.9%).
- Expanded the CPU test matrix to 60 CTest targets; the focused v0.9.4 regression target passes in the normal build.

## 0.9.3 - self-consistent electromagnetic PIC baseline

- Added periodic 1-D/3V self-consistent electromagnetic PIC with Yee-style transverse field update, optional longitudinal Poisson field, Vay particle push, field gather and source deposition.
- Added spectral charge-conserving longitudinal current reconstruction from old/new CIC charge density plus transverse CIC current deposition for particle velocities.
- Coupled deterministic Monte-Carlo elastic/ionization collisions into EM-PIC stepping so secondary macro-particles continue in subsequent self-consistent field updates.
- Added `particle-em-pic1d` CLI/CTest smoke coverage and a focused regression target for continuity residuals, source-free wave energy, beam-current Ampere response and ionization source coupling.
- Expanded Phase 12 from 29/50 to 31/51 validated capabilities and the whole tracker to 497/639 (77.8%).
- Expanded the CPU test matrix to 58 CTest targets; all pass in the debug/CI-speed validation build. A focused v0.9.3 ASan+UBSan run with leak detection passes for the new EM-PIC slice.

## 0.9.2 - SI/PI/EMC workflows, Vay/MCC particles and heterogeneous bioheat

- Added NRZ eye-diagram post-processing with one/zero level statistics, eye height, eye width and Gaussian BER estimate.
- Added PDN impedance, greedy decoupling-capacitor selection and loaded resistive-grid IR-drop solving.
- Added normalized double-exponential and damped-sine EMC waveform sources plus probe peak/RMS/impulse/energy metrics.
- Added a Vay relativistic particle pusher for high-gamma crossed-field regimes.
- Added deterministic Monte-Carlo neutral collision handling with elastic scattering, ionization loss and secondary macro-particle emission.
- Added heterogeneous voxel tissue materials, voxel SAR conversion, local mass-averaged SAR and mass-weighted SAR projection into Pennes bioheat meshes.
- Expanded Phase 12 from 22/50 to 29/50 validated capabilities and the whole tracker to 495/638 (77.6%).
- Expanded the CPU test matrix to 56 CTest targets; all pass. A focused v0.9.2 ASan+UBSan run with leak detection passes for the new workflow slice.

## 0.9.1 - 3-D edge FEM, adaptive RF, cables, relativistic particles and nonlinear magnetics

- Added reusable complex sparse systems through real-block CSR + ILU(0)-GMRES with a sparse GMRES fallback for indefinite harmonic systems.
- Added lowest-order first-kind Tet4 Nedelec/Whitney geometry, basis and curl operators.
- Added a sparse driven 3-D frequency-domain Maxwell FEM with PEC boundary edges, volumetric impressed current and directed edge-current/lumped excitation.
- Added reconstructed 3-D E/H fields and resonator electric/magnetic energy, dielectric loss, conductor-wall loss and Q extraction.
- Added rectangular PEC TE/TM wave-port mode extraction with cutoff, propagation constant, modal impedance and 1 W power normalization.
- Added face-jump RF indicators, Dorfler marking and conforming longest-edge edge-star Tet4 refinement.
- Added stable common-pole rational N-port MOR fitting on top of the adaptive RF sweep.
- Added multiconductor cable/harness frequency-domain RLCG propagation with full matrix terminations, shield transfer impedance and effective-height field coupling.
- Added relativistic Boris particle advance, absorbing/specular wall models, secondary-emission yield, resonator wake generation, bunch-wake convolution and wake impedance.
- Added nonlinear piecewise B-H magnetostatic iteration, stranded-coil flux-linkage/inductance coupling helpers and Maxwell-stress force/torque integration.
- Expanded Phase 12 from 8/47 to 22/50 validated capabilities and the whole tracker to 488/638 (76.5%).
- Expanded the CPU test matrix to 55 CTest targets; all pass. A focused v0.9.1 ASan+UBSan run with leak detection also passes.

## 0.9.0 - CST-class frequency-domain EM, particles/PIC and bioheat

- Added a driven complex 1-D frequency-domain Maxwell/Helmholtz reference solver with PEC boundaries, conductive loss and manufactured-solution validation.
- Added generalized dielectric-loaded 1-D PEC electromagnetic cavity eigenmodes with analytical frequency regression.
- Added reusable lowest-order first-kind Nedelec Tri3 basis/curl utilities and a driven 2-D edge-element Maxwell solver with deterministic global-edge orientation and PEC tangential constraints.
- Added mesh-refinement validation for the 2-D Nedelec Maxwell manufactured solution.
- Added a non-relativistic 3-D Boris charged-particle pusher with long-run uniform-magnetic-field energy conservation validation.
- Added periodic 1-D electrostatic PIC with cloud-in-cell charge deposition, spectral Poisson electric fields, field interpolation and particle stepping.
- Added RMS E-field to SAR conversion and an implicit 2-D Pennes bioheat solver with conduction, perfusion, metabolic heat and spatial SAR sources.
- Added CLI/CTest smoke cases for frequency-domain EM, edge-element EM, electrostatic PIC and SAR/bioheat.
- Added `docs/CST_COVERAGE_AUDIT.md` and Phase 12 to distinguish current solver coverage from still-open CST-class domains such as 3-D RF FEM, MLFMM/SBR, cable/harness EMC, electromagnetic PIC/plasma and wakefields.
- Reconciled stale Phase 10/11 tracker entries against already-tested RF/SPICE implementations; tracker now reports 473/634 capabilities (74.6%).
- Full CPU validation passes 52/52 CTest targets; a focused standalone ASan+UBSan build of the new frequency-domain EM/Nedelec/PIC/bioheat slice passes with leak detection enabled.

## 0.8.2 - unequal-step Gear/BDF2 and arbitrary-orientation wire MoM

- Extended adaptive circuit transient integration with variable-step BDF2/Gear coefficients after a backward-Euler bootstrap interval.
- Added second-order Richardson LTE scaling/controller behavior and an analytical RC refinement regression.
- Added arbitrary-orientation disjoint straight-wire MoM using the projected dyadic free-space Green kernel.
- Added an ideal lossless transformer convenience model from a power-conserving VCVS/CCCS pair.
- Added focused validation for lossy sampled TEM-line attenuation and reconciled already-tested microstrip/stripline/CPW/TE10 AC circuit elements in the tracker.
- Kept the parallel-wire API source-compatible by delegating it to the new orientation-independent solver.
- Added rotational-invariance validation for a rigidly rotated coupled two-wire problem and explicit rejection of touching/crossing wires until junction basis functions are implemented.
- Preserved the full existing v0.8.1 CPU test matrix while advancing the RF/SPICE tracker.

## 0.8.1 - adaptive circuit analysis and coupled-wire RF

- Added adaptive backward-Euler transient integration with Richardson step-doubling local-truncation-error control, accepted/rejected-step accounting and bounded variable timesteps.
- Added periodic steady-state shooting-by-settling with phase-aligned cycle convergence and a low-order complex rational-fit pole/zero baseline.
- Generalized the thin-wire MoM core from a single center-fed straight wire to coupled parallel z-directed wires with independent segment counts, multiple delta-gap feeds and lumped series loads, while retaining the original dipole API.
- Added dedicated analytical/reciprocity regressions and `spice-adaptive`, `spice-pss-pz` and `rf-multiwire` CLI smoke cases.
- Full GCC/OpenMP CPU validation passes 47/47 CTest targets.
- Integration tracker advances to 433/577 validated capabilities (75.0%): Phase 10 is 31/63 and Phase 11 is 64/96.
- Deliberately still open: higher-order variable-step Gear/BDF, exact descriptor generalized-eigenvalue pole/zero, harmonic balance/Newton shooting, NEC-grade arbitrary connected-wire geometry and canonical antenna benchmark convergence.

## 0.8.0 - RF/antenna and SPICE-class solver expansion

- Integrated the latest remote CPU continuation with the locally validated advanced FDTD/multiphysics work without dropping either feature set.
- Added generic two-port/N-port RF conversions, Touchstone SnP, real-reference renormalization, mixed-mode conversion, reference-plane shifts, passivity/reciprocity, K/mu stability, stability circles, transducer gain and load-pull sampling.
- Added microstrip, stripline, CPW, coax, TEM-line and rectangular-waveguide TE10 analytical/network baselines plus antenna-array/polarization utilities.
- Added thin-wire center-fed MoM and PEEC filament extraction baselines, including skin-depth/AC resistance and SPICE coupled-inductor export.
- Added a clean-room MNA circuit engine with R/C/L, independent/controlled sources, switches, mutual inductance, diode, MOS Level-1, BJT, JFET and sampled RF N-port devices.
- Added DC/AC/transient analyses, backward-Euler/trapezoidal/BDF2 integration, homotopy/PN limiting, sweeps, thermal noise integration, sensitivity, Monte Carlo, Fourier/THD and circuit S-parameter extraction.
- Extended SPICE parsing with engineering suffixes, `.MODEL`, `.PARAM`, expressions, nested `.SUBCKT`, scoped local models, `.FUNC`, `.INCLUDE` and `.LIB`.
- Added an OSDI/OpenVAF dynamic-library discovery/ABI-version seam while deliberately leaving descriptor evaluation and full DAE integration for the next compact-model phase.
- Added `docs/RF_CIRCUIT_RESEARCH.md` and expanded the authoritative tracker with Palace/OpenSEMBA/OpenNEC and ngspice/Xyce/QucsatorRF capability families.
- Final tracker state for this checkpoint is 429/576 validated capabilities (74.5%).
- GCC/OpenMP clean matrix passes 43/43 tests; focused Clang and serial RF/SPICE tests pass; focused ASan+UBSan runs pass for SPICE hierarchy/file IO, OSDI loading, PEEC/MoM and nonlinear circuit paths.

## 0.7.1 - CPML, TFSF and lumped FDTD elements

- Added a clean-room 1-D convolutional PML with graded sigma/kappa profiles and convolution memory variables.
- Added a +x 1-D total-field/scattered-field source with Yee temporal and half-cell spatial staggering.
- Added field-coupled parallel lumped R/L/C elements; resistor and capacitor contributions modify local Ampere coefficients while inductor current is advanced as an auxiliary state.
- Added CPML attenuation, homogeneous TFSF leakage, analytical inductor-ramp, resistor-dissipation and capacitor-stability regressions.
- Added CLI/CTest and sanitizer-smoke coverage for the new FDTD paths.
- Integration tracker advanced to 152/398 validated capabilities (38.2%), with openEMS-class FDTD at 17/26 (65.4%).

## 0.7.0 - dispersive FDTD and coupled-multiphysics foundation

- Added Debye, Drude and Lorentz 1-D FDTD auxiliary-differential-equation material models.
- Added 3-D PMC boundary baseline, TEM-like 1-D wave-port decomposition and complex S-parameter extraction.
- Added portable legacy ASCII VTK export for 3-D electric and magnetic fields.
- Added unit/location/topology-aware multiphysics field registry.
- Added conservative exact-overlap cell remapping and conservative face-flux-to-cell transfer.
- Added partitioned fixed-point coupling with Aitken relaxation.
- Added shared-mesh DC-conduction -> Joule-heating -> transient-FEM thermal coupling.
- Added user-visible FDTD/multiphysics CLI smoke cases and sanitizer coverage.
- Reworked the high-level roadmap to use the same phase numbering as the authoritative integration tracker.
- Upgraded `scripts/integration_status.py` with explicit tracker/manifest selection and JSON output.
- Integration tracker advanced to 149/398 validated capabilities (37.4%), with openEMS-class FDTD at 14/26 (53.8%) and coupled multiphysics at 6/15 (40.0%).

## 0.6.1 - nonlinear/adaptive FEM, porous flow, magnetostatics and modal analysis

- Added a reusable Newton solver with backtracking line search and ILU(0)-GMRES Jacobian solves.
- Added a nonlinear Tri3 Poisson/reaction manufactured-solution regression using the shared Newton stack.
- Added residual/flux-jump Poisson error estimation, Dorfler marking and conforming longest-edge Tri3 refinement.
- Added saturated isotropic Darcy porous-flow FEM with analytical channel validation.
- Added 2-D magnetostatics through the out-of-plane magnetic vector potential and reconstructed magnetic flux density.
- Added a sparse generalized eigenvalue solver using M-orthogonal deflated inverse iteration and shared PCG inner solves.
- Added a fixed-free Line2 bar modal/eigenfrequency application with consistent mass and analytical mode validation.
- Added CLI/CTest smoke cases and sanitizer coverage for the new FEM capabilities.
- Integration tracker advanced to 137/398 validated capabilities (34.4%), with Elmer-class FEM at 27/46 (58.7%).

## 0.6.0 - FEM/FDTD/optics breadth and release traceability

- Completed the Phase-0 release/benchmark traceability checklist with machine-readable benchmark NDJSON, history summarization, and a tag-driven artifact workflow.
- Added a reusable FEM reference-element catalogue for Line2, Tri3, Quad4, Tet4, Hex8, Prism6 and Pyramid5 with shape functions, gradients, quadrature, isoparametric mapping and Jacobians.
- Added reusable Tri3 CSR/matrix-free Laplace assembly and mixed Dirichlet/Neumann/Robin scalar-diffusion boundaries.
- Added 2-D Poisson, heat, linear elasticity, electrostatics and DC conduction applications using shared FEM infrastructure.
- Added axisymmetric linear elasticity with cylindrical weighting and hoop strain.
- Added a convergent 3-D Tet4 Poisson baseline.
- Added a 3-D Maxwell FDTD baseline plus 1-D material preprocessing, hard/soft sources, first-order Mur boundaries, time probes and DFT monitors.
- Added sequential and paraxial optics, Sellmeier materials, Jones/Stokes polarization, Fresnel coefficients and normal-incidence multilayer transfer matrices.
- Added CLI/CTest smoke cases for 3-D FEM, Mur FDTD and sequential/paraxial optics.
- Expanded sanitizer smoke coverage across the new FEM/FDTD/optics paths.
- Integration tracker advanced to 132/397 validated capabilities (33.2%).

## 0.5.0 - tracked integration, sparse numerics and electrochemistry foundation

- Added a 386-item checkable integration tracker spanning the OpenFOAM, FluidX3D, Elmer, openEMS, Optiland and chemistry/electrochemistry target capabilities.
- Added pinned upstream provenance manifest and integration-status reporter.
- Added CSR sparse matrix assembly/SpMV, Jacobi-PCG, BiCGStab, ILU(0) and restarted GMRES.
- Replaced the default collocated momentum fixed-point sweeps with an ILU(0)-GMRES full-momentum solve while retaining the legacy fallback for A/B validation.
- Added bounded limited-linear MUSCL/Barth-Jespersen-style finite-volume reconstruction.
- Added elementary Arrhenius/mass-action reaction networks.
- Added Nernst, Butler-Volmer and Faraday electrochemistry helpers.
- Added nonlinear 1-D corrosion current/recession solver with electrolyte ohmic feedback.
- Added conservative 1-D multi-species Nernst-Planck transport with diffusion, advection and electromigration.
- Added arbitrary-`PolyMesh` Nernst-Planck transport with electroneutral and Poisson electrostatic formulations.
- Added Scharfetter-Gummel exponential-fit drift-diffusion fluxes.
- Added Butler-Volmer/Faradaic reactive electrode species-flux boundaries and anodic dissolution coupling.
- Added a multi-reaction galvanic/mixed-potential solver.
- Added implicit generic `PolyMesh` scalar advection-diffusion transport through the shared CSR/ILU-GMRES layer.
- Added chemistry/electrochemistry research and validation documentation plus CLI/CTest smoke cases.


## 0.4.2 - Phase 4C collocated pressure/velocity coupling

### Added

- patch-based fixed-value/zero-gradient/slip velocity boundary conditions;
- fixed-value/zero-gradient pressure boundary conditions;
- orthogonal/non-orthogonal face-area decomposition;
- clean-room Rhie-Chow-style face-flux reconstruction using direct pressure differences;
- arbitrary-`PolyMesh` collocated momentum predictor retaining `HbyA` and `V/aP`;
- matrix-free non-orthogonal pressure correction using shared CG;
- SIMPLE, PISO and PIMPLE-style coupling loops;
- deterministic affine-sheared Cartesian-hexa mesh generator;
- checkerboard-pressure, sheared-mesh, Poiseuille-channel and lid-cavity regressions;
- `fvm-collocated-channel`, `fvm-collocated-cavity`, and `fvm-collocated-skew` CLI cases.

### Validation

- Rhie-Chow face flux responds to an alternating collocated pressure mode;
- sheared-mesh PISO reduces continuity L2 by more than three orders of magnitude and PIMPLE improves it further;
- 32x16 pressure-driven channel agrees with analytical Poiseuille profile within the 2% regression gate (about 0.5% in the development-container run);
- 20x20 moving-lid cavity develops the expected clockwise primary recirculation and conservative face flux;
- GCC/OpenMP local CTest matrix expanded to nine tests.

## 0.4.1 - Phase 4B transient incompressible baseline

### Added

- reusable linear and first-order upwind scalar/vector face interpolation on `PolyMesh`;
- conservative face-flux scalar/vector convective-divergence operators;
- periodic staggered transient incompressible Navier-Stokes solver with conservative momentum convection, centered viscosity and shared-CG pressure projection;
- Taylor-Green analytical velocity-error measurement and refinement regression;
- `fvm-taylor-green2d` CLI case;
- sanitizer coverage for face schemes and the transient solver.

### Validation

- constant upwind transport and uniform momentum conservation regressions pass;
- linear scalar advection is exact to floating-point tolerance on the Cartesian-hexa validation mesh using linear face interpolation;
- 24x24 -> 48x48 Taylor-Green refinement reduces velocity error by the required regression factor;
- the 64x64 CLI case at `t=0.01` reports divergence near machine precision and velocity RMS error around `1.5e-7` in the development-container run.

## 0.4.0 - Phase 4A polyhedral FVM foundation

### Added

- clean-room owner/neighbour `PolyMesh` topology with boundary patches and oriented face-area vectors;
- Cartesian hexahedral mesh generator for validation;
- OpenMP cell-parallel Gauss scalar gradient and vector divergence;
- orthogonal two-point scalar Laplacian;
- reusable matrix-free conjugate-gradient core solver;
- periodic staggered finite-volume pressure projection using CG;
- manufactured FVM operator regressions and projection regressions;
- `fvm-operators3d` and `fvm-projection2d` CLI cases;
- sanitizer smoke coverage for the new mesh/operator/projection code.

### Validation

- linear Gauss gradient and linear-vector divergence are exact to floating-point tolerance on the Cartesian hexa validation mesh;
- interior quadratic Laplacian is exact to floating-point tolerance;
- the 96x80 projection CLI reduces divergence RMS from approximately `3.14e-1` to `2.85e-11` and converges its pressure solve in 7 CG iterations in the development-container case.

## 0.3.1 - Phase 3B device-resident distributed transport

### Added

- D3Q19/D3Q27 selective face/edge/corner halo plans that transmit only crossing populations.
- Shared `constexpr` host/device halo linearization.
- Persistent compact MPI point-to-point requests.
- Device-resident distributed D3Q19/D3Q27 SYCL two-grid pull reference block.
- Device halo pack/unpack kernels.
- Pinned-host SYCL USM staging transport for ordinary MPI.
- Explicit direct device-buffer `--gpu-aware-mpi` transport.
- MPI+SYCL staged distributed parity tests when both backends are enabled.
- `cfd-halo-plan` exact halo-volume reporting CLI.
- `scripts/benchmark_distributed.sh` strong/weak-scaling driver with machine-readable CSV output.
- aggregate halo message/byte counters in `cfd-distributed` benchmark reporting.
- D3Q19/D3Q27 selection and CPU/SYCL backend selection in `cfd-distributed`.
- Distributed SYCL checkpoint/restart using the Phase-3 rank-local format.

### Changed

- The CPU MPI CLI now uses the selective persistent halo exchange rather than the Phase-3A full-cell nonpersistent baseline.
- CMake project version is 0.3.1.

### Validation

- Selective virtual D3Q19/D3Q27 exchange is bitwise-identical to full-cell halo exchange on the uneven 17x15x13 2x2x2 decomposition.
- A 64^3 brick reduces exact FP32 halo payload by 74.33% for D3Q19 and 67.35% for D3Q27.
- GCC/OpenMP regressions pass.
- MPI-only and combined MPI+SYCL Phase-3B translation units pass local API syntax checks.
- Real MPI+GPU execution remains a hardware/CI validation item because this development environment has neither MPI nor AdaptiveCpp/GPU runtimes.
## 0.3.0 - Phase 3A distributed CPU/hybrid runtime

### Added

- Aspect-aware Cartesian process-grid selection and uneven brick decomposition.
- 26-neighbour face/edge/corner haloed SoA storage and virtual-rank exchange.
- Distributed D3Q19/D3Q27 two-grid pull blocks with Guo forcing.
- OpenMP collision work inside each distributed rank.
- Strict-interior / boundary-shell split for communication overlap.
- MPI Cartesian runtime, node-local rank discovery, and nonblocking `MPI_Irecv`/`MPI_Isend` halo exchange.
- Deterministic local-rank -> visible-SYCL-accelerator mapping.
- Rank-local versioned binary checkpoint/restart.
- `cfd-distributed` MPI CLI.
- Four-rank OpenMPI CI job.
- Distributed runtime and Phase-3 result documentation.

### Validation

- D3Q19 and D3Q27 2x2x2 virtual-distributed runs match the monolithic reference on an uneven 17x15x13 periodic domain.
- Face, edge, and corner ghost regions are exercised by the parity tests.
- Checkpoint reload is exact and continuation remains bitwise-identical.
- Distributed virtual halo code is covered by the ASan+UBSan smoke target.
- MPI translation units are syntax-checked locally; real MPI execution remains a CI/hardware validation item because MPI is not installed in this development container.

## 0.2.1 - Phase 2B physical LBM validation

### Added

- Guo-style body acceleration with half-step velocity correction in D2Q9, generic two-grid pull, single-grid CPU, and optional single-grid SYCL kernels.
- Halfway bounce-back stationary walls for the D2Q9 CPU reference solver.
- Moving-wall bounce-back correction and lid-driven cavity case.
- Zou-He left velocity inlet and right pressure/density outlet.
- Force-driven Poiseuille analytical regression.
- Velocity-inlet/pressure-outlet channel regression.
- Taylor-Green analytical decay/convergence regression.
- Generic templated FP32/FP64 periodic pull reference solver.
- FP32/FP64 solution-parity regression.
- `lbm-poiseuille`, `lbm-cavity`, and `lbm-channel-io` CLI cases.
- `docs/LBM_BOUNDARIES.md` and `docs/PHASE2B_RESULTS.md`.

### Validation

- D2Q9/D3Q19/D3Q27 forced in-place/two-grid parity: pass.
- Poiseuille analytical velocity profile: pass.
- Lid-driven cavity vortex-sign regression: pass.
- Zou-He inlet/outlet throughput and outlet-density checks: pass.
- Taylor-Green refinement trend: pass.
- FP32/FP64 maximum field difference gate: pass.
- Stationary-wall API reset semantics regression: pass.
- Dedicated ASan+UBSan boundary/precision smoke target: pass.
- Optional SYCL force path and forced CPU/SYCL parity tests implemented; hardware validation remains pending.

## 0.2.0 - Phase 2A fluid performance baseline

### Added

- D2Q9, D3Q19 and D3Q27 compile-time lattice descriptors with opposite-direction pairing.
- Clean-room Esoteric-Pull single-grid BGK solver for D2Q9/D3Q19/D3Q27.
- Generic two-grid pull solver for numerical/performance comparison.
- Optional AdaptiveCpp/SYCL single-grid kernels for D2Q9/D3Q19/D3Q27.
- CPU/SYCL parity tests gated by `CFD_ENABLE_SYCL`.
- D3Q19/D3Q27 descriptor-isotropy, mass, uniform-flow and Taylor-Green regression tests.
- `cfd-bench` with MLUPS, estimated DDF traffic, allocated population memory and relative mass drift.
- Side-by-side `--streaming both` benchmark mode.
- Esoteric-Pull implementation notes, GPU backend strategy and Phase-2 validation report.

### Changed

- Existing SYCL D2Q9 queue is explicitly in-order so USM kernel dependencies do not rely on implementation scheduling.
- CLI exposes D2Q9 single-grid, D3Q19 and D3Q27 examples.
- CMake project version is now 0.2.0.
