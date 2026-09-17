# cfd_solvers 0.9.3

## Scope

v0.9.3 continues the CST-class particle/plasma expansion by adding the first self-consistent electromagnetic PIC baseline. It is intentionally scoped as a deterministic 1-D/3V periodic reference that connects field update, charge/current deposition, relativistic particle push and ionization source coupling.

## Electromagnetic PIC core

`cfd::particle::ElectromagneticPic1D` stores particles with one periodic spatial coordinate and three velocity components. The field grid contains nodal `Ex/Ey/Ez` and half-cell transverse `By/Bz` samples.

The step sequence is:

1. advance transverse magnetic fields by a half Yee step;
2. gather `E/B` at each particle position;
3. push particles with the Vay relativistic update;
4. reconstruct longitudinal current from old/new charge densities and deposit transverse CIC current from midpoint particle motion;
5. update transverse electric fields through Ampere's law;
6. finish the magnetic half step;
7. update longitudinal `Ex` from periodic Poisson when enabled;
8. optionally apply the Monte-Carlo neutral collision/ionization model.

This gives a true field-particle feedback loop rather than a post-processing particle tracker.

## Charge-conserving current baseline

`deposit_charge_conserving_current_1d` deposits old and new CIC charge density, then reconstructs the longitudinal current in Fourier space so

```text
(rho_new - rho_old) / dt + div(Jx) = 0
```

is satisfied to numerical precision on the periodic grid. The focused regression checks the returned continuity residual directly.

## Plasma source coupling

The existing deterministic Monte-Carlo neutral collision model is now callable from the EM-PIC step. Ionization events remove threshold energy and append secondary macro-particles; those particles remain in the solver state and participate in subsequent field gathers, deposition and pushes.

This is not yet a multi-species chemistry solver. It is the first validated ionization source coupling seam for later plasma/breakdown work.

## Validation status

- focused `cfd-v093-em-pic-tests`: pass;
- `particle-em-pic1d` CLI smoke case: pass;
- normal CPU CTest matrix in debug/CI-speed build: **58/58 pass**;
- focused v0.9.3 ASan+UBSan executable: **pass** with leak detection for the EM-PIC slice;
- integration tracker: **497/639 (77.8%)**;
- Phase 12: **31/51 (60.8%)**.

## Explicit limitations

The EM-PIC baseline is periodic 1-D/3V only. The next production steps are 2-D/3-D charge-conserving deposition, metallic/open/absorbing boundaries, domain decomposition, accelerator kernels, field smoothing/filtering, particle load balancing, multi-species plasma chemistry, multipactor and gas-breakdown threshold workflows.
