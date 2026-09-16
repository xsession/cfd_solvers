# Phase 2 fluid baseline results

## Correctness result

The single-grid Esoteric-Pull implementation is checked against a separately implemented conventional two-grid pull solver for D2Q9, D3Q19 and D3Q27. Tests exercise both odd and even Esoteric-Pull parity layouts after Taylor-Green evolution and compare density plus all velocity components.

Additional tests cover descriptor isotropy, uniform-flow invariance, mass conservation and viscous Taylor-Green energy decay. See `VALIDATION.md` for the complete build/test matrix.

## Memory result

Population storage is exactly halved relative to the conventional two-grid reference:

```text
single grid: Q * N * sizeof(float)
two grids : 2 * Q * N * sizeof(float)
```

For the benchmark cases used in this phase:

- D3Q19 `64^3`: 19.0 MiB instead of 38.0 MiB;
- D3Q27 `48^3`: 11.39 MiB instead of 22.78 MiB.

This is the main Phase-2 performance result because it is deterministic and hardware-independent.

## Representative CPU throughput

On the development container (5 vCPUs from an Intel Xeon Platinum 8370C VM, GCC 14.2 Release/OpenMP), a longer 200-step snapshot gave:

| Lattice | Streaming | MLUPS | Estimated DDF GB/s |
|---|---|---:|---:|
| D3Q19 `64^3` | in-place | 18.03 | 2.74 |
| D3Q19 `64^3` | two-grid pull | 17.98 | 2.73 |
| D3Q27 `48^3` | in-place | 12.03 | 2.60 |
| D3Q27 `48^3` | two-grid pull | 11.68 | 2.52 |

These are regression/smoke measurements from a virtualized shared environment, not portable hardware claims. The important observation is that halving DDF allocation did not impose a material throughput penalty in this sample.

## SYCL result

D2Q9/D3Q19/D3Q27 SYCL kernels and CPU/SYCL parity tests are present, but this development container does not contain AdaptiveCpp or a GPU runtime. Accelerator execution is therefore not marked release-validated yet.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFD_ENABLE_OPENMP=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure

./build/cfd-bench --lattice d3q19 --streaming both \
  --nx 64 --ny 64 --nz 64 --warmup 20 --steps 200 --threads 5
```
