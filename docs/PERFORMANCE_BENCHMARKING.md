# Performance benchmarking

The repository now has a repeatable CPU benchmark path for deciding whether a
solver change is worth keeping. It measures correctness and performance in the
same run, records build and host provenance, stores repeated samples, and can
gate a change against a saved baseline.

## What is measured

`cfd-performance-bench` currently compares:

- a regular-grid seven-point stencil against an equivalent CSR matrix;
- an explicitly separated stencil application plus residual pass against a
  fused application/residual pass;
- the previous distributed CSR halo lookup path (`lower_bound` per off-rank
  entry) against the remapped, parallel path now in
  `DistributedCsrPartition::multiply()`.

Every kernel emits:

- elapsed time and median-ready `ms_per_step`;
- cells/second, algorithmic estimated GB/s, and estimated FLOP/s;
- dimensions, nonzeros, effective thread count, warmup and measured steps;
- a checksum to detect optimized-away or invalid work;
- an independent oracle error where two implementations are compared;
- residual norm where the benchmark produces a residual;
- a Linux resident-memory snapshot (`working_set_bytes`).

The traffic and FLOP rates are model-based estimates, not hardware-counter
measurements. They are useful for comparing implementations with the same
operation model; they should not be presented as measured DRAM bandwidth.
Resident memory is a process snapshot rather than a peak allocation trace.

## One benchmark run

Build the benchmark with the same flags used for the solver under evaluation:

```bash
cmake --build build-review --target cfd-performance-bench -j2
./build-review/cfd-performance-bench \
  --nx 64 --ny 64 --nz 64 \
  --warmup 5 --steps 30 --threads 4
```

The machine-readable lines begin with `CFD_PERF` and
`CFD_PERF_SETUP`. Human-readable compiler or shell output can therefore be
filtered without changing the benchmark schema.

## Reproducible matrix and history

The Python runner executes a size/thread/repetition matrix and writes one JSON
record containing raw samples, medians, spread, setup costs, git state, host
affinity, compiler, and CMake feature flags:

```bash
python3 scripts/run_performance_benchmarks.py \
  --build-dir build-review \
  --output build-review/performance/current.json \
  --history build-review/performance/history.jsonl \
  --size 32x32x32 --size 64x64x64 \
  --threads 1 --threads 2 --threads 4 \
  --warmup 5 --steps 30 --repeats 5
```

For a fast smoke matrix:

```bash
python3 scripts/run_performance_benchmarks.py \
  --build-dir build-review \
  --output build-review/performance/quick.json \
  --quick --threads 1 --threads 2 --repeats 3
```

The output schema is `cfd_solvers.performance_suite.v1`. Keep result files
with the commit that produced them; the runner records dirty-state status but
does not assume that a dirty tree is invalid.

## Baseline and regression gate

Save a trusted result as a baseline, then compare a later build using the same
matrix:

```bash
python3 scripts/run_performance_benchmarks.py \
  --build-dir build-review \
  --output build-review/performance/candidate.json \
  --baseline build-review/performance/baseline.json \
  --check --max-regression-percent 5 \
  --size 32x32x32 --size 64x64x64 \
  --threads 1 --threads 2 --threads 4 \
  --warmup 5 --steps 30 --repeats 5
```

The gate compares median elapsed time for matching benchmark, implementation,
dimensions, and effective thread count. A positive `change_percent` means the
candidate is slower; `speedup > 1` means it is faster. A missing key is marked
`new` rather than silently treated as a pass.

Use the following policy when deciding whether to keep an optimization:

1. First require all oracle errors to remain at numerical tolerance and all
   correctness tests to pass.
2. Compare medians, not the fastest single sample. Inspect stdev and the raw
   samples when a change is close to the threshold.
3. Use at least one cache-sized case and one case large enough to expose
   memory traffic. Repeat on a pinned or otherwise stable CPU allocation when
   possible.
4. Treat setup time separately from steady-state time. A structured operator
   may be cheaper to construct while a CSR path may be reusable across many
   solves; the right choice depends on reuse count.
5. Do not claim a general speedup from one size, one thread count, or a quick
   run. Record the JSON result and report the matrix that supports the claim.

The latest representative run in this workspace used 32³ and 64³ grids,
1/2/4 threads, five warmups, 30 measured steps, and three repeats. The
following ratios are median elapsed time for the first implementation divided
by the second implementation, so values above 1 are faster:

| Grid | Threads | Structured stencil / CSR | Remapped distributed / legacy | Fused residual / explicit |
|---|---:|---:|---:|---:|
| 32³ | 1 | 2.49x | 1.39x | 0.99x |
| 32³ | 2 | 2.92x | 2.24x | 1.26x |
| 32³ | 4 | 4.08x | 3.07x | 1.31x |
| 64³ | 1 | 4.61x | 1.23x | 1.11x |
| 64³ | 2 | 4.88x | 1.96x | 0.98x |
| 64³ | 4 | 5.26x | 2.57x | 2.42x |

These figures are measurements on the recorded host/build, not portable
guarantees. The fused residual result is intentionally shown as conditional:
it is not a universal win at every matrix/thread point, which is exactly the
kind of trade-off this benchmark should expose.

## Adding a new optimization

Keep the benchmark workload and oracle fixed, add a new implementation label,
and preserve the existing implementation labels as controls. Then:

1. run the correctness test and baseline matrix before the change;
2. implement the optimization without changing the mathematical workload;
3. run the same matrix with the same build flags and affinity;
4. compare median time, spread, memory, setup cost, and oracle error;
5. retain the change only when the improvement is repeatable and its memory or
   setup trade-off is acceptable for the intended solver use.

This makes results auditable across the structured-stencil work inspired by the
research repositories and the existing unstructured/CSR solvers, without
pretending that a regular-grid kernel replaces a general sparse operator.
