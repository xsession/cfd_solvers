#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${1:-build}"
THREADS="${CFD_BENCH_THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"
STEPS="${CFD_BENCH_STEPS:-100}"
WARMUP="${CFD_BENCH_WARMUP:-10}"

bench() {
    "$BUILD_DIR/cfd-bench" --lattice "$1" --nx "$2" --ny "$3" --nz "$4" \
        --steps "$STEPS" --warmup "$WARMUP" --threads "$THREADS" --streaming both --csv
}

echo "# CPU LBM benchmark: threads=$THREADS steps=$STEPS warmup=$WARMUP" >&2
bench d2q9 512 512 1
bench d3q19 64 64 64
bench d3q27 56 56 56
