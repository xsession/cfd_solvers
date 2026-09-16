#!/usr/bin/env bash
set -euo pipefail

MPIEXEC=${MPIEXEC:-mpiexec}
EXE=${EXE:-./build-mpi/cfd-distributed}
RANKS_LIST=${RANKS_LIST:-"1 2 4 8"}
REPEATS=${REPEATS:-3}
BACKEND=${BACKEND:-cpu}
LATTICE=${LATTICE:-d3q19}
MODE=${MODE:-strong}
WEAK_AXIS=${WEAK_AXIS:-x}
NX=${NX:-256}
NY=${NY:-256}
NZ=${NZ:-256}
STEPS=${STEPS:-200}
TAU=${TAU:-0.7}
GPU_AWARE_MPI=${GPU_AWARE_MPI:-0}

case "$MODE" in
    strong|weak) ;;
    *) echo "MODE must be strong or weak" >&2; exit 2 ;;
esac
case "$WEAK_AXIS" in
    x|y|z) ;;
    *) echo "WEAK_AXIS must be x, y, or z" >&2; exit 2 ;;
esac

extra=(--csv)
if [[ "$GPU_AWARE_MPI" == "1" ]]; then
    extra+=(--gpu-aware-mpi)
fi

printf '%s\n' 'mode,repeat,backend,lattice,ranks,pg_x,pg_y,pg_z,nx,ny,nz,steps,seconds,mlups,halo_messages_per_step,halo_bytes_per_step,halo_gib,mass'
for ranks in $RANKS_LIST; do
    run_nx=$NX
    run_ny=$NY
    run_nz=$NZ
    if [[ "$MODE" == "weak" ]]; then
        case "$WEAK_AXIS" in
            x) run_nx=$((NX * ranks)) ;;
            y) run_ny=$((NY * ranks)) ;;
            z) run_nz=$((NZ * ranks)) ;;
        esac
    fi

    for ((repeat=1; repeat<=REPEATS; ++repeat)); do
        line=$($MPIEXEC -n "$ranks" "$EXE" \
            --backend "$BACKEND" --lattice "$LATTICE" \
            --nx "$run_nx" --ny "$run_ny" --nz "$run_nz" --steps "$STEPS" --tau "$TAU" \
            "${extra[@]}" | tail -n 1)
        printf '%s,%s,%s\n' "$MODE" "$repeat" "$line"
    done
done
