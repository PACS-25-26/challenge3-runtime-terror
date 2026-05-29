#!/usr/bin/env bash

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
ROOT_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
DATA_DIR="$ROOT_DIR/test/data"
OUTPUT_CSV="$DATA_DIR/scalability_results.csv"
HW_INFO="$DATA_DIR/hw.info"

GRID_SIZES=${GRID_SIZES:-"16 32 64 128 256"}
MPI_PROCESSES=${MPI_PROCESSES:-"1 2 4"}
THREADS=${THREADS:-2}

run_test() {
    local processes=$1
    local n=$2

    echo "Running case: n=$n | MPI_Ranks=$processes | OMP_Threads=$THREADS"
    
    local output
    output=$(OMP_NUM_THREADS="$THREADS" mpirun --oversubscribe -np "$processes" "$ROOT_DIR/main.out" "$n" 2>&1)

    local ser_time=$(printf '%s\n' "$output" | grep "\[SERIAL\] Time:" | awk '{print $3}')
    local ser_error=$(printf '%s\n' "$output" | grep "\[SERIAL\] Time:" | awk '{print $7}')

    local par_time=$(printf '%s\n' "$output" | grep "\[PARALLEL\] Time:" | awk '{print $3}')
    local par_error=$(printf '%s\n' "$output" | grep "\[PARALLEL\] Time:" | awk '{print $7}')
    local par_speedup=$(printf '%s\n' "$output" | grep "Speedup:" | sed -n '1p' | awk '{print $2}')

    local sch_time=$(printf '%s\n' "$output" | grep "\[SCHWARZ\] Time:" | awk '{print $3}')
    local sch_error=$(printf '%s\n' "$output" | grep "\[SCHWARZ\] Time:" | awk '{print $7}')
    local sch_speedup=$(printf '%s\n' "$output" | grep "Speedup:" | sed -n '2p' | awk '{print $2}')

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$n" "$processes" "$THREADS" \
        "$ser_time" "$par_time" "$par_speedup" "$par_error" \
        "$sch_time" "$sch_speedup" "$sch_error" \
        >> "$OUTPUT_CSV"
}

# MAIN

mkdir -p "$DATA_DIR"

echo "Saving hardware info..."
cat /proc/cpuinfo > "$HW_INFO"

echo "Building project..."
make -C "$ROOT_DIR" clean all

echo "Initializing CSV..."
printf 'n,mpi_ranks,omp_threads,serial_time_s,parallel_time_s,speedup,l2_error,schwarz_time_s,schwarz_speedup,schwarz_error\n' > "$OUTPUT_CSV"


for n in $GRID_SIZES; do
    for processes in $MPI_PROCESSES; do
        run_test "$processes" "$n"
    done
done

echo "Benchmarking complete. Results saved in $OUTPUT_CSV"