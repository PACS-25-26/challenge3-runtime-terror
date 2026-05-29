#!/usr/bin/env bash

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
ROOT_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
DATA_DIR="$ROOT_DIR/test/data"
OUTPUT_CSV="$DATA_DIR/scalability_results.csv"
HW_INFO="$DATA_DIR/hw.info"

GRID_SIZES=${GRID_SIZES:-"128 256 512"}
MPI_PROCESSES=${MPI_PROCESSES:-"1 2 4"}
THREADS=${THREADS:-2}

extract_metric() {
    local output="$1"
    local pattern="$2"
    local column="$3"
    printf '%s\n' "$output" | grep "$pattern" | awk -v col="$column" '{print $col}' || echo "NaN"
}

run_test() {
    local processes=$1
    local n=$2

    echo "Running case: n=$n | MPI_Ranks=$processes | OMP_Threads=$THREADS"
    
    local output
    output=$(OMP_NUM_THREADS="$THREADS" mpirun --oversubscribe -np "$processes" "$ROOT_DIR/main.out" "$n" 2>&1)

    local serial_time
    local parallel_time
    local speedup
    local error

    serial_time=$(extract_metric "$output" "\[SERIAL\] Time:" 3)
    parallel_time=$(extract_metric "$output" "\[PARALLEL\] Time:" 3)
    speedup=$(extract_metric "$output" "Speedup:" 2)
    error=$(extract_metric "$output" "\[PARALLEL\] Time:" 7) 

    printf '%s,%s,%s,%s,%s,%s,%s\n' \
        "$n" "$processes" "$THREADS" "$serial_time" "$parallel_time" "$speedup" "$error" \
        >> "$OUTPUT_CSV"
}

# MAIN

mkdir -p "$DATA_DIR"

echo "Saving hardware info..."
cat /proc/cpuinfo > "$HW_INFO"

echo "Building project..."
make -C "$ROOT_DIR" clean all

echo "Initializing CSV..."
printf 'n,mpi_ranks,omp_threads,serial_time_s,parallel_time_s,speedup,l2_error\n' > "$OUTPUT_CSV"


for n in $GRID_SIZES; do
    for processes in $MPI_PROCESSES; do
        run_test "$processes" "$n"
    done
done

echo "Benchmarking complete. Results saved in $OUTPUT_CSV"