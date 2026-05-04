#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

OUTPUT_DIR="result/output"
LOG_DIR="result/generation"
TOUR_DIR="result/tour"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$TOUR_DIR"

for seed in $(seq 16 30); do
  echo "=== Run ${seed}/30 (seed=${seed}) ==="
  make run-eax_tabu_combined_tabu_list ARGS="--file xpr2308.tsp --ps 200 --children 100 --trials 1 --seed ${seed} --selection greedy --eax-type EAX_5_AB --output ../../${OUTPUT_DIR}/xpr2308_greedy_tabu_s${seed}.md --log ../../${LOG_DIR}/xpr2308_greedy_tabu_s${seed}.txt --tour ../../${TOUR_DIR}/xpr2308_s${seed}.tour"
done

echo "All runs completed."
