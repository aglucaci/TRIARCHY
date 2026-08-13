#!/usr/bin/env bash
set -Eeuo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly DATASETS=(
    TRI-01 TRI-02 TRI-03 TRI-05 TRI-06 TRI-07
    TRI-10 TRI-11 TRI-12 TRI-13 TRI-14 TRI-15 TRI-16
    TRI-17 TRI-18 TRI-19
    TRI-20 TRI-21 TRI-30 TRI-31
)

run_nextflow=false
if [[ "${1:-}" == "--nextflow" ]]; then
    run_nextflow=true
elif (($#)); then
    printf 'Usage: scripts/run_benchmark_suite.sh [--nextflow]\n' >&2
    exit 2
fi

cd -- "${PROJECT_DIR}"
python3 scripts/generate_benchmark_datasets.py
python3 -m unittest discover -s tests -p 'test_*.py'

if [[ "${run_nextflow}" == true ]]; then
    for dataset in "${DATASETS[@]}"; do
        printf '==> Stub-running %s\n' "${dataset}"
        ./run_triarchy.sh --stub \
            --samples "${PROJECT_DIR}/benchmarks/${dataset}/metadata/samples.tsv" \
            --controls "${PROJECT_DIR}/benchmarks/${dataset}/metadata/control_relationships.tsv" \
            --outdir "${PROJECT_DIR}/results/benchmarks/${dataset}"
    done
fi
