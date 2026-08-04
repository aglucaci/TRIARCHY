#!/usr/bin/env bash
set -uo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly DATASETS=(
    TRI-01 TRI-02 TRI-03 TRI-05 TRI-06 TRI-07
    TRI-10 TRI-11 TRI-12 TRI-13 TRI-14 TRI-15 TRI-16
    TRI-17 TRI-18 TRI-19 TRI-20 TRI-21 TRI-30
)

mode="stub"
generate=true
resume=true
results_root=""

usage() {
    cat <<'EOF'
Usage: scripts/run_all_test_datasets.sh [options]

Generate and run every executable TRIARCHY software-contract test dataset.

Options:
  --mode MODE         stub (default); production fixtures are rejected
  --no-generate       Use existing benchmark files without regenerating them
  --no-resume         Do not pass -resume to Nextflow
  --results-dir DIR   Output root (default: results/all-tests-MODE)
  -h, --help          Show this help

Examples:
  scripts/run_all_test_datasets.sh
  scripts/run_all_test_datasets.sh --no-generate

Repository fixtures exercise the complete deterministic stub graph. They are
not biologically valid production inputs. Use benchmark_resources.py with a
runtime benchmark such as TRI-T30/TRI-T60 for production resource profiling.
EOF
}

while (($#)); do
    case "$1" in
        --mode)
            (($# >= 2)) || { printf 'ERROR: --mode requires stub or production\n' >&2; exit 2; }
            mode="$2"
            shift 2
            ;;
        --no-generate)
            generate=false
            shift
            ;;
        --no-resume)
            resume=false
            shift
            ;;
        --results-dir)
            (($# >= 2)) || { printf 'ERROR: --results-dir requires a directory\n' >&2; exit 2; }
            results_root="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ "${mode}" == "stub" || "${mode}" == "production" ]] || {
    printf 'ERROR: --mode must be stub or production\n' >&2
    exit 2
}

if [[ "${mode}" == "production" ]]; then
    cat >&2 <<'EOF'
ERROR: repository test fixtures cannot be run as production sequencing data.
Their tiny deterministic read pairs validate workflow contracts but are not
guaranteed to merge, classify, or support DADA2 error learning.

Run all repository tests with:
  scripts/run_all_test_datasets.sh --mode stub --no-generate

For production profiling, install TRI-T30/TRI-T60 and use:
  python3 scripts/benchmark_resources.py --mode production --dataset TRI-T30 --dataset TRI-T60
EOF
    exit 2
fi

if [[ -z "${results_root}" ]]; then
    results_root="${PROJECT_DIR}/results/all-tests-${mode}"
elif [[ "${results_root}" != /* ]]; then
    results_root="${PROJECT_DIR}/${results_root}"
fi

mkdir -p "${results_root}/logs"
cd -- "${PROJECT_DIR}"

if [[ "${generate}" == true ]]; then
    printf '==> Generating all deterministic benchmark datasets\n'
    python3 scripts/generate_benchmark_datasets.py || exit $?
fi

printf '==> Running repository tests\n'
python3 -m unittest discover -s tests -p 'test_*.py' || exit $?

failed_count=0
failed_summary=""
for dataset in "${DATASETS[@]}"; do
    samples="${PROJECT_DIR}/benchmarks/${dataset}/metadata/samples.tsv"
    controls="${PROJECT_DIR}/benchmarks/${dataset}/metadata/control_relationships.tsv"
    outdir="${results_root}/${dataset}"
    logfile="${results_root}/logs/${dataset}.log"

    printf '==> Running %s (%s); log: %s\n' "${dataset}" "${mode}" "${logfile}"
    command=(
        "${PROJECT_DIR}/run_triarchy.sh"
        --samples "${samples}"
        --controls "${controls}"
        --outdir "${outdir}"
    )
    if [[ "${mode}" == "stub" ]]; then
        command+=(--stub)
    fi
    if [[ "${resume}" == true ]]; then
        command+=(-resume)
    fi

    if "${command[@]}" >"${logfile}" 2>&1; then
        printf '    PASS %s\n' "${dataset}"
    else
        status=$?
        failed_count=$((failed_count + 1))
        if [[ -n "${failed_summary}" ]]; then
            failed_summary="${failed_summary} "
        fi
        failed_summary="${failed_summary}${dataset}:${status}"
        printf '    FAIL %s (exit %d; continuing)\n' "${dataset}" "${status}" >&2
    fi
done

printf '\nTRIARCHY all-test summary\n'
printf '  mode:     %s\n' "${mode}"
printf '  datasets: %d\n' "${#DATASETS[@]}"
printf '  results:  %s\n' "${results_root}"

if ((failed_count)); then
    printf '  failed:   %s\n' "${failed_summary}" >&2
    printf 'Inspect per-dataset logs under %s/logs\n' "${results_root}" >&2
    exit 1
fi

printf '  failed:   0\n'
printf 'All datasets completed successfully.\n'
