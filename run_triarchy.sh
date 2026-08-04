#!/usr/bin/env bash
set -Eeuo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly VERSION_FILE="${PROJECT_DIR}/.nextflow-version"
readonly LOCAL_NEXTFLOW="${PROJECT_DIR}/.tools/bin/nextflow"

usage() {
    cat <<'EOF'
Usage:
  ./run_triarchy.sh --stub [NEXTFLOW_OPTIONS]
  ./run_triarchy.sh [NEXTFLOW_OPTIONS] --samples FILE --controls FILE

Wrapper for `nextflow run main.nf`.

Wrapper options:
  --stub             Run the complete deterministic TRI-00 stub workflow
  -h, --help         Show this help
  --                 Pass all following arguments directly to Nextflow

Examples:
  ./run_triarchy.sh --stub
  ./run_triarchy.sh --stub -resume --outdir results/stub
  ./run_triarchy.sh -profile local \
    --samples /data/study/samples.tsv \
    --controls /data/study/control_relationships.tsv \
    --host_index /refs/host/grch38 \
    --kraken2_db /refs/kraken2 \
    --s16_taxonomy /refs/silva/train.fa.gz \
    --its_taxonomy /refs/unite/train.fa.gz \
    --outdir results/study

Environment:
  TRIARCHY_SAMPLES     Default production sample manifest
  TRIARCHY_CONTROLS    Default production control-relationship table
  TRIARCHY_OUTDIR      Default output directory
  NXF_HOME             Nextflow cache (default: PROJECT/.nextflow)
  NXF_VER              Nextflow version (default: .nextflow-version)

All unrecognized options are passed directly to Nextflow.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

stub_run=false
forwarded=()

while (($#)); do
    case "$1" in
        --stub)
            stub_run=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            forwarded+=("$@")
            break
            ;;
        *)
            forwarded+=("$1")
            shift
            ;;
    esac
done

[[ -f "${VERSION_FILE}" ]] || die "missing ${VERSION_FILE}; run ./install.sh"

export PATH="${PROJECT_DIR}/.tools/bin:${PATH}"
if [[ -d "${PROJECT_DIR}/.tools/env/bin" ]]; then
    export PATH="${PROJECT_DIR}/.tools/env/bin:${PATH}"
fi
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]] &&
    command -v brew >/dev/null 2>&1; then
    export PATH="$(brew --prefix)/bin:${PATH}"
fi
export NXF_HOME="${NXF_HOME:-${PROJECT_DIR}/.nextflow}"
export NXF_VER="${NXF_VER:-$(tr -d '[:space:]' <"${VERSION_FILE}")}"

if [[ ! -x "${LOCAL_NEXTFLOW}" ]] && ! command -v nextflow >/dev/null 2>&1; then
    die "Nextflow is not installed; run ./install.sh first"
fi

nextflow_bin="${LOCAL_NEXTFLOW}"
if [[ ! -x "${nextflow_bin}" ]]; then
    nextflow_bin="$(command -v nextflow)"
fi

has_option() {
    local wanted="$1"
    local argument
    for argument in ${forwarded[@]+"${forwarded[@]}"}; do
        [[ "${argument}" == "${wanted}" ]] && return 0
    done
    return 1
}

if [[ "${stub_run}" == true ]]; then
    has_option --samples ||
        forwarded+=(--samples "${PROJECT_DIR}/benchmarks/TRI-00/valid/samples.tsv")
    has_option --controls ||
        forwarded+=(--controls "${PROJECT_DIR}/benchmarks/TRI-00/valid/control_relationships.tsv")
    has_option --outdir ||
        forwarded+=(--outdir "${TRIARCHY_OUTDIR:-${PROJECT_DIR}/results/stub}")
    forwarded=(-stub-run "${forwarded[@]}")
else
    reference_config="${PROJECT_DIR}/conf/references.config"
    if [[ -f "${reference_config}" ]] && ! has_option -c; then
        forwarded=(-c "${reference_config}" "${forwarded[@]}")
    fi
    if ! has_option --samples && [[ -n "${TRIARCHY_SAMPLES:-}" ]]; then
        forwarded+=(--samples "${TRIARCHY_SAMPLES}")
    fi
    if ! has_option --controls && [[ -n "${TRIARCHY_CONTROLS:-}" ]]; then
        forwarded+=(--controls "${TRIARCHY_CONTROLS}")
    fi
    if ! has_option --outdir && [[ -n "${TRIARCHY_OUTDIR:-}" ]]; then
        forwarded+=(--outdir "${TRIARCHY_OUTDIR}")
    fi
    has_option --samples ||
        die "production mode requires --samples FILE or TRIARCHY_SAMPLES"
    has_option --controls ||
        die "production mode requires --controls FILE or TRIARCHY_CONTROLS"
fi

printf 'TRIARCHY: Nextflow %s\n' "${NXF_VER}"
printf 'TRIARCHY: output mode: %s\n' "$([[ "${stub_run}" == true ]] && printf stub || printf production)"

exec "${nextflow_bin}" run "${PROJECT_DIR}/main.nf" "${forwarded[@]}"
