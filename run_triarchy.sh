#!/usr/bin/env bash
set -Eeuo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly VERSION_FILE="${PROJECT_DIR}/.nextflow-version"
readonly LOCAL_NEXTFLOW="${PROJECT_DIR}/.tools/bin/nextflow"

usage() {
    cat <<'EOF'
Usage:
  ./run_triarchy.sh --samples FILE --controls FILE [NEXTFLOW_OPTIONS]

Run the TRIARCHY production workflow on real sequencing data. All unrecognized
options are passed directly to Nextflow.

Required:
  --samples FILE       Production sample manifest
  --controls FILE      Production control-relationship table

Recommended:
  -profile laptop      Constrained production profile for a laptop
  -profile local       Default local production resources
  --outdir DIR         Production output directory
  -resume              Resume a previous production run

Environment:
  TRIARCHY_SAMPLES     Default production sample manifest
  TRIARCHY_CONTROLS    Default production control-relationship table
  TRIARCHY_OUTDIR      Default production output directory
  NXF_HOME             Nextflow cache (default: PROJECT/.nextflow)
  NXF_VER              Nextflow version (default: .nextflow-version)
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

activate_triarchy_environment() {
    command -v mamba >/dev/null 2>&1 ||
        die "Mamba is required; install Mamba and create the TRIARCHY environment"
    if [[ "${CONDA_DEFAULT_ENV:-}" != "TRIARCHY" ]]; then
        eval "$(mamba shell hook --shell bash)"
        set +u
        if ! mamba activate TRIARCHY; then
            set -u
            die "Mamba environment TRIARCHY is unavailable; create it with: mamba env create --file environment.yml"
        fi
        set -u
    fi
}

forwarded=()
while (($#)); do
    case "$1" in
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

has_option() {
    local wanted="$1"
    local argument
    for argument in ${forwarded[@]+"${forwarded[@]}"}; do
        [[ "${argument}" == "${wanted}" ]] && return 0
    done
    return 1
}

[[ -f "${VERSION_FILE}" ]] || die "missing ${VERSION_FILE}; run ./install.sh"
activate_triarchy_environment

export PATH="${PROJECT_DIR}/.tools/bin:${PATH}"
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]] &&
    command -v brew >/dev/null 2>&1; then
    export PATH="$(brew --prefix)/bin:${PATH}"
fi
export NXF_HOME="${NXF_HOME:-${PROJECT_DIR}/.nextflow}"
export NXF_VER="${NXF_VER:-$(tr -d '[:space:]' <"${VERSION_FILE}")}"

nextflow_bin="${LOCAL_NEXTFLOW}"
if [[ ! -x "${nextflow_bin}" ]]; then
    command -v nextflow >/dev/null 2>&1 || die "Nextflow is not installed; run ./install.sh first"
    nextflow_bin="$(command -v nextflow)"
fi

reference_config="${PROJECT_DIR}/conf/references.config"
if [[ -f "${reference_config}" ]] && ! has_option -c; then
    forwarded=(-c "${reference_config}" ${forwarded[@]+"${forwarded[@]}"})
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
has_option --samples || die "production run requires --samples FILE or TRIARCHY_SAMPLES"
has_option --controls || die "production run requires --controls FILE or TRIARCHY_CONTROLS"

printf 'TRIARCHY: production run with Nextflow %s\n' "${NXF_VER}"
exec "${nextflow_bin}" run "${PROJECT_DIR}/main.nf" ${forwarded[@]+"${forwarded[@]}"}
