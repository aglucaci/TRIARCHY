#!/usr/bin/env bash
set -Eeuo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly INSTALLER="${PROJECT_DIR}/install.sh"
readonly RUNNER="${PROJECT_DIR}/run_triarchy.sh"

usage() {
    cat <<'EOF'
Usage:
  ./setup_and_run_triarchy.sh [--demo|--lite|--stub|--production] [options] [-- NEXTFLOW_OPTIONS]

Install TRIARCHY dependencies/references, then run the workflow.

Modes:
  --demo             Install the lightweight runtime and run the TRI-31 showcase
                    deterministic workflow. This is the default when no mode
                    or production inputs are supplied.
  --lite             Install the lightweight runtime and run the laptop-friendly
                    deterministic workflow.
  --stub             Install the minimal dependency/database profile and run the
                    complete deterministic stub workflow.
  --production       Install the production dependency/database profile and run
                    on real input manifests.

Installer options:
  --prefix DIR
  --env-prefix DIR
  --environment-manager CMD
  --database-dir DIR
  --database-profile minimal|standard
  --kraken-archive FILE
  --cpus N
  --skip-nextflow
  --skip-dependencies
  --skip-databases
  --skip-tests
  --yes

Workflow options:
  --samples FILE
  --controls FILE
  --outdir DIR
  --resume           Alias for Nextflow -resume
  --                 Pass all following arguments directly to run_triarchy.sh

Examples:
  ./setup_and_run_triarchy.sh
  ./setup_and_run_triarchy.sh --lite --resume
  ./setup_and_run_triarchy.sh --stub --database-profile minimal --outdir results/stub
  ./setup_and_run_triarchy.sh --production \
    --database-dir /refs/triarchy \
    --database-profile standard \
    --cpus 8 \
    --yes \
    -- \
    -profile local \
    --samples /data/study/samples.tsv \
    --controls /data/study/control_relationships.tsv \
    --outdir results/study

The standard production database profile can exceed 400 GB and take many hours.
Completed database stages are checkpointed by install.sh and reused on reruns.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

info() {
    printf '==> %s\n' "$*"
}

has_option() {
    local wanted="$1"
    shift
    local argument
    for argument in "$@"; do
        [[ "${argument}" == "${wanted}" ]] && return 0
    done
    return 1
}

append_install_value_option() {
    (($# == 2)) || die "internal parser error for installer option"
    install_args+=("$1" "$2")
}

append_run_value_option() {
    (($# == 2)) || die "internal parser error for workflow option"
    run_args+=("$1" "$2")
}

mode=""
explicit_database_profile=false
explicit_skip_dependencies=false
explicit_skip_databases=false
explicit_skip_tests=false
explicit_yes=false
saw_production_input=false
install_args=()
run_args=()

while (($#)); do
    case "$1" in
        --demo|--lite|--stub|--production)
            [[ -z "${mode}" ]] || die "choose only one mode"
            mode="${1#--}"
            shift
            ;;
        --prefix|--env-prefix|--environment-manager|--database-dir|--kraken-archive|--cpus)
            (($# >= 2)) || die "$1 requires a value"
            append_install_value_option "$1" "$2"
            shift 2
            ;;
        --database-profile)
            (($# >= 2)) || die "$1 requires minimal or standard"
            explicit_database_profile=true
            append_install_value_option "$1" "$2"
            shift 2
            ;;
        --skip-nextflow)
            install_args+=("$1")
            shift
            ;;
        --skip-dependencies)
            explicit_skip_dependencies=true
            install_args+=("$1")
            shift
            ;;
        --skip-databases)
            explicit_skip_databases=true
            install_args+=("$1")
            shift
            ;;
        --skip-tests)
            explicit_skip_tests=true
            install_args+=("$1")
            shift
            ;;
        --yes)
            explicit_yes=true
            install_args+=("$1")
            shift
            ;;
        --samples|--controls)
            (($# >= 2)) || die "$1 requires a file"
            saw_production_input=true
            append_run_value_option "$1" "$2"
            shift 2
            ;;
        --outdir)
            (($# >= 2)) || die "$1 requires a directory"
            append_run_value_option "$1" "$2"
            shift 2
            ;;
        --resume)
            run_args+=("-resume")
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            run_args+=("$@")
            break
            ;;
        *)
            run_args+=("$1")
            shift
            ;;
    esac
done

if [[ -z "${mode}" ]]; then
    if [[ "${saw_production_input}" == true ||
        -n "${TRIARCHY_SAMPLES:-}" ||
        -n "${TRIARCHY_CONTROLS:-}" ]] ||
        has_option --samples "${run_args[@]}" ||
        has_option --controls "${run_args[@]}"; then
        mode="production"
    else
        mode="demo"
    fi
fi

[[ -f "${INSTALLER}" ]] || die "missing installer: ${INSTALLER}"
[[ -f "${RUNNER}" ]] || die "missing workflow wrapper: ${RUNNER}"

case "${mode}" in
    demo)
        [[ "${explicit_skip_dependencies}" == true ]] || install_args+=(--skip-dependencies)
        [[ "${explicit_skip_databases}" == true ]] || install_args+=(--skip-databases)
        [[ "${explicit_skip_tests}" == true ]] || install_args+=(--skip-tests)
        [[ "${explicit_yes}" == true ]] || install_args+=(--yes)
        run_args=(--demo "${run_args[@]}")
        ;;
    lite)
        [[ "${explicit_skip_dependencies}" == true ]] || install_args+=(--skip-dependencies)
        [[ "${explicit_skip_databases}" == true ]] || install_args+=(--skip-databases)
        [[ "${explicit_skip_tests}" == true ]] || install_args+=(--skip-tests)
        [[ "${explicit_yes}" == true ]] || install_args+=(--yes)
        run_args=(--lite "${run_args[@]}")
        ;;
    stub)
        [[ "${explicit_database_profile}" == true ]] || install_args+=(--database-profile minimal)
        [[ "${explicit_yes}" == true ]] || install_args+=(--yes)
        run_args=(--stub "${run_args[@]}")
        ;;
    production)
        [[ "${explicit_database_profile}" == true ]] || install_args+=(--database-profile standard)
        if ! has_option --samples "${run_args[@]}" && [[ -z "${TRIARCHY_SAMPLES:-}" ]]; then
            die "production mode requires --samples FILE or TRIARCHY_SAMPLES"
        fi
        if ! has_option --controls "${run_args[@]}" && [[ -z "${TRIARCHY_CONTROLS:-}" ]]; then
            die "production mode requires --controls FILE or TRIARCHY_CONTROLS"
        fi
        ;;
    *)
        die "unsupported mode: ${mode}"
        ;;
esac

info "Installing TRIARCHY prerequisites for ${mode} mode"
bash "${INSTALLER}" "${install_args[@]}"

info "Running TRIARCHY workflow in ${mode} mode"
bash "${RUNNER}" "${run_args[@]}"
