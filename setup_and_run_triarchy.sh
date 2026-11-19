#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly PROJECT_DIR

usage() {
    cat <<'EOF'
Usage:
  ./setup_and_run_triarchy.sh [INSTALLER_OPTIONS] -- [TRIARCHY_RUN_OPTIONS]

Install the full TRIARCHY production environment and standard references, then
run real sequencing data. Supply --samples and --controls after `--`.

Example:
  ./setup_and_run_triarchy.sh --database-dir /Volumes/TRIARCHY/refs --cpus 8 --yes -- \
    --profile laptop \
    --samples /data/study/samples.tsv \
    --controls /data/study/control_relationships.tsv \
    --outdir results/study
EOF
}

install_args=()
run_args=()
while (($#)); do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --) shift; run_args+=("$@"); break ;;
        --prefix|--environment-manager|--database-dir|--kraken-archive|--cpus)
            (($# >= 2)) || { printf 'ERROR: %s requires a value\n' "$1" >&2; exit 2; }
            install_args+=("$1" "$2"); shift 2 ;;
        --skip-nextflow|--skip-dependencies|--skip-databases|--locked-dependencies|--yes)
            install_args+=("$1"); shift ;;
        *) run_args+=("$1"); shift ;;
    esac
done

bash "${PROJECT_DIR}/install.sh" ${install_args[@]+"${install_args[@]}"}
"${PROJECT_DIR}/triarchy" run ${run_args[@]+"${run_args[@]}"}
