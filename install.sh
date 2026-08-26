#!/usr/bin/env bash
set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly VERSION_FILE="${SCRIPT_DIR}/.nextflow-version"
readonly DATABASE_DIR_FILE="${SCRIPT_DIR}/references/.database-dir"

PREFIX="${SCRIPT_DIR}/.tools/bin"
readonly ENV_NAME="TRIARCHY"
ENV_PREFIX=""
DATABASE_DIR="${SCRIPT_DIR}/references/data"
DATABASE_DIR_EXPLICIT=false
KRAKEN_ARCHIVE=""
INSTALL_CPUS=4
ENVIRONMENT_MANAGER="${TRIARCHY_ENV_MANAGER:-mamba}"
SKIP_NEXTFLOW=false
SKIP_DEPENDENCIES=false
LOCKED_DEPENDENCIES=false
SKIP_DATABASES=false
ASSUME_YES=false
SETUP_STARTED="${SECONDS}"

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    COLOR_CYAN=$'\033[1;36m'
    COLOR_GREEN=$'\033[1;32m'
    COLOR_YELLOW=$'\033[1;33m'
    COLOR_RED=$'\033[1;31m'
    COLOR_RESET=$'\033[0m'
else
    COLOR_CYAN=""
    COLOR_GREEN=""
    COLOR_YELLOW=""
    COLOR_RED=""
    COLOR_RESET=""
fi

usage() {
    cat <<'EOF'
Usage: ./install.sh [options]

Set up a complete local TRIARCHY environment.

Options:
  --prefix DIR          Install Nextflow in DIR (default: .tools/bin)
  --environment-manager CMD
                        mamba, micromamba, conda, or auto (default: mamba)
  --database-dir DIR    Reference destination (default: last successful location,
                        otherwise references/data)
  --kraken-archive FILE Use an existing PlusPFP .tar.gz instead of downloading
  --cpus N              Threads for the locally built host index (default: 4)
  --skip-nextflow       Do not download or validate Nextflow
  --skip-dependencies   Do not create/update the Mamba environment
  --locked-dependencies Install dependencies from environment.lock.txt instead
                        of solving environment.yml
  --skip-databases      Do not download or build references
  --yes                 Accept the large standard database installation
  -h, --help            Show this help

Environment:
  NXF_VER            Override the version pinned in .nextflow-version
  NXF_HOME           Override Nextflow's cache directory
  TRIARCHY_ENV_MANAGER
                     Override the environment manager (default: mamba)
  NO_COLOR           Disable ANSI colors

The complete production reference set can exceed 400 GB. Downloads are
resumable and validated before the generated configuration is activated.
EOF
}

die() {
    printf '%bERROR:%b %s\n' "${COLOR_RED}" "${COLOR_RESET}" "$*" >&2
    exit 1
}

info() {
    printf '%b==>%b %s\n' "${COLOR_CYAN}" "${COLOR_RESET}" "$*"
}

success() {
    printf '%b✓%b %s\n' "${COLOR_GREEN}" "${COLOR_RESET}" "$*"
}

warning() {
    printf '%bWARNING:%b %s\n' "${COLOR_YELLOW}" "${COLOR_RESET}" "$*" >&2
}

format_duration() {
    local seconds="$1"
    printf '%02d:%02d:%02d' \
        $((seconds / 3600)) $(((seconds % 3600) / 60)) $((seconds % 60))
}

while (($#)); do
    case "$1" in
        --prefix)
            (($# >= 2)) || die "--prefix requires a directory"
            PREFIX="$2"
            shift 2
            ;;
        --environment-manager)
            (($# >= 2)) || die "--environment-manager requires a command"
            ENVIRONMENT_MANAGER="$2"
            shift 2
            ;;
        --database-dir)
            (($# >= 2)) || die "--database-dir requires a directory"
            DATABASE_DIR="$2"
            DATABASE_DIR_EXPLICIT=true
            shift 2
            ;;
        --kraken-archive)
            (($# >= 2)) || die "--kraken-archive requires a file"
            KRAKEN_ARCHIVE="$2"
            shift 2
            ;;
        --cpus)
            (($# >= 2)) || die "--cpus requires an integer"
            INSTALL_CPUS="$2"
            shift 2
            ;;
        --skip-nextflow)
            SKIP_NEXTFLOW=true
            shift
            ;;
        --skip-dependencies)
            SKIP_DEPENDENCIES=true
            shift
            ;;
        --locked-dependencies)
            LOCKED_DEPENDENCIES=true
            shift
            ;;
        --skip-databases)
            SKIP_DATABASES=true
            shift
            ;;
        --yes)
            ASSUME_YES=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown option: $1"
            ;;
    esac
done

# Reuse the last successful database destination. Older installations did not
# write the state file, so recover the data root from references.config once.
if [[ "${DATABASE_DIR_EXPLICIT}" == false ]]; then
    if [[ -s "${DATABASE_DIR_FILE}" ]]; then
        DATABASE_DIR="$(head -n 1 "${DATABASE_DIR_FILE}")"
    elif [[ -s "${SCRIPT_DIR}/conf/references.config" ]]; then
        configured_kraken_dir="$(
            awk -F "'" '/^[[:space:]]*kraken2_db[[:space:]]*=/ { print $2; exit }' \
                "${SCRIPT_DIR}/conf/references.config"
        )"
        if [[ -n "${configured_kraken_dir}" ]]; then
            DATABASE_DIR="$(dirname -- "${configured_kraken_dir}")"
        fi
    fi
fi

if [[ "${SKIP_DATABASES}" == false && "${DATABASE_DIR}" == /Volumes/* ]]; then
    database_volume="${DATABASE_DIR#/Volumes/}"
    database_volume="${database_volume%%/*}"
    [[ -d "/Volumes/${database_volume}" ]] ||
        die "reference volume /Volumes/${database_volume} is not mounted"
fi

[[ "${INSTALL_CPUS}" =~ ^[1-9][0-9]*$ ]] || die "--cpus must be a positive integer"
case "${ENVIRONMENT_MANAGER}" in
    mamba|micromamba|conda|auto) ;;
    *) die "--environment-manager must be mamba, micromamba, conda, or auto" ;;
esac

command -v bash >/dev/null 2>&1 || die "bash is required"
command -v java >/dev/null 2>&1 || die "Java 17 or later is required"
command -v python3 >/dev/null 2>&1 || die "Python 3 is required for repository checks"

java_version="$(
    java -version 2>&1 |
        awk -F '"' '/version/ { print $2; exit }'
)"
[[ -n "${java_version}" ]] || die "could not determine the Java version"
java_major="${java_version%%.*}"
if [[ "${java_major}" == "1" ]]; then
    java_major="$(cut -d. -f2 <<<"${java_version}")"
fi
[[ "${java_major}" =~ ^[0-9]+$ ]] || die "unrecognized Java version: ${java_version}"
((java_major >= 17)) || die "Java 17+ is required; found ${java_version}"
info "Java ${java_version}"

mkdir -p \
    "${SCRIPT_DIR}/conf" \
    "${SCRIPT_DIR}/modules" \
    "${SCRIPT_DIR}/results" \
    "${PREFIX}"

if [[ "${SKIP_NEXTFLOW}" == false ]]; then
    command -v curl >/dev/null 2>&1 || die "curl is required to install Nextflow"
    command -v install >/dev/null 2>&1 || die "the POSIX install utility is required"
    [[ -f "${VERSION_FILE}" ]] || die "missing ${VERSION_FILE}"

    nextflow_version="${NXF_VER:-$(tr -d '[:space:]' <"${VERSION_FILE}")}"
    [[ "${nextflow_version}" =~ ^[0-9]{2}\.[0-9]{2}\.[0-9]+([.-][A-Za-z0-9]+)?$ ]] ||
        die "invalid Nextflow version: ${nextflow_version}"

    nextflow_bin="${PREFIX}/nextflow"
    needs_install=true
    if [[ -x "${nextflow_bin}" ]]; then
        installed_version="$(
            NXF_VER="${nextflow_version}" \
                NXF_HOME="${NXF_HOME:-${SCRIPT_DIR}/.nextflow}" \
                "${nextflow_bin}" -version 2>/dev/null |
                awk '/version/ { print $2; exit }'
        )"
        if [[ "${installed_version}" == "${nextflow_version}" ]]; then
            needs_install=false
            info "Nextflow ${installed_version} is already installed"
        fi
    fi

    if [[ "${needs_install}" == true ]]; then
        tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/triarchy-nextflow.XXXXXX")"
        cleanup() {
            rm -rf -- "${tmp_dir}"
        }
        trap cleanup EXIT

        info "Installing Nextflow ${nextflow_version} from get.nextflow.io"
        curl --fail --silent --show-error --location \
            https://get.nextflow.io \
            --output "${tmp_dir}/nextflow"
        (
            cd -- "${tmp_dir}"
            NXF_VER="${nextflow_version}" \
                NXF_HOME="${NXF_HOME:-${SCRIPT_DIR}/.nextflow}" \
                bash "${tmp_dir}/nextflow"
        )
        install -m 0755 "${tmp_dir}/nextflow" "${nextflow_bin}"
        installed_version="$(
            NXF_VER="${nextflow_version}" \
                NXF_HOME="${NXF_HOME:-${SCRIPT_DIR}/.nextflow}" \
                "${nextflow_bin}" -version |
                awk '/version/ { print $2; exit }'
        )"
        [[ "${installed_version}" == "${nextflow_version}" ]] ||
            die "installed Nextflow ${installed_version}; expected ${nextflow_version}"
        info "Installed Nextflow ${installed_version} at ${nextflow_bin}"
    fi
fi

if [[ "${ENVIRONMENT_MANAGER}" == "auto" ]]; then
    if command -v mamba >/dev/null 2>&1; then
        environment_manager="mamba"
    elif command -v micromamba >/dev/null 2>&1; then
        environment_manager="micromamba"
    elif command -v conda >/dev/null 2>&1; then
        environment_manager="conda"
    else
        die "Mamba, Micromamba, or Conda is required"
    fi
else
    environment_manager="${ENVIRONMENT_MANAGER}"
    command -v "${environment_manager}" >/dev/null 2>&1 ||
        die "${environment_manager} is required; install it or use --environment-manager auto"
fi

environment_file="${SCRIPT_DIR}/environment.yml"
environment_lock="${SCRIPT_DIR}/environment.lock.txt"
if [[ "${SKIP_DEPENDENCIES}" == false ]]; then
    info "Installing the named ${ENV_NAME} production environment with ${environment_manager}"
    if [[ "${LOCKED_DEPENDENCIES}" == true ]]; then
        [[ -f "${environment_lock}" ]] ||
            die "missing ${environment_lock}; run without --locked-dependencies first"
        dependency_args=(install --name "${ENV_NAME}" --file "${environment_lock}")
    else
        [[ -f "${environment_file}" ]] || die "missing ${environment_file}"
        dependency_args=(env update --name "${ENV_NAME}" --file "${environment_file}" --prune)
    fi
    [[ "${ASSUME_YES}" == true ]] && dependency_args+=(--yes)
    "${environment_manager}" "${dependency_args[@]}"
fi

ENV_PREFIX="$(
    "${environment_manager}" env list --json |
        python3 -c 'import json,sys; envs=json.load(sys.stdin)["envs"]; matches=[p for p in envs if p.rstrip("/").endswith("/TRIARCHY")]; print(matches[0] if matches else "")'
)"
[[ -x "${ENV_PREFIX}/bin/python3" ]] ||
    die "named environment ${ENV_NAME} is missing; run: mamba env create --file environment.yml"
ln -sfn "${SCRIPT_DIR}/triarchy" "${ENV_PREFIX}/bin/triarchy"
"${environment_manager}" list --name "${ENV_NAME}" --explicit >"${environment_lock}"

if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    native_args=(--env-prefix "${ENV_PREFIX}")
    if [[ "${SKIP_DEPENDENCIES}" == false ]]; then
        native_args+=(--install)
    fi
    info "Configuring native Apple-silicon Bowtie2, Kraken2, Bracken, and DIAMOND"
    "${SCRIPT_DIR}/scripts/configure_conda_native_tools.sh" \
        "${native_args[@]}"
fi

if [[ "${SKIP_DATABASES}" == false ]]; then
    database_args=(
        --database-dir "${DATABASE_DIR}"
        --env-prefix "${ENV_PREFIX}"
        --cpus "${INSTALL_CPUS}"
    )
    if [[ "${ASSUME_YES}" == true ]]; then
        database_args+=(--yes)
    fi
    if [[ -n "${KRAKEN_ARCHIVE}" ]]; then
        database_args+=(--kraken-archive "${KRAKEN_ARCHIVE}")
    fi
    "${SCRIPT_DIR}/scripts/install_databases.sh" "${database_args[@]}"
    mkdir -p "${SCRIPT_DIR}/references"
    printf '%s\n' "${DATABASE_DIR}" >"${DATABASE_DIR_FILE}"
fi

if [[ -f "${SCRIPT_DIR}/main.nf" && -f "${SCRIPT_DIR}/nextflow.config" ]]; then
    info "Nextflow workflow entrypoint detected"
else
    printf '\n'
    warning "Setup completed, but main.nf and nextflow.config are not both present."
    printf 'The repository checks are usable; the assay workflow is not yet runnable.\n' >&2
fi

success "TRIARCHY setup complete in $(format_duration $((SECONDS - SETUP_STARTED)))"
cat <<EOF

Add the local Nextflow installation to this shell:
  mamba activate "${ENV_NAME}"
  export PATH="${PREFIX}:\$PATH"
  export NXF_HOME="${NXF_HOME:-${SCRIPT_DIR}/.nextflow}"
  export NXF_VER="${nextflow_version:-$(tr -d '[:space:]' <"${VERSION_FILE}")}"

Then inspect the installation with:
  nextflow -version
EOF
