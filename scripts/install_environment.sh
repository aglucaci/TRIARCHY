#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly SCRIPT_DIR
readonly VERSION_FILE="${SCRIPT_DIR}/.nextflow-version"

PREFIX="${SCRIPT_DIR}/.tools/bin"
readonly ENV_NAME="TRIARCHY"
ENV_PREFIX=""
ENVIRONMENT_MANAGER="${TRIARCHY_ENV_MANAGER:-mamba}"
SKIP_NEXTFLOW=false
SKIP_DEPENDENCIES=false
LOCKED_DEPENDENCIES=false
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
Usage: scripts/install_environment.sh [options]

Install the TRIARCHY software: the named TRIARCHY conda environment, the pinned
Nextflow release, and the `triarchy` command on the environment's PATH.
Reference databases are installed separately with `triarchy databases install`.

Options:
  --prefix DIR          Install Nextflow in DIR (default: .tools/bin)
  --environment-manager CMD
                        mamba, micromamba, conda, or auto (default: mamba)
  --skip-nextflow       Do not download or validate Nextflow
  --skip-dependencies   Do not create/update the Mamba environment
  --locked-dependencies Install the committed lock for this platform
                        (locks/environment.<platform>.lock.txt) instead of
                        solving environment.yml
  --yes                 Do not prompt during environment creation
  -h, --help            Show this help

Environment:
  NXF_VER            Override the version pinned in .nextflow-version
  NXF_HOME           Override Nextflow's cache directory
  TRIARCHY_ENV_MANAGER
                     Override the environment manager (default: mamba)
  NO_COLOR           Disable ANSI colors
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
case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) conda_platform="linux-64" ;;
    Linux-aarch64) conda_platform="linux-aarch64" ;;
    Darwin-arm64) conda_platform="osx-arm64" ;;
    Darwin-x86_64) conda_platform="osx-64" ;;
    *) conda_platform="unknown" ;;
esac
committed_lock="${SCRIPT_DIR}/locks/environment.${conda_platform}.lock.txt"
if [[ "${SKIP_DEPENDENCIES}" == false ]]; then
    info "Installing the named ${ENV_NAME} production environment with ${environment_manager}"
    if [[ "${LOCKED_DEPENDENCIES}" == true ]]; then
        [[ -f "${committed_lock}" ]] ||
            die "no committed lock for ${conda_platform} (${committed_lock}); run without --locked-dependencies"
        if "${environment_manager}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
            dependency_args=(install --name "${ENV_NAME}" --file "${committed_lock}")
        else
            dependency_args=(create --name "${ENV_NAME}" --file "${committed_lock}")
        fi
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
for required_tool in cutadapt itsxpress vsearch Rscript; do
    [[ -x "${ENV_PREFIX}/bin/${required_tool}" ]] ||
        die "${required_tool} is missing from ${ENV_NAME}; rerun the installer without --skip-dependencies"
done
ln -sfn "${SCRIPT_DIR}/triarchy" "${ENV_PREFIX}/bin/triarchy"
# Record the exact packages installed on this machine (local, not committed).
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

if [[ -f "${SCRIPT_DIR}/main.nf" && -f "${SCRIPT_DIR}/nextflow.config" ]]; then
    info "Nextflow workflow entrypoint detected"
else
    printf '\n'
    warning "Setup completed, but main.nf and nextflow.config are not both present."
    printf 'The repository checks are usable; the assay workflow is not yet runnable.\n' >&2
fi

success "TRIARCHY setup complete in $(format_duration $((SECONDS - SETUP_STARTED)))"
cat <<EOF

Next, install the reference databases (several hundred GB; resumable):
  triarchy databases install --destination /path/to/references --yes

Add the local Nextflow installation to this shell:
  mamba activate "${ENV_NAME}"
  export PATH="${PREFIX}:\$PATH"
  export NXF_HOME="${NXF_HOME:-${SCRIPT_DIR}/.nextflow}"
  export NXF_VER="${nextflow_version:-$(tr -d '[:space:]' <"${VERSION_FILE}")}"

Then inspect the installation with:
  nextflow -version
EOF
