#!/usr/bin/env bash
set -Eeuo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"

DATABASE_DIR="${PROJECT_DIR}/references/data"
ENV_PREFIX="${PROJECT_DIR}/.tools/env"
CPUS="${TRIARCHY_INSTALL_CPUS:-4}"
readonly PROFILE="standard"
KRAKEN_ARCHIVE=""
SILVA_ARCHIVE=""
ASSUME_YES=false
CURRENT_STEP=0
TOTAL_STEPS=0
INSTALL_STARTED="${SECONDS}"
ACTIVE_STEP_KEY=""
ACTIVE_STEP_LABEL=""
ACTIVE_STEP_STARTED=0

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    COLOR_BLUE=$'\033[1;34m'
    COLOR_CYAN=$'\033[1;36m'
    COLOR_GREEN=$'\033[1;32m'
    COLOR_YELLOW=$'\033[1;33m'
    COLOR_RED=$'\033[1;31m'
    COLOR_DIM=$'\033[2m'
    COLOR_RESET=$'\033[0m'
else
    COLOR_BLUE=""
    COLOR_CYAN=""
    COLOR_GREEN=""
    COLOR_YELLOW=""
    COLOR_RED=""
    COLOR_DIM=""
    COLOR_RESET=""
fi

usage() {
    cat <<'EOF'
Usage: scripts/install_databases.sh [options]

Download and install TRIARCHY production reference databases.

Options:
  --database-dir DIR   Destination (default: references/data)
  --env-prefix DIR     Dependency environment (default: .tools/env)
  --cpus N             Build threads (default: 4)
  --silva-archive FILE Use an existing SILVA 138.2 DADA2 .fa.gz archive
  --kraken-archive FILE
                       Use an existing PlusPFP .tar.gz instead of downloading
  --yes                Accept the large download without prompting
  -h, --help           Show this help

The installer provisions the full host, SILVA, ITS, MetaPhlAn,
HUMAnN ChocoPhlAn/UniRef90, and Kraken2 PlusPFP/Bracken production references.

Overrides:
  TRIARCHY_SILVA_URL   Additional checksum-identical SILVA mirror
  TRIARCHY_ITS_URL     DADA2-compatible ITS taxonomy archive
  TRIARCHY_HUMAN_ACCESSION  NCBI assembly accession
  TRIARCHY_KRAKEN_URL  Prebuilt Kraken2 archive URL
  TRIARCHY_KRAKEN_MD5_URL  Prebuilt Kraken2 checksum-manifest URL
  TRIARCHY_KRAKEN_SIZE  Expected archive size in bytes
  TRIARCHY_DOWNLOAD_ATTEMPTS  Retry attempts per download stage (default: 20)
  TRIARCHY_MIRROR_ATTEMPTS  Retry attempts per SILVA mirror (default: 3)
  NO_COLOR             Disable ANSI colors
EOF
}

die() {
    printf '%bERROR:%b %s\n' "${COLOR_RED}" "${COLOR_RESET}" "$*" >&2
    if [[ -n "${ACTIVE_STEP_LABEL}" && -n "${CHECKPOINT_DIR:-}" ]]; then
        printf '%bRESUME:%b Stage "%s" is incomplete; prior validated stages remain in %s\n' \
            "${COLOR_YELLOW}" "${COLOR_RESET}" \
            "${ACTIVE_STEP_LABEL}" "${CHECKPOINT_DIR}" >&2
    fi
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

checkpoint_file() {
    printf '%s/%s.done' "${CHECKPOINT_DIR}" "$1"
}

checkpoint_exists() {
    [[ -f "$(checkpoint_file "$1")" ]]
}

write_checkpoint() {
    local key="$1"
    local label="$2"
    local destination
    local temporary
    destination="$(checkpoint_file "${key}")"
    temporary="${destination}.tmp"
    {
        printf 'checkpoint=%s\n' "${key}"
        printf 'label=%s\n' "${label}"
        printf 'profile=%s\n' "${PROFILE}"
        printf 'completed_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } >"${temporary}"
    mv "${temporary}" "${destination}"
}

progress_prefix() {
    local completed="$1"
    local percent=$((completed * 100 / TOTAL_STEPS))
    printf '[%02d/%02d %3d%%]' "${CURRENT_STEP}" "${TOTAL_STEPS}" "${percent}"
}

stage_begin() {
    ACTIVE_STEP_KEY="$1"
    ACTIVE_STEP_LABEL="$2"
    CURRENT_STEP=$((CURRENT_STEP + 1))
    ACTIVE_STEP_STARTED="${SECONDS}"
    printf '%b%s START%b %s\n' \
        "${COLOR_BLUE}" "$(progress_prefix $((CURRENT_STEP - 1)))" \
        "${COLOR_RESET}" "${ACTIVE_STEP_LABEL}"
    if checkpoint_exists "${ACTIVE_STEP_KEY}"; then
        warning "Checkpoint exists but validation failed; rebuilding ${ACTIVE_STEP_LABEL}"
    fi
}

stage_complete() {
    local elapsed=$((SECONDS - ACTIVE_STEP_STARTED))
    write_checkpoint "${ACTIVE_STEP_KEY}" "${ACTIVE_STEP_LABEL}"
    printf '%b%s DONE%b  %s %b(elapsed %s)%b\n' \
        "${COLOR_GREEN}" "$(progress_prefix "${CURRENT_STEP}")" \
        "${COLOR_RESET}" "${ACTIVE_STEP_LABEL}" "${COLOR_DIM}" \
        "$(format_duration "${elapsed}")" "${COLOR_RESET}"
}

stage_skip() {
    local key="$1"
    local label="$2"
    local detail="${3:-validated output exists}"
    CURRENT_STEP=$((CURRENT_STEP + 1))
    write_checkpoint "${key}" "${label}"
    printf '%b%s SKIP%b  %s %b(%s)%b\n' \
        "${COLOR_GREEN}" "$(progress_prefix "${CURRENT_STEP}")" \
        "${COLOR_RESET}" "${label}" "${COLOR_DIM}" "${detail}" \
        "${COLOR_RESET}"
}

resume_kraken_partial_downloads() {
    local root="${KRAKEN_PARTIAL_ROOT:-}"
    local partial final relative url partial_bytes
    [[ -n "${root}" && -d "${root}" ]] || return 0
    while IFS= read -r -d '' partial; do
        final="${partial%.tmp}"
        relative="${final#"${root}/"}"
        url="https://ftp.ncbi.nlm.nih.gov/${relative}"
        if gzip -t "${partial}" 2>/dev/null; then
            info "Promoting complete Kraken2 partial: ${relative}"
            mv "${partial}" "${final}"
            continue
        fi
        partial_bytes="$(wc -c <"${partial}" | tr -d '[:space:]')"
        info "Resuming interrupted Kraken2 file: ${relative} (${partial_bytes} bytes retained)"
        if ! curl \
            --fail \
            --location \
            --retry 20 \
            --retry-all-errors \
            --retry-delay 5 \
            --continue-at - \
            --output "${partial}" \
            "${url}"; then
            return 1
        fi
        if ! gzip -t "${partial}"; then
            warning "Resumed file failed gzip validation: ${partial}"
            return 1
        fi
        success "Recovered and validated ${relative}"
        mv "${partial}" "${final}"
    done < <(
        find "${root}" -type f -name '*.gz.tmp' ! -name '._*' -print0
    )
}

report_kraken_library_progress() {
    local root="${KRAKEN_PARTIAL_ROOT:-}"
    local manifest expected completed percent
    [[ -n "${root}" ]] || return 0
    manifest="${root}/manifest.txt"
    [[ -f "${manifest}" ]] || return 0
    expected="$(wc -l <"${manifest}" | tr -d '[:space:]')"
    [[ "${expected}" =~ ^[1-9][0-9]*$ ]] || return 0
    completed="$(
        find "${root}" -type f -name '*_genomic.fna.gz' -size +0c |
            wc -l |
            tr -d '[:space:]'
    )"
    percent=$((completed * 100 / expected))
    info "Library files: ${completed}/${expected} downloaded (${percent}%); completed files are retained"
}

retry_download() {
    local attempt=1
    local max_attempts="${TRIARCHY_DOWNLOAD_ATTEMPTS:-20}"
    local delay
    resume_kraken_partial_downloads ||
        info "Initial partial-file recovery did not complete; k2 will retry it"
    while true; do
        report_kraken_library_progress
        info "Download attempt ${attempt}/${max_attempts}"
        if "$@"; then
            return 0
        fi
        if ((attempt >= max_attempts)); then
            die "download failed after ${max_attempts} attempts: $*"
        fi
        resume_kraken_partial_downloads ||
            info "Partial-file recovery did not complete; the next k2 attempt will retry it"
        delay=$((attempt * 10))
        warning "Download interrupted after attempt ${attempt}/${max_attempts}"
        info "Retrying in ${delay}s; completed files remain checkpointed"
        sleep "${delay}"
        attempt=$((attempt + 1))
    done
}

while (($#)); do
    case "$1" in
        --database-dir)
            (($# >= 2)) || die "--database-dir requires a directory"
            DATABASE_DIR="$2"
            shift 2
            ;;
        --env-prefix)
            (($# >= 2)) || die "--env-prefix requires a directory"
            ENV_PREFIX="$2"
            shift 2
            ;;
        --cpus)
            (($# >= 2)) || die "--cpus requires an integer"
            CPUS="$2"
            shift 2
            ;;
        --silva-archive)
            (($# >= 2)) || die "--silva-archive requires a file"
            SILVA_ARCHIVE="$2"
            shift 2
            ;;
        --kraken-archive)
            (($# >= 2)) || die "--kraken-archive requires a file"
            KRAKEN_ARCHIVE="$2"
            shift 2
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

[[ "${CPUS}" =~ ^[1-9][0-9]*$ ]] || die "--cpus must be a positive integer"
[[ -x "${ENV_PREFIX}/bin/python3" ]] ||
    die "dependency environment is missing; run ./install.sh first"

if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    command -v brew >/dev/null 2>&1 ||
        die "Apple-silicon installation requires Homebrew; install it, then run: brew install bowtie2 bracken"
    HOMEBREW_PREFIX="$(brew --prefix)"
    export PATH="${HOMEBREW_PREFIX}/bin:${ENV_PREFIX}/bin:${PATH}"
    for native_command in bowtie2-build; do
        command -v "${native_command}" >/dev/null 2>&1 ||
            die "missing native ${native_command}; run: brew install bowtie2 bracken"
    done
else
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi
for command in datasets bowtie2-build \
    humann_databases metaphlan curl gzip unzip 7z sha256sum tar; do
    command -v "${command}" >/dev/null 2>&1 ||
        die "required dependency not found: ${command}"
done

bowtie2-build --version >/dev/null 2>&1 ||
    die "bowtie2-build is installed but cannot execute on this platform"

if [[ "${PROFILE}" == "standard" && "${ASSUME_YES}" == false ]]; then
    if [[ ! -t 0 ]]; then
        die "standard databases require substantial disk and network use; rerun with --yes"
    fi
    printf '%s\n' \
        "The standard reference installation can exceed 500 GB and take many hours." \
        "It downloads the prebuilt Kraken2 PlusPFP/Bracken database and other third-party data."
    read -r -p "Continue? [y/N] " answer
    [[ "${answer}" =~ ^[Yy]$ ]] || die "database installation cancelled"
fi

mkdir -p "${DATABASE_DIR}"/{downloads,host,kraken2,metaphlan,humann,silva,its}
DATABASE_DIR="$(cd -- "${DATABASE_DIR}" && pwd -P)"
kraken_db_dir="${DATABASE_DIR}/kraken2"
mkdir -p "${kraken_db_dir}"
TOTAL_STEPS=8
CHECKPOINT_DIR="${DATABASE_DIR}/.triarchy-install/checkpoints/${PROFILE}"
mkdir -p "${CHECKPOINT_DIR}"
readonly CHECKPOINT_DIR

installation_failed() {
    local status="$?"
    trap - ERR
    warning "Installation stopped during: ${ACTIVE_STEP_LABEL:-initialization}"
    info "Validated checkpoints remain in ${CHECKPOINT_DIR}"
    info "Rerun the same command to resume; completed stages will be skipped"
    exit "${status}"
}
trap installation_failed ERR

info "Database destination: ${DATABASE_DIR}"
info "Checkpoint directory: ${CHECKPOINT_DIR}"
info "Plan: ${TOTAL_STEPS} validated stages; existing outputs will be adopted"
checkpoint_count="$(
    find "${CHECKPOINT_DIR}" -maxdepth 1 -type f -name '*.done' |
        wc -l |
        tr -d '[:space:]'
)"
info "Checkpoint receipts found: ${checkpoint_count}; outputs will be revalidated"

if [[ "${PROFILE}" == "standard" ]]; then
    minimum_free_gb="${TRIARCHY_MIN_FREE_GB:-600}"
    available_kb="$(df -Pk "${DATABASE_DIR}" | awk 'NR == 2 { print $4 }')"
    required_kb=$((minimum_free_gb * 1024 * 1024))
    if ((available_kb < required_kb)) &&
        [[ "${TRIARCHY_ALLOW_LOW_DISK:-0}" != "1" ]]; then
        available_gb=$((available_kb / 1024 / 1024))
        if ((checkpoint_count > 0)) ||
            [[ -s "${DATABASE_DIR}/kraken2/taxonomy/nodes.dmp" ||
               -s "${DATABASE_DIR}/host/${TRIARCHY_HUMAN_ACCESSION:-GCF_000001405.40}_genomic.fna" ||
               -s "${DATABASE_DIR}/host/grch38.1.bt2" ||
               -s "${DATABASE_DIR}/host/grch38.1.bt2l" ]]; then
            warning "Resuming with ${available_gb} GB free (new installs require ${minimum_free_gb} GB)"
            warning "Monitor free space; the Kraken2 archive and extracted database coexist during installation"
        else
            die "standard profile requires at least ${minimum_free_gb} GB free; ${available_gb} GB available. Set TRIARCHY_ALLOW_LOW_DISK=1 to override."
        fi
    fi
fi
mkdir -p "${PROJECT_DIR}/references" "${PROJECT_DIR}/conf"
readonly MANIFEST="${PROJECT_DIR}/references/manifest.tsv"
MANIFEST_WORK="${CHECKPOINT_DIR}/manifest.tsv.partial"
printf 'resource\tversion\tsource\tpath\tsha256\tinstalled_at\n' >"${MANIFEST_WORK}"

record() {
    local resource="$1"
    local version="$2"
    local source="$3"
    local path="$4"
    local digest
    if [[ -f "${path}" ]]; then
        digest="$(sha256sum "${path}" | awk '{print $1}')"
    else
        digest="directory"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${resource}" "${version}" "${source}" "${path}" "${digest}" \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${MANIFEST_WORK}"
}

download() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"
    if [[ -s "${destination}" ]]; then
        info "Using existing $(basename "${destination}")"
        return
    fi
    info "Downloading ${url}"
    curl --fail --location \
        --retry "${TRIARCHY_DOWNLOAD_ATTEMPTS:-20}" \
        --retry-all-errors \
        --retry-delay 10 \
        --connect-timeout 60 \
        --continue-at - --output "${partial}" "${url}"
    mv "${partial}" "${destination}"
}

download_pinned_md5() {
    local url="$1"
    local destination="$2"
    local expected_md5="$3"
    local actual_md5=""

    if [[ -s "${destination}" ]]; then
        actual_md5="$(md5sum "${destination}" | awk '{print $1}')"
        if [[ "${actual_md5}" == "${expected_md5}" ]]; then
            info "Using verified $(basename "${destination}")"
            return 0
        fi
        warning "Discarding corrupt $(basename "${destination}") (checksum ${actual_md5})"
        rm -f -- "${destination}" "${destination}.part"
    fi

    download "${url}" "${destination}"
    actual_md5="$(md5sum "${destination}" | awk '{print $1}')"
    if [[ "${actual_md5}" != "${expected_md5}" ]]; then
        rm -f -- "${destination}" "${destination}.part"
        die "checksum mismatch for $(basename "${destination}"): ${actual_md5}"
    fi
    success "Verified $(basename "${destination}")"
}

VERIFIED_DOWNLOAD_SOURCE=""
download_verified_md5() {
    local destination="$1"
    local expected_md5="$2"
    shift 2
    local mirror_count="$#"
    local mirror_index=0
    local url
    local partial
    local source_marker
    local partial_source
    local actual_md5

    ((mirror_count > 0)) || die "no download mirrors configured for ${destination}"

    # Older TRIARCHY releases used one unlabelled partial. It belongs to a
    # different SILVA artifact and must not be combined with the new mirrors.
    if [[ -e "${destination}.part" ]]; then
        warning "Discarding an unverified partial from the legacy downloader"
        rm -f -- "${destination}.part"
    fi

    for url in "$@"; do
        mirror_index=$((mirror_index + 1))
        partial="${destination}.part.${mirror_index}"
        source_marker="${partial}.url"
        partial_source=""
        if [[ -s "${source_marker}" ]]; then
            partial_source="$(head -n 1 "${source_marker}")"
        fi
        if [[ -e "${partial}" && "${partial_source}" != "${url}" ]]; then
            warning "Resetting partial data because mirror ${mirror_index} changed"
            rm -f -- "${partial}" "${source_marker}"
        fi
        printf '%s\n' "${url}" >"${source_marker}"

        info "Downloading mirror ${mirror_index}/${mirror_count}: ${url}"
        if curl \
            --fail \
            --location \
            --http1.1 \
            --retry "${TRIARCHY_MIRROR_ATTEMPTS:-3}" \
            --retry-all-errors \
            --retry-delay 5 \
            --connect-timeout 20 \
            --speed-limit 1024 \
            --speed-time 120 \
            --continue-at - \
            --output "${partial}" \
            "${url}"; then
            actual_md5="$(md5sum "${partial}" | awk '{print $1}')"
            if [[ "${actual_md5}" == "${expected_md5}" ]]; then
                mv "${partial}" "${destination}"
                rm -f -- "${source_marker}"
                VERIFIED_DOWNLOAD_SOURCE="${url}"
                success "Verified $(basename "${destination}") from mirror ${mirror_index}"
                return 0
            fi
            warning "Mirror ${mirror_index} returned checksum ${actual_md5}; expected ${expected_md5}"
            rm -f -- "${partial}" "${source_marker}"
        else
            warning "Mirror ${mirror_index}/${mirror_count} is unavailable; trying the next source"
        fi
    done

    die "all ${mirror_count} verified mirrors failed for $(basename "${destination}")"
}

bowtie_index_complete() {
    local index_base="$1"
    local extension
    local suffix
    if [[ -s "${index_base}.1.bt2l" ]]; then
        extension="bt2l"
    elif [[ -s "${index_base}.1.bt2" ]]; then
        extension="bt2"
    else
        return 1
    fi
    for suffix in 1 2 3 4 rev.1 rev.2; do
        [[ -s "${index_base}.${suffix}.${extension}" ]] ||
            return 1
    done
}

build_bowtie2_index() {
    local reference="$1"
    local index_base="$2"
    if (bowtie2-build --help 2>&1 || true) | grep -q -- '--threads'; then
        bowtie2-build --threads "${CPUS}" "${reference}" "${index_base}"
    else
        warning "bowtie2-build does not support --threads; building with default single-threaded mode"
        bowtie2-build "${reference}" "${index_base}"
    fi
}

    human_accession="${TRIARCHY_HUMAN_ACCESSION:-GCF_000001405.40}"
    human_source="https://www.ncbi.nlm.nih.gov/datasets/genome/${human_accession}/"
    human_zip="${DATABASE_DIR}/downloads/${human_accession}.zip"
    human_fasta="${DATABASE_DIR}/host/${human_accession}_genomic.fna"
    host_index_base="${DATABASE_DIR}/host/grch38"
    host_reference_label="Human host reference ${human_accession}"
    host_index_label="Bowtie2 GRCh38 host index"

if [[ -s "${human_fasta}" ]]; then
    stage_skip host_reference "${host_reference_label}"
else
    stage_begin host_reference "${host_reference_label}"
        info "Downloading human reference ${human_accession} from NCBI"
        datasets download genome accession "${human_accession}" \
            --include genome --filename "${human_zip}"
        unzip -o "${human_zip}" -d "${DATABASE_DIR}/host/ncbi_dataset"
        source_fasta="$(
            find "${DATABASE_DIR}/host/ncbi_dataset" -type f \
                -name '*_genomic.fna' -print -quit
        )"
        [[ -n "${source_fasta}" ]] || die "NCBI package did not contain a genomic FASTA"
        cp "${source_fasta}" "${human_fasta}"
    stage_complete
fi

if bowtie_index_complete "${host_index_base}"; then
    stage_skip host_bowtie2 "${host_index_label}"
else
    stage_begin host_bowtie2 "${host_index_label}"
    info "Building Bowtie2 host index"
    build_bowtie2_index "${human_fasta}" "${host_index_base}"
    bowtie_index_complete "${host_index_base}" ||
        die "Bowtie2 build exited without a complete six-file index"
    stage_complete
fi
record host "${human_accession}" "${human_source}" "${human_fasta}"

silva_primary_url="https://ftp.arb-silva.de/release_138.2/DADA2/1.36.0/SSU/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_mirror_url="https://zenodo.org/records/20955974/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1"
silva_urls=()
if [[ -n "${TRIARCHY_SILVA_URL:-}" ]]; then
    silva_urls+=("${TRIARCHY_SILVA_URL}")
fi
silva_urls+=("${silva_primary_url}" "${silva_mirror_url}")
silva_fasta="${DATABASE_DIR}/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_expected_md5="41f1d231ce5079b88f4c01742c288ba9"
silva_legacy_md5="1764e2a36b4500ccb1c7d5261948a414"
silva_source="validated existing SILVA 138.2 artifact"
silva_md5=""
if [[ -n "${SILVA_ARCHIVE}" ]]; then
    [[ -f "${SILVA_ARCHIVE}" ]] ||
        die "SILVA archive does not exist: ${SILVA_ARCHIVE}"
    SILVA_ARCHIVE="$(cd -- "$(dirname -- "${SILVA_ARCHIVE}")" && pwd -P)/$(basename -- "${SILVA_ARCHIVE}")"
    silva_source="${SILVA_ARCHIVE}"
fi
if [[ -s "${silva_fasta}" ]]; then
    silva_md5="$(md5sum "${silva_fasta}" | awk '{print $1}')"
fi
if [[ "${silva_md5}" == "${silva_expected_md5}" || "${silva_md5}" == "${silva_legacy_md5}" ]]; then
    stage_skip silva_taxonomy "SILVA 138.2 DADA2 taxonomy"
else
    [[ -e "${silva_fasta}" ]] &&
        die "Existing SILVA file failed validation: ${silva_fasta}"
    stage_begin silva_taxonomy "SILVA 138.2 DADA2 taxonomy"
    if [[ -n "${SILVA_ARCHIVE}" ]]; then
        silva_archive_md5="$(md5sum "${SILVA_ARCHIVE}" | awk '{print $1}')"
        [[ "${silva_archive_md5}" == "${silva_expected_md5}" || "${silva_archive_md5}" == "${silva_legacy_md5}" ]] ||
            die "SILVA archive checksum mismatch: ${silva_archive_md5}"
        info "Installing verified SILVA archive from ${SILVA_ARCHIVE}"
        cp "${SILVA_ARCHIVE}" "${silva_fasta}.part"
        mv "${silva_fasta}.part" "${silva_fasta}"
    else
        download_verified_md5 "${silva_fasta}" "${silva_expected_md5}" "${silva_urls[@]}"
        silva_source="${VERIFIED_DOWNLOAD_SOURCE}"
    fi
    silva_md5="$(md5sum "${silva_fasta}" | awk '{print $1}')"
    [[ "${silva_md5}" == "${silva_expected_md5}" || "${silva_md5}" == "${silva_legacy_md5}" ]] ||
        die "SILVA checksum mismatch: ${silva_md5}"
    stage_complete
fi
record silva "138.2" "${silva_source}" "${silva_fasta}"

its_url="${TRIARCHY_ITS_URL:-https://sisu.ut.ee/wp-content/uploads/sites/643/DADA2_EUK_ITS_99_v2.0.zip}"
its_zip="${DATABASE_DIR}/downloads/DADA2_EUK_ITS_99_v2.0.zip"
its_fasta="$(
    find "${DATABASE_DIR}/its" -type f \
        \( -name '*.fa' -o -name '*.fasta' -o -name '*.fa.gz' -o -name '*.fasta.gz' \) \
        -print -quit
)"
if [[ -n "${its_fasta}" && -s "${its_fasta}" ]]; then
    stage_skip its_taxonomy "EUKARYOME ITS taxonomy"
else
    stage_begin its_taxonomy "EUKARYOME ITS taxonomy"
    download "${its_url}" "${its_zip}"
    its_fasta="$(
        find "${DATABASE_DIR}/its" -type f \
            \( -name '*.fa' -o -name '*.fasta' -o -name '*.fa.gz' -o -name '*.fasta.gz' \) \
            -print -quit
    )"
    if [[ -z "${its_fasta}" || ! -s "${its_fasta}" ]]; then
        if ! unzip -o "${its_zip}" -d "${DATABASE_DIR}/its"; then
            info "System unzip cannot read this archive; retrying with 7-Zip"
            7z x -y -o"${DATABASE_DIR}/its" "${its_zip}"
        fi
    fi
    its_fasta="$(
        find "${DATABASE_DIR}/its" -type f \
            \( -name '*.fa' -o -name '*.fasta' -o -name '*.fa.gz' -o -name '*.fasta.gz' \) \
            -print -quit
    )"
    [[ -n "${its_fasta}" && -s "${its_fasta}" ]] ||
        die "ITS archive did not contain a nonempty FASTA file"
    stage_complete
fi
record its_taxonomy "EUKARYOME-2.0-99pct" "${its_url}" "${its_fasta}"

kraken_taxonomy_complete() {
    local taxonomy="${kraken_db_dir}/taxonomy"
    [[ -s "${taxonomy}/nodes.dmp" &&
       -s "${taxonomy}/names.dmp" ]]
}

kraken_library_complete() {
    local library="$1"
    local directory="${kraken_db_dir}/library/${library}"
    [[ -s "${directory}/library.fna" &&
       -s "${directory}/prelim_map.txt" ]]
}

kraken_database_complete() {
    [[ -s "${kraken_db_dir}/hash.k2d" &&
       -s "${kraken_db_dir}/opts.k2d" &&
       -s "${kraken_db_dir}/taxo.k2d" ]]
}

bracken_database_complete() {
    local distribution
    distribution="$(
        find "${kraken_db_dir}" -maxdepth 1 \
            -name 'database*mers.kmer_distrib' -type f -size +0c \
            -print -quit
    )"
    [[ -n "${distribution}" ]]
}

kraken_prebuilt_complete() {
    [[ -f "${kraken_db_dir}/.triarchy-pluspfp-20260626.complete" ]] &&
        kraken_database_complete &&
        bracken_database_complete
}


    kraken_version="PlusPFP-20260626"
    kraken_url="${TRIARCHY_KRAKEN_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20260626.tar.gz}"
    kraken_md5_url="${TRIARCHY_KRAKEN_MD5_URL:-https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20260626/pluspfp.md5}"
    kraken_expected_size="${TRIARCHY_KRAKEN_SIZE:-184467543914}"
    if [[ -n "${KRAKEN_ARCHIVE}" ]]; then
        [[ -s "${KRAKEN_ARCHIVE}" ]] ||
            die "Kraken2 archive does not exist or is empty: ${KRAKEN_ARCHIVE}"
        kraken_archive="$(cd -- "$(dirname -- "${KRAKEN_ARCHIVE}")" && pwd -P)/$(basename -- "${KRAKEN_ARCHIVE}")"
    else
        kraken_archive="${DATABASE_DIR}/downloads/k2_pluspfp_20260626.tar.gz"
    fi
    kraken_md5="${DATABASE_DIR}/downloads/pluspfp_20260626.md5"

    if kraken_prebuilt_complete; then
        stage_skip kraken_prebuilt \
            "Kraken2 ${kraken_version} database and Bracken distributions"
    else
        stage_begin kraken_prebuilt \
            "Kraken2 ${kraken_version} database and Bracken distributions"
        download "${kraken_md5_url}" "${kraken_md5}"
        if [[ -n "${KRAKEN_ARCHIVE}" ]]; then
            info "Using existing $(basename "${kraken_archive}")"
        else
            download "${kraken_url}" "${kraken_archive}"
        fi

        actual_archive_size="$(
            wc -c <"${kraken_archive}" | tr -d '[:space:]'
        )"
        if [[ "${actual_archive_size}" != "${kraken_expected_size}" ]]; then
            die "Kraken2 archive is incomplete: ${actual_archive_size}/${kraken_expected_size} bytes. Allow the separate download to finish, then rerun this command."
        fi

        expected_archive_md5="$(
            awk '$2 == "k2_pluspfp_20260626.tar.gz" { print $1; exit }' \
                "${kraken_md5}"
        )"
        [[ "${expected_archive_md5}" =~ ^[0-9a-fA-F]{32}$ ]] ||
            die "Kraken2 checksum manifest does not contain the archive checksum"
        actual_archive_md5="$(md5sum "${kraken_archive}" | awk '{print $1}')"
        [[ "${actual_archive_md5}" == "${expected_archive_md5}" ]] ||
            die "Kraken2 archive checksum mismatch: ${actual_archive_md5}"

        info "Extracting the prebuilt Kraken2/Bracken database"
        tar -xzf "${kraken_archive}" -C "${kraken_db_dir}"
        (
            cd -- "${kraken_db_dir}"
            awk '$2 != "k2_pluspfp_20260626.tar.gz"' "${kraken_md5}" |
                md5sum --check --strict -
        )
        kraken_database_complete ||
            die "Prebuilt archive did not contain a complete Kraken2 database"
        bracken_database_complete ||
            die "Prebuilt archive did not contain Bracken distributions"
        touch "${kraken_db_dir}/.triarchy-pluspfp-20260626.complete"
        stage_complete
    fi

metaphlan_index="mpa_vJun23_CHOCOPhlAnSGB_202307"
metaphlan_db_dir="${DATABASE_DIR}/metaphlan"
metaphlan_base_url="https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"
metaphlan_main_url="${metaphlan_base_url}/${metaphlan_index}.tar"
metaphlan_bt2_url="${metaphlan_base_url}/bowtie2_indexes/${metaphlan_index}_bt2.tar"
metaphlan_main_archive="${DATABASE_DIR}/downloads/${metaphlan_index}.tar"
metaphlan_bt2_archive="${DATABASE_DIR}/downloads/${metaphlan_index}_bt2.tar"
metaphlan_main_md5="f883708ea3c26b59d5f3410b72256bb2"
metaphlan_bt2_md5="d306a65273bb9b00ab6977c1e264cdfc"

metaphlan_database_complete() {
    local index_base="${metaphlan_db_dir}/${metaphlan_index}"
    local extension=""
    local forward_size
    local reverse_size
    local size_difference

    [[ -s "${metaphlan_db_dir}/${metaphlan_index}.pkl" ]] || return 1
    bowtie_index_complete "${index_base}" || return 1
    if [[ -s "${index_base}.1.bt2l" ]]; then
        extension="bt2l"
    else
        extension="bt2"
    fi
    forward_size="$(wc -c <"${index_base}.1.${extension}" | tr -d '[:space:]')"
    reverse_size="$(wc -c <"${index_base}.rev.1.${extension}" | tr -d '[:space:]')"
    if ((forward_size >= reverse_size)); then
        size_difference=$((forward_size - reverse_size))
    else
        size_difference=$((reverse_size - forward_size))
    fi
    ((size_difference <= 1000))
}

if metaphlan_database_complete; then
    stage_skip metaphlan_markers \
        "MetaPhlAn ${metaphlan_index} marker database"
else
    stage_begin metaphlan_markers \
        "MetaPhlAn ${metaphlan_index} marker database"
    download_pinned_md5 \
        "${metaphlan_main_url}" "${metaphlan_main_archive}" "${metaphlan_main_md5}"
    download_pinned_md5 \
        "${metaphlan_bt2_url}" "${metaphlan_bt2_archive}" "${metaphlan_bt2_md5}"
    info "Extracting verified MetaPhlAn database archives"
    tar -xf "${metaphlan_bt2_archive}" -C "${metaphlan_db_dir}"
    tar -xf "${metaphlan_main_archive}" -C "${metaphlan_db_dir}"
    metaphlan_database_complete ||
        die "MetaPhlAn archives did not produce a complete usable database"
    printf '%s\n' "${metaphlan_index}" >"${metaphlan_db_dir}/mpa_latest"
    stage_complete
fi

    humann_chocophlan_file="$(
        find "${DATABASE_DIR}/humann" -type f -name '*.ffn.gz' \
            -size +0c -print -quit
    )"
    if [[ -n "${humann_chocophlan_file}" ]]; then
        stage_skip humann_chocophlan "HUMAnN ChocoPhlAn database"
    else
        stage_begin humann_chocophlan "HUMAnN ChocoPhlAn database"
        info "Downloading full HUMAnN ChocoPhlAn"
        humann_databases --download chocophlan full \
            "${DATABASE_DIR}/humann" --update-config no
        humann_chocophlan_file="$(
            find "${DATABASE_DIR}/humann" -type f -name '*.ffn.gz' \
                -size +0c -print -quit
        )"
        [[ -n "${humann_chocophlan_file}" ]] ||
            die "HUMAnN ChocoPhlAn download did not pass validation"
        stage_complete
    fi

    humann_uniref_file="$(
        find "${DATABASE_DIR}/humann" -type f -name '*.dmnd' \
            -size +0c -print -quit
    )"
    if [[ -n "${humann_uniref_file}" ]]; then
        stage_skip humann_uniref90 "HUMAnN UniRef90 Diamond database"
    else
        stage_begin humann_uniref90 "HUMAnN UniRef90 Diamond database"
        info "Downloading HUMAnN UniRef90 Diamond"
        humann_databases --download uniref uniref90_diamond \
            "${DATABASE_DIR}/humann" --update-config no
        humann_uniref_file="$(
            find "${DATABASE_DIR}/humann" -type f -name '*.dmnd' \
                -size +0c -print -quit
        )"
        [[ -n "${humann_uniref_file}" ]] ||
            die "HUMAnN UniRef90 download did not pass validation"
        stage_complete
    fi

humann_nucleotide="$(dirname -- "${humann_chocophlan_file}")"
humann_protein="$(dirname -- "${humann_uniref_file}")"
    record kraken2 "${kraken_version}" "${kraken_url}" \
        "${kraken_db_dir}"
record humann_chocophlan "${PROFILE}" "humann_databases" "${humann_nucleotide}"
record humann_uniref "${PROFILE}" "humann_databases" "${humann_protein}"
record metaphlan "${metaphlan_index}" "${metaphlan_base_url}" "${metaphlan_db_dir}"

config="${PROJECT_DIR}/conf/references.config"
config_work="${CHECKPOINT_DIR}/references.config.partial"
cat >"${config_work}" <<EOF
params {
    host_index = '${host_index_base}'
    kraken2_db = '${kraken_db_dir}'
    kraken2_memory_mapping = true
    humann_db = '${humann_nucleotide}'
    humann_protein_db = '${humann_protein}'
    metaphlan_db = '${metaphlan_db_dir}'
    metaphlan_index = '${metaphlan_index}'
    s16_taxonomy = '${silva_fasta}'
    its_taxonomy = '${its_fasta}'
}
EOF

mv "${config_work}" "${config}"
mv "${MANIFEST_WORK}" "${MANIFEST}"
success "Reference installation complete in $(format_duration $((SECONDS - INSTALL_STARTED)))"
info "Progress: ${TOTAL_STEPS}/${TOTAL_STEPS} stages validated (100%)"
printf 'Manifest: %s\nConfig:   %s\n' "${MANIFEST}" "${config}"
