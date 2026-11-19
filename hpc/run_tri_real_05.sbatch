#!/bin/bash
#SBATCH --job-name=triarchy-tri-real-05
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=triarchy-tri-real-05-%j.log
#
# Nextflow head job for the TRI-REAL-05 cohort. It needs little CPU or memory
# but must stay alive for the whole run; Nextflow submits every workflow task
# as its own SLURM job using the `slurm` profile plus hpc/site.config.
#
# Before submitting (from the TRIARCHY repository root):
#   1. triarchy databases check                                   -> "ok": true
#   2. triarchy real-world --dataset TRI-REAL-05 --download-only  (login node)
#   3. cp hpc/site.config.example hpc/site.config                 and edit it
#
# Submit:
#   sbatch hpc/run_tri_real_05.sbatch
#
# Arguments after the script name are passed to `triarchy real-world`. Pass the
# same cohort options used for the download, for example:
#   sbatch hpc/run_tri_real_05.sbatch --max-biological-per-plate 20
#
# Resubmit the same command to resume; finished tasks are not rerun.

set -Eeo pipefail

# Shared scratch visible from every compute node; override with WORK_DIR=...
WORK_DIR="${WORK_DIR:-/scratch/${USER}/triarchy-tri-real-05-work}"
OUTDIR="${OUTDIR:-results/TRI-REAL-05}"
SITE_CONFIG="${SITE_CONFIG:-hpc/site.config}"

# Batch jobs do not read ~/.bashrc, so load mamba's shell hook explicitly.
if mamba shell hook --shell bash >/dev/null 2>&1; then
    eval "$(mamba shell hook --shell bash)"          # mamba 2.x
else
    source "$(conda info --base)/etc/profile.d/conda.sh"   # mamba 1.x
    source "$(conda info --base)/etc/profile.d/mamba.sh"
fi
mamba activate TRIARCHY
# Conda's activation scripts are not nounset-safe; enable it afterwards.
set -u

cd "${SLURM_SUBMIT_DIR}"
if [[ ! -f "${SITE_CONFIG}" ]]; then
    echo "ERROR: ${SITE_CONFIG} not found; copy hpc/site.config.example and set your account and partition." >&2
    exit 2
fi

triarchy real-world --dataset TRI-REAL-05 --skip-download \
    --profile slurm \
    --site-config "${SITE_CONFIG}" \
    --work-dir "${WORK_DIR}" \
    --outdir "${OUTDIR}" \
    "$@"
