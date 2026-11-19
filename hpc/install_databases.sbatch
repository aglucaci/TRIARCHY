#!/bin/bash
#SBATCH --job-name=triarchy-databases
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=triarchy-databases-%j.log
#
# Install the TRIARCHY production reference databases (several hundred GB).
# This job needs internet access. If your compute nodes have none, run the
# same `triarchy databases install` command on a login or transfer node
# inside tmux or screen instead.
#
# Submit from the TRIARCHY repository root:
#   sbatch hpc/install_databases.sbatch /path/to/project/storage/triarchy-refs
#
# Resubmit the same command to resume an interrupted installation.

set -Eeo pipefail

REFS_DIR="${1:?usage: sbatch hpc/install_databases.sbatch /path/to/reference/directory}"

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
mkdir -p "${REFS_DIR}"
triarchy databases install --destination "${REFS_DIR}" --cpus "${SLURM_CPUS_PER_TASK:-8}" --yes
triarchy databases check
