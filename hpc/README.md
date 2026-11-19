# Running TRIARCHY on a SLURM cluster

Templates for running the `TRI-REAL-05` cohort (and, by extension, any
TRIARCHY study) on SLURM. Run every command from the TRIARCHY repository root.

| File | Purpose |
| --- | --- |
| [`site.config.example`](site.config.example) | Your SLURM account, partition and concurrency limits. Copy to `hpc/site.config` (ignored by git). |
| [`install_databases.sbatch`](install_databases.sbatch) | Installs the reference databases as a batch job. |
| [`run_tri_real_05.sbatch`](run_tri_real_05.sbatch) | Nextflow head job for `TRI-REAL-05`. |

## Requirements

- Mamba (1.x or 2.x), Java 17+ and `git` on the login node.
- Shared POSIX storage visible at the same path from login and compute nodes,
  for the repository, references, reads, Nextflow work directory and results.
  Plan for about 1 TB: over 500 GB of references, 6 GB of reads for the default
  cohort, and work files.
- Compute nodes able to run 12-CPU, 96 GB jobs for 48 h (the `high` label).

## 1. Install the software (once)

```bash
git clone https://github.com/aglucaci/TRIARCHY.git
cd TRIARCHY
scripts/install_environment.sh --yes
mamba activate TRIARCHY
triarchy --version
```

This creates the `TRIARCHY` environment, installs the pinned Nextflow into
`.tools/bin`, and puts `triarchy` on the environment's PATH. If `triarchy` is
still not found, run `ln -sfn "$PWD/triarchy" "$CONDA_PREFIX/bin/triarchy"`.

## 2. Install the reference databases (once, several hours)

The installer needs internet access. If compute nodes have it, submit:

```bash
sbatch hpc/install_databases.sbatch /path/to/project/storage/triarchy-refs
```

Otherwise run it on a login or transfer node inside `tmux` or `screen`:

```bash
triarchy databases install --destination /path/to/project/storage/triarchy-refs --cpus 8 --yes
```

Rerun the same command to resume an interruption. When it finishes, this must
report `"ok": true`:

```bash
triarchy databases check
```

## 3. Download the cohort (login or transfer node)

```bash
triarchy real-world --dataset TRI-REAL-05 --download-only
```

Reads, metadata and the mock truth table land in `data/real_world/TRI-REAL-05/`
(about 6.1 GB for the default 190 samples). If any file fails, the command
lists it in `data/real_world/TRI-REAL-05/metadata/download_failures.tsv` and
exits with an error; rerun it until that file is gone. To change the cohort
size, see [Scaling the cohort](../docs/datasets.md#scaling-the-cohort), and
pass the same options again in step 5.

## 4. Configure your site

```bash
cp hpc/site.config.example hpc/site.config
```

Set `queue` and `clusterOptions` (account) in `hpc/site.config`. Add
`beforeScript = 'module load ...'` if your compute nodes need modules.

## 5. Submit the workflow

```bash
sbatch hpc/run_tri_real_05.sbatch
```

The job is the Nextflow head process: 2 CPUs, 8 GB, up to 7 days. It submits
every workflow task as its own SLURM job. Defaults can be overridden with
environment variables, and any arguments are passed to `triarchy real-world`:

```bash
WORK_DIR=/scratch/$USER/tri05-work OUTDIR=/project/results/TRI-REAL-05 \
  sbatch hpc/run_tri_real_05.sbatch --max-biological-per-plate 20
```

| Variable | Default |
| --- | --- |
| `WORK_DIR` | `/scratch/$USER/triarchy-tri-real-05-work` (must be shared storage) |
| `OUTDIR` | `results/TRI-REAL-05` |
| `SITE_CONFIG` | `hpc/site.config` |

Monitor with `squeue -u $USER` and `tail -f triarchy-tri-real-05-<jobid>.log`.
If the job stops or hits its time limit, submit the same command again:
`--resume` is on by default, so finished tasks are reused.

## 6. Results

| Path | Contents |
| --- | --- |
| `results/TRI-REAL-05/integrated/taxon_evidence.tsv` | Per-sample, per-taxon evidence: contamination score, detection power, cross-assay status, authenticity |
| `results/TRI-REAL-05/mock_evaluation/mock_evaluation.json` | ZymoBIOMICS truth evaluation: recall by assay, abundance and dilution; false support of unexpected taxa |
| `results/TRI-REAL-05/mock_evaluation/mock_taxon_calls.tsv` | One row per mock sample and expected genus |
| `results/TRI-REAL-05/pipeline_info/` | Nextflow trace, timeline, report and provenance |

See [docs/datasets.md](../docs/datasets.md) for the cohort design and
[docs/outputs.md](../docs/outputs.md) for every output file.

## Other studies

For your own data, copy `run_tri_real_05.sbatch` and replace the
`triarchy real-world ...` command with:

```bash
triarchy run --profile slurm \
    --site-config "${SITE_CONFIG}" \
    --work-dir "${WORK_DIR}" \
    --samples /path/to/samples.tsv \
    --controls /path/to/control_relationships.tsv \
    --outdir "${OUTDIR}" \
    --resume
```

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| `triarchy: command not found` | `mamba activate TRIARCHY`; if still missing, `ln -sfn "$PWD/triarchy" "$CONDA_PREFIX/bin/triarchy"` |
| `mamba activate` fails in the job | The script loads mamba's shell hook itself; check `mamba --version` works on compute nodes |
| `reference manifest not found` | Step 2 has not completed in this checkout |
| `hpc/site.config not found` | Step 4 |
| Tasks pending forever | Wrong `queue` or account in `hpc/site.config`, or the partition cannot satisfy 96 GB / 48 h |
| Nextflow cannot find work files on compute nodes | `WORK_DIR` is not on shared storage |
