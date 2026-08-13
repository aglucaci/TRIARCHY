# Running TRIARCHY

## Installation

For a lightweight local check that avoids the full Conda environment and all
reference databases, install only Nextflow and run the deterministic lite mode:

```bash
./install.sh --skip-dependencies --skip-databases --skip-tests --yes
./run_triarchy.sh --lite
```

Lite mode writes to `results/lite` and is intended for laptop/WSL smoke tests,
report inspection, and development. It is not suitable for biological
interpretation; see [lightweight.md](lightweight.md).

On an Apple-silicon Mac, `install.sh` runs `brew install bowtie2 bracken`
(Bracken brings Kraken2), then installs Conda activation hooks that prepend
the native formula binaries while the TRIARCHY environment is active. This
avoids incompatible Bioconda executables without overwriting Conda-managed
files. From the repository root, install the complete production environment
and database set. Mamba is the default environment manager. The standard
reference set can exceed 400 GB and take many hours:

```bash
conda activate TRIARCHY

./install.sh \
  --env-prefix "$CONDA_PREFIX" \
  --environment-manager mamba \
  --database-dir /Volumes/ExternalHD/TRIARCHY/references/data \
  --database-profile standard \
  --cpus 8 \
  --yes
```

To reuse an archive downloaded separately, add:

```bash
--kraken-archive /path/to/k2_pluspfp_20260626.tar.gz
```

Wait for a separately managed archive download to finish before running this
command. The installer rejects a partial archive by exact byte size before
performing the full MD5 validation.

Do not use `--skip-dependencies`. The installer automatically gives native
Homebrew executables priority over incompatible Bioconda binaries on
Apple-silicon systems.

### Database checkpoints and progress

Rerun the exact same installation command after an interruption. The database
installer validates existing outputs and skips completed work; it does not
restart from the human reference or redownload completed Kraken libraries.
Durable checkpoint receipts are stored outside the repository at:

```text
/Volumes/ExternalHD/TRIARCHY/references/data/.triarchy-install/checkpoints/standard/
```

The terminal reports every standard-profile phase as `[step/7 percent]`,
followed by `START`, `DONE`, or `SKIP` and elapsed time for completed phases.
There are separate checkpoints for the host reference, Bowtie2 index, SILVA,
ITS, the prebuilt Kraken2/Bracken database, and both HUMAnN databases. A
checkpoint is trusted only when its corresponding output also passes
validation. Existing valid output without a receipt is adopted and receives
one automatically.

ANSI colors are enabled when output is attached to a terminal. Set
`NO_COLOR=1` to disable them. Set `TRIARCHY_DOWNLOAD_ATTEMPTS` to change the
default eight retry attempts for a database download phase.

The standard profile downloads the Langmead Lab June 2026 PlusPFP archive over
resumable HTTPS. The installer verifies the official MD5 checksum before
extraction, then verifies every extracted Kraken2 and Bracken file. PlusPFP is
broader than Kraken2's preset named `standard`: it includes archaeal,
bacterial, viral, fungal, protozoan, plant, human, plasmid, and UniVec_Core
libraries so fungal WMS-ITS evidence is not silently omitted.

For an environment that already exists, configure only the native-tool bridge:

```bash
scripts/configure_conda_native_tools.sh \
  --env-prefix "$CONDA_PREFIX" \
  --install

conda deactivate
conda activate TRIARCHY
```

## Shell environment

From the repository root:

```bash
export PATH="$PWD/.tools/env/bin:$PWD/.tools/bin:$PATH"
export NXF_HOME="$PWD/.nextflow"
export NXF_VER="$(<.nextflow-version)"
```

## Complete stub execution

TRI-00 contains valid metadata but deliberately does not contain FASTQs. Use it
to verify every Nextflow process and integration output:

```bash
./run_triarchy.sh --stub
```

Stub profiles are deterministic orchestration fixtures. They are not simulated
sequencing results and must not be used for biological interpretation.

## Production prerequisites

The installer creates `.tools/env` with these programs:

- fastp;
- Bowtie2;
- Kraken2 and Bracken;
- HUMAnN;
- cutadapt;
- ITSxpress;
- R with the DADA2 package;
- Python 3.

`environment.lock.txt` records the exact resolved packages after installation.
The `docker`, `singularity`, and `slurm` profiles remain available for
deployment-specific container and cluster configurations.

## Required production inputs

- `--samples`: sample manifest defined in
  [data-contracts.md](data-contracts.md);
- `--controls`: sample-control edge table;
- `--host_index`: Bowtie2 host-index prefix;
- `--kraken2_db`: Kraken2/Bracken database directory;
- `--humann_db`: HUMAnN ChocoPhlAn directory;
- `--humann_protein_db`: HUMAnN UniRef Diamond directory;
- `--s16_taxonomy`: DADA2-compatible bacterial/archaeal taxonomy FASTA;
- `--its_taxonomy`: DADA2-compatible fungal taxonomy FASTA;
- 16S and ITS forward primers, supplied through the cutadapt parameters below.

Recommended: `--taxonomy_crosswalk`, an explicit canonical taxonomy crosswalk.
Optional: `--assay_performance`, a TSV overriding the versioned assay/kingdom
reliability, minimum depth, reference coverage, and false-positive parameters.

Optional: `--evidence_policy`, a JSON file overriding inference thresholds,
the prior, model version, and calibration mapping. The shipped mapping is
synthetic-only; fit a study-appropriate calibration with:

```bash
python3 workflow/scripts/triarchy_integrate.py fit-calibration \
  --evidence results/integrated/taxon_evidence.tsv \
  --truth benchmark_truth.tsv \
  --output fitted_calibration.json
```

Relative FASTQ paths are resolved relative to the sample manifest. Reference
paths are resolved before processes enter their isolated work directories.

## Preflight validation

Run preflight checks before scheduling a production job:

```bash
bin/triarchy validate \
  --input /data/study/samples.tsv \
  --paired-controls /data/study/control_relationships.tsv \
  --databases conf/databases.example.yml \
  --database-destination /refs/triarchy \
  --report results/preflight.json
```

The command validates the metadata contract, verifies FASTQ paths, checks
database paths, and verifies file checksums when `sha256` values are present in
the manifest. Use `--skip-fastq-check` only for deterministic stub fixtures
whose FASTQ paths are intentionally synthetic.

## Production command

```bash
./run_triarchy.sh \
  -profile local \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --host_index /refs/host/grch38 \
  --kraken2_db /refs/kraken2 \
  --humann_db /refs/humann/chocophlan \
  --humann_protein_db /refs/humann/uniref90 \
  --s16_taxonomy /refs/silva/silva_train_set.fa.gz \
  --its_taxonomy /refs/unite/unite_train_set.fa.gz \
  --cutadapt_16s_forward CCTACGGGNGGCWGCAG \
  --cutadapt_16s_reverse GACTACHVGGGTATCTAATCC \
  --cutadapt_its_forward GTGAATCATCGAATCTTTGAA \
  --cutadapt_its_reverse TCCTCCGCTTATTGATATGC \
  --taxonomy_crosswalk /refs/triarchy/taxonomy_crosswalk.tsv \
  --outdir results/study
```

Primer sequences above are examples and must be replaced with those actually
used in the study.

## Execution profiles

```bash
# Local host tools
nextflow run main.nf -profile local ...

# Enable Docker for a configuration that assigns process containers
nextflow run main.nf -profile docker -c conf/containers.config ...

# Enable Singularity/Apptainer
nextflow run main.nf -profile singularity -c conf/containers.config ...

# Submit processes to Slurm with site-specific settings
nextflow run main.nf -profile slurm -c conf/slurm.config ...
```

The repository does not ship fictional universal cluster or image settings.
Add deployment-specific configuration under `conf/` and record it with the run.

## Database manifests

`conf/databases.example.yml` documents the reference-manifest shape expected by
`bin/triarchy databases check`. Production releases should replace the example
paths with the exact submitted reference locations and checksums. A manifest
entry should include `id`, `kind`, `version`, `path`, `required`, and, for
single files, `sha256`.

## Outputs

```text
results/
|-- qc/
|   `-- library_qc.tsv
|-- profiles/
|   |-- control_annotated/
|   |   |-- taxa.tsv
|   |   `-- functions.tsv
|   |-- high_confidence/taxa.tsv
|   `-- sample_confidence.tsv
|-- integrated/
|   |-- harmonized_taxa.tsv
|   |-- taxonomy_crosswalk.tsv
|   |-- taxon_evidence.tsv
|   |-- evidence_model.json
|   `-- cross_assay_concordance.tsv
|-- diversity/
|   |-- alpha_diversity.tsv
|   `-- beta_diversity.tsv
|-- pipeline_info/
|   |-- provenance.json
|   |-- trace.tsv
|   |-- timeline.html
|   |-- execution_report.html
|   `-- dag.html
`-- triarchy_report.html
```

Native assay outputs remain in the Nextflow work directory in the current
implementation. Production releases should add explicit publication rules for
the native QC, taxonomy, ASV, and HUMAnN artifacts appropriate to their storage
budget.

## Resume and troubleshooting

Resume after correcting a failed input, tool, or reference:

```bash
./run_triarchy.sh -resume [the same parameters]
```

Inspect failed tasks with:

```bash
nextflow log
nextflow log last -f process,hash,status,exit
```

Nextflow prints the failed work directory. Its `.command.sh`, `.command.err`,
and `.command.log` files contain the exact command and diagnostics.
