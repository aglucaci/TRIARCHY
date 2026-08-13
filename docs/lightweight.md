# Lightweight Local Run

Use this mode when you want TRIARCHY to run on a laptop or WSL Ubuntu install
without downloading production reference databases or creating the full
bioinformatics Conda environment.

This run exercises the full Nextflow graph with deterministic fixture outputs.
It is useful for installation checks, report inspection, documentation, and
development. It is not a biological analysis mode.

## Requirements

- Linux, macOS, or WSL Ubuntu;
- Bash;
- Java 17 or later;
- Python 3;
- `curl`.

You do not need Mamba, Conda, Kraken2, HUMAnN, DADA2, SILVA, UNITE, or the WMS
reference databases for the lightweight run.

On WSL Ubuntu:

```bash
sudo apt update
sudo apt install -y openjdk-17-jre python3 curl bash
```

## Install The Lightweight Runtime

From the repository root:

```bash
./install.sh --skip-dependencies --skip-databases --skip-tests --yes
```

This installs only the pinned local Nextflow runtime under `.tools/bin`.

## Run TRIARCHY Lite

```bash
./run_triarchy.sh --lite
```

Results are written to:

```text
results/lite
```

The main report is:

```text
results/lite/triarchy_report.html
```

On WSL, open the results folder from Ubuntu with:

```bash
explorer.exe results/lite
```

Resume an interrupted lightweight run with:

```bash
./run_triarchy.sh --lite -resume
```

## What Lite Mode Does

`--lite` is a convenience alias for:

- Nextflow `-stub-run`;
- the `lite` resource profile;
- TRI-00 validated metadata fixtures;
- `results/lite` as the default output directory.

The workflow still builds the QC, control annotation, harmonization,
cross-assay evidence, diversity, provenance, and HTML report outputs. Assay
tool steps emit deterministic mock profiles instead of running production
classifiers and reference-backed pipelines.

## Production Is Separate

For biological analysis, use the normal production install and production run
commands in `docs/running.md`. Production mode requires real FASTQs, primers,
assay tools, and reference databases.
