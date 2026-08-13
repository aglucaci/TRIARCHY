# TRIARCHY

![TRIARCHY logo](logo/triarchy_logo.png)

**TRIARCHY is a cross-assay evidence framework for matched
whole-metagenome shotgun sequencing (WMS), 16S rRNA gene sequencing, and fungal
ITS sequencing in environmental microbiome studies.**

TRIARCHY treats WMS, 16S, and ITS as complementary measurements of the same
environmental system, not as interchangeable abundance tables. It processes
each assay on its own scale, links every library to the controls that are
experimentally relevant to it, harmonizes taxonomic calls with provenance, and
then explains where the assays agree, disagree, or lack enough power to support
a conclusion.

The core question TRIARCHY answers is:

> Given matched WMS, 16S, ITS, and stage-specific control data, which
> sample-level taxa and ecological conclusions are robust, which are likely
> contaminants, and which disagreements can be explained by assay scope,
> depth, taxonomy, reference coverage, or experimental design?

## Contents

- [What This Repository Contains](#what-this-repository-contains)
- [Why TRIARCHY Exists](#why-triarchy-exists)
- [What TRIARCHY Does](#what-triarchy-does)
- [What TRIARCHY Does Not Do](#what-triarchy-does-not-do)
- [Quick Start](#quick-start)
- [System Requirements](#system-requirements)
- [Installation Modes](#installation-modes)
- [Installer Reference](#installer-reference)
- [Running The Workflow](#running-the-workflow)
- [Workflow Wrapper Reference](#workflow-wrapper-reference)
- [Production Parameter Reference](#production-parameter-reference)
- [Input Data Model](#input-data-model)
- [Workflow Overview](#workflow-overview)
- [Evidence Model](#evidence-model)
- [Output Files](#output-files)
- [Output Table Schemas](#output-table-schemas)
- [Inspecting Results](#inspecting-results)
- [Benchmark Suite](#benchmark-suite)
- [Real-World Validation Downloads](#real-world-validation-downloads)
- [Repository Layout](#repository-layout)
- [Command Line Helpers](#command-line-helpers)
- [Configuration And References](#configuration-and-references)
- [Database Manifest Format](#database-manifest-format)
- [Reference Installation Details](#reference-installation-details)
- [Development And Testing](#development-and-testing)
- [Troubleshooting](#troubleshooting)
- [Glossary](#glossary)
- [Project Status](#project-status)
- [License](#license)

## What This Repository Contains

This repository contains an executable Nextflow workflow plus the scripts,
configuration, documentation, and benchmark fixtures needed to develop and test
TRIARCHY.

The repo includes:

- a Nextflow workflow entrypoint: `main.nf`;
- local Nextflow and environment installation helpers: `install.sh`;
- a workflow wrapper for lite, stub, and production runs: `run_triarchy.sh`;
- a combined setup-and-run entrypoint: `setup_and_run_triarchy.sh`;
- process modules for WMS, 16S, ITS, metadata validation, and integration;
- Python and R scripts for metadata validation, evidence integration,
  reporting, and benchmark generation;
- configuration files for assay performance, evidence policy, references, and
  execution profiles;
- deterministic benchmark fixtures and truth tables;
- documentation for running, data contracts, methods, results, and release
  readiness.

TRIARCHY is designed to run in a POSIX environment: Linux, macOS, or WSL
Ubuntu on Windows.

## Why TRIARCHY Exists

Environmental microbiome studies often combine multiple sequencing assays:

- **WMS** can observe bacteria, archaea, fungi, viruses, and functional
  potential, but sensitivity depends heavily on depth, biomass, host content,
  and reference representation.
- **16S** preferentially observes bacterial and archaeal marker genes, with
  primer bias and limited taxonomic resolution.
- **ITS** preferentially observes fungi, with its own primer, copy-number, and
  taxonomic constraints.

A common failure mode is to collapse these assays into one merged abundance
table or to treat a missing call in one assay as a biological absence. That can
misclassify contaminants, overstate disagreement, or hide useful evidence from
an assay that is simply outside scope or underpowered.

TRIARCHY instead keeps the assay outputs separate and integrates evidence only
after recording:

- which assays were actually measured for each biological sample;
- whether each assay had enough depth and quality to detect the taxon;
- whether the assay is biologically appropriate for the kingdom being tested;
- whether reference coverage supports a defensible negative call;
- whether linked negative controls suggest contamination;
- whether disagreement is only a taxonomic-rank mismatch.

The result is an auditable evidence layer, not a forced multi-assay abundance
matrix.

## What TRIARCHY Does

TRIARCHY implements a multi-stage workflow:

1. **Validate the experimental design.**
   It reads a sample manifest and a sample-control edge table, then validates
   assay labels, sample roles, required batches, FASTQ paths, paired controls,
   and cross-assay consistency.

2. **Build a sample-control graph.**
   Controls are stage-linked: field, extraction, PCR, library, and sequencing
   controls remain connected to the biological libraries they are relevant to.

3. **Process each assay independently.**
   WMS, 16S, and ITS are processed with assay-appropriate tools and retained
   on assay-native scales.

4. **Annotate contamination without overwriting raw values.**
   TRIARCHY emits unfiltered, control-annotated, and high-confidence profiles.
   The linked-control evidence remains inspectable.

5. **Harmonize taxonomy with provenance.**
   Source names and identifiers are mapped to canonical taxon records while
   retaining original database identifiers, names, ranks, mapping methods, and
   confidence.

6. **Estimate authenticity and explain cross-assay evidence.**
   Bacterial and archaeal calls compare WMS with 16S. Fungal calls compare WMS
   with ITS. TRIARCHY first estimates biological authenticity, then assigns an
   explanation for concordance, assay limitation, taxonomic limitation,
   contamination, discordance, or indeterminate evidence.

7. **Calculate assay-specific ecology.**
   Alpha and beta diversity are calculated within assay scope. Cross-assay
   concordance compares detection, rank order, and ecological conclusions
   rather than raw relative abundance.

8. **Generate an analytical report.**
   The HTML report summarizes QC, sample suitability, taxon-level evidence,
   control evidence, cross-assay concordance, diversity, and provenance.

## What TRIARCHY Does Not Do

TRIARCHY is intentionally narrower than a universal microbiome platform.

It does not:

- claim to be the first WMS, 16S, or ITS processing workflow;
- average WMS, 16S, and ITS relative abundances into a single abundance table;
- interpret a missing assay as a measured negative result;
- treat an underpowered assay as evidence of absence;
- treat functional WMS pathway evidence as if it were independently validated
  by 16S or ITS;
- provide clinically calibrated diagnostic probabilities;
- replace physical mock communities or real environmental validation data.

Existing classifiers, marker-gene reconciliation methods, per-assay
decontamination tools, global-control approaches, and multi-table methods are
comparators or upstream evidence sources. TRIARCHY's contribution is the
stage-linked, control-aware integration layer that explains agreement and
disagreement across assays.

## Quick Start

### Requirements

For the lightweight local mode:

- Linux, macOS, or WSL Ubuntu;
- Bash;
- Java 17 or later;
- Python 3;
- `curl`.

For production analysis, TRIARCHY also needs Mamba, Micromamba, or Conda plus
the production bioinformatics tools and references installed by `install.sh`.

### One-Command Setup And Run

Use the combined wrapper when you want a single command to install what the
selected mode needs and then launch the workflow:

```bash
cd /path/to/TRIARCHY
bash setup_and_run_triarchy.sh --demo
```

For a real production run, provide the production inputs after `--`:

```bash
bash setup_and_run_triarchy.sh --production \
  --database-dir /refs/triarchy \
  --database-profile standard \
  --cpus 8 \
  --yes \
  -- \
  -profile local \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --outdir results/study
```

The standard production database profile can exceed 400 GB and take many
hours. Completed database stages are checkpointed and reused on reruns.

### Lightest Local Run

Use this mode first on a laptop or WSL Ubuntu machine. It installs only the
local Nextflow runtime and runs the complete workflow graph with deterministic
stub outputs.

```bash
cd /path/to/TRIARCHY
./install.sh --skip-dependencies --skip-databases --skip-tests --yes
./run_triarchy.sh --lite
```

Lite results are written to:

```text
results/lite
```

The main report is:

```text
results/lite/triarchy_report.html
```

Lite mode is useful for installation checks, documentation, and development.
It is not suitable for biological analysis because it does not run the assay
classifiers or reference-backed pipelines. See [docs/lightweight.md](docs/lightweight.md).

### Showcase Demo Run

Use this mode to run the TRI-31 showcase dataset, which implements the
30-minute demo design with 12 biological samples, seven controls, WMS, 16S,
ITS, host contamination, reagent contamination, batch contamination, primer
bias, rank ambiguity, fungal reference limitation, expected tables, expected
figures, mini function references, and a precomputed narrative report.

```bash
./run_triarchy.sh --demo
```

Demo results are written to:

```text
results/showcase_30min
```

The dataset itself lives under:

```text
benchmarks/TRI-31
```

The committed FASTQs are scaled down for Git and CI. The dataset metadata
records the full-scale target of 2,526,000 paired-end reads for a future
publication-sized bundle.

### Stub Run With Minimal Install

The deterministic stub workflow can also be launched after a minimal database
profile installation:

```bash
cd /path/to/TRIARCHY
./install.sh --database-profile minimal --yes
./run_triarchy.sh --stub
```

Stub results are written to:

```text
results/stub
```

Resume an interrupted stub or lite run with `-resume`:

```bash
./run_triarchy.sh --demo -resume
./run_triarchy.sh --lite -resume
./run_triarchy.sh --stub -resume
```

### Windows And WSL Note

Run TRIARCHY from Ubuntu/WSL, not from PowerShell. For better Nextflow
performance, copy the repository into the Ubuntu filesystem, for example
`~/TRIARCHY`, rather than running large workflows directly under `/mnt/c`.

If Ubuntu reports `/usr/bin/env: 'bash\r': No such file or directory`, the
script has Windows CRLF line endings. This repository includes `.gitattributes`
to keep shell scripts LF-normalized; rerun Git checkout or convert the scripts
with `dos2unix` if needed.

## System Requirements

### Lightweight And Stub Requirements

The lightweight and stub workflows run the full TRIARCHY orchestration graph
with deterministic mock assay outputs. They require only:

| Requirement | Notes |
| --- | --- |
| POSIX shell | Linux, macOS, or WSL Ubuntu |
| Bash | Used by `install.sh`, `run_triarchy.sh`, and helper scripts |
| Java 17 or later | Required by Nextflow |
| Python 3 | Required by validation and integration scripts |
| `curl` | Required to install the pinned local Nextflow runtime |

On WSL Ubuntu:

```bash
sudo apt update
sudo apt install -y openjdk-17-jre python3 curl bash
```

### Production Requirements

Production runs additionally require the tools and references installed by the
full installer:

| Component | Purpose |
| --- | --- |
| Mamba, Micromamba, or Conda | Creates the local production environment |
| fastp | WMS read QC |
| Bowtie2 | Host read removal |
| Kraken2 and Bracken | WMS taxonomic profiling |
| HUMAnN, MetaPhlAn, DIAMOND | WMS functional profiling |
| cutadapt | 16S and ITS primer/read trimming |
| ITSxpress | ITS region extraction |
| R and DADA2 | 16S and ITS ASV inference and taxonomy |
| Python 3 packages from `environment.yml` | Integration, reporting, utilities |
| Host, Kraken2/Bracken, HUMAnN, SILVA, and ITS references | Production biological analysis |

The standard reference set can exceed 400 GB. Plan for additional space for
FASTQs, Nextflow `work/`, intermediate files, and final results. On WSL, a
native Ubuntu filesystem path is usually faster and less error-prone than
`/mnt/c` for large production jobs.

## Installation Modes

### Lightweight Runtime Only

This installs only local Nextflow under `.tools/bin`.

```bash
./install.sh --skip-dependencies --skip-databases --skip-tests --yes
```

Use this with:

```bash
./run_triarchy.sh --lite
```

### Minimal Installation Check

This creates the Conda/Mamba environment and installs the minimal reference
profile used for installation checks and deterministic stub runs.

```bash
./install.sh --database-profile minimal --yes
```

The minimal profile installs a tiny host Bowtie2 index, a tiny Kraken2/Bracken
test database, HUMAnN demo databases, and marker taxonomy files. It is not
suitable for biological interpretation.

### Production Installation

For a complete local production installation, including standard references:

```bash
./install.sh --yes
```

The production installer:

- installs the pinned Nextflow runtime under `.tools/bin`;
- creates `.tools/env` with WMS, 16S, ITS, R, and Python dependencies;
- installs fastp, Bowtie2, Kraken2, Bracken, HUMAnN, MetaPhlAn, DIAMOND,
  cutadapt, ITSxpress, R, and DADA2 through the configured environment;
- downloads and validates WMS, host, HUMAnN, SILVA, and ITS references;
- generates `conf/references.config`;
- records package and database provenance;
- validates the TRI-00 metadata fixture and runs repository tests.

The standard reference installation can exceed 400 GB and take many hours.
Database installation is checkpointed by validated stage, so rerunning the same
command resumes completed work where possible.

Use a larger reference destination when needed:

```bash
./install.sh \
  --database-profile standard \
  --database-dir /path/to/large/reference/storage \
  --cpus 8 \
  --yes
```

Use a separately downloaded Kraken archive with:

```bash
./install.sh \
  --database-profile standard \
  --kraken-archive /path/to/k2_pluspfp_20260626.tar.gz \
  --yes
```

## Installer Reference

`install.sh` sets up local runtime pieces. It is safe to rerun the same install
command after an interruption; completed reference stages are validated and
skipped where possible.

Usage:

```bash
./install.sh [options]
```

Common options:

| Option | Meaning |
| --- | --- |
| `--prefix DIR` | Install the local Nextflow executable in `DIR`; default `.tools/bin` |
| `--env-prefix DIR` | Install the Conda/Mamba environment in `DIR`; default `.tools/env` |
| `--environment-manager CMD` | Use `mamba`, `micromamba`, `conda`, or `auto` |
| `--database-dir DIR` | Reference destination; default `references/data` or last successful destination |
| `--database-profile minimal` | Install tiny references for checks and stub-compatible setup |
| `--database-profile standard` | Install production references |
| `--kraken-archive FILE` | Use an existing PlusPFP archive rather than downloading it |
| `--cpus N` | Threads for local indexing/build phases |
| `--skip-nextflow` | Do not install or validate local Nextflow |
| `--skip-dependencies` | Do not create or update the Conda/Mamba environment |
| `--skip-databases` | Do not install references |
| `--skip-tests` | Do not run metadata validation and unit tests |
| `--yes` | Accept non-interactive prompts, including large standard references |
| `-h`, `--help` | Print installer help |

Environment variables:

| Variable | Meaning |
| --- | --- |
| `NXF_VER` | Override the version pinned in `.nextflow-version` |
| `NXF_HOME` | Override Nextflow cache location |
| `TRIARCHY_ENV_MANAGER` | Default environment manager for install |
| `TRIARCHY_DOWNLOAD_ATTEMPTS` | Override database download retry attempts |
| `TRIARCHY_MINIMAL_HOST_ACCESSION` | Override the minimal-profile host FASTA accession label |
| `TRIARCHY_MINIMAL_HOST_URL` | Override the minimal-profile host FASTA URL |
| `TRIARCHY_MINIMAL_KRAKEN_MAX_DB_SIZE` | Override the minimal Kraken2 hash-size cap |
| `NO_COLOR=1` | Disable colored installer output |

Recommended install commands:

```bash
# Laptop/WSL smoke test only
./install.sh --skip-dependencies --skip-databases --skip-tests --yes

# Minimal local check
./install.sh --database-profile minimal --yes

# Full production install
./install.sh --database-profile standard --database-dir /refs/triarchy --cpus 8 --yes
```

## Running The Workflow

The wrapper script configures `PATH`, `NXF_HOME`, and `NXF_VER`, then runs
Nextflow against `main.nf`.

### Lite Mode

```bash
./run_triarchy.sh --lite
```

Equivalent behavior:

- enables Nextflow `-stub-run`;
- uses the `lite` Nextflow profile;
- uses the TRI-00 validated metadata fixture;
- writes to `results/lite` by default.

### Stub Mode

```bash
./run_triarchy.sh --stub
```

Stub mode verifies the complete graph with deterministic outputs. The TRI-00
FASTQ paths are intentionally synthetic, so this fixture must use Nextflow
stub execution.

### Production Mode

Production mode runs without `-stub-run` and requires real FASTQs, primers,
assay tools, taxonomy references, and WMS databases.

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

Primer sequences above are examples. Use the primers actually used in the
study.

See [docs/running.md](docs/running.md) for production parameters, execution
profiles, restart behavior, and output inspection.

## Workflow Wrapper Reference

`run_triarchy.sh` is the preferred entrypoint for local use. It configures the
local Nextflow installation and forwards unrecognized options directly to
Nextflow.

Usage:

```bash
./run_triarchy.sh --lite [NEXTFLOW_OPTIONS]
./run_triarchy.sh --stub [NEXTFLOW_OPTIONS]
./run_triarchy.sh [NEXTFLOW_OPTIONS] --samples FILE --controls FILE
```

Wrapper options:

| Option | Meaning |
| --- | --- |
| `--demo` | Run the TRI-31 showcase deterministic workflow with `-stub-run`, `-profile demo`, and `results/showcase_30min` |
| `--lite` | Run a laptop-friendly deterministic workflow using `-stub-run`, `-profile lite`, TRI-00 metadata, and `results/lite` |
| `--stub` | Run the complete deterministic stub workflow using TRI-00 metadata and `results/stub` |
| `-h`, `--help` | Print wrapper help |
| `--` | Pass all following arguments directly to Nextflow |

Environment variables accepted by the wrapper:

| Variable | Meaning |
| --- | --- |
| `TRIARCHY_SAMPLES` | Default production sample manifest |
| `TRIARCHY_CONTROLS` | Default production control edge table |
| `TRIARCHY_OUTDIR` | Default output directory |
| `NXF_HOME` | Nextflow cache; default `.nextflow` under the project |
| `NXF_VER` | Nextflow version; default from `.nextflow-version` |

Examples:

```bash
# Lightweight laptop check
./run_triarchy.sh --lite

# Showcase demo dataset
./run_triarchy.sh --demo

# Reuse successful tasks after interruption
./run_triarchy.sh --lite -resume

# Stub a benchmark dataset through the graph
./run_triarchy.sh --stub \
  --samples benchmarks/TRI-12/metadata/samples.tsv \
  --controls benchmarks/TRI-12/metadata/control_relationships.tsv \
  --outdir results/benchmarks/TRI-12

# Production local run
./run_triarchy.sh -profile local \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --outdir results/study
```

## Production Parameter Reference

Required production parameters:

| Parameter | Required | Meaning |
| --- | --- | --- |
| `--samples FILE` | Yes | Sample manifest TSV |
| `--controls FILE` | Yes | Sample-control edge TSV |
| `--host_index PREFIX` | Yes for WMS | Bowtie2 host index prefix |
| `--kraken2_db DIR` | Yes for WMS | Kraken2/Bracken database directory |
| `--humann_db DIR` | Yes for HUMAnN | HUMAnN ChocoPhlAn nucleotide database |
| `--humann_protein_db DIR` | Yes for HUMAnN | HUMAnN UniRef DIAMOND database |
| `--s16_taxonomy FILE` | Yes for 16S | DADA2-compatible bacterial/archaeal taxonomy FASTA |
| `--its_taxonomy FILE` | Yes for ITS | DADA2-compatible fungal taxonomy FASTA |
| `--cutadapt_16s_forward SEQ` | Yes for 16S | Forward 16S primer used by the study |
| `--cutadapt_16s_reverse SEQ` | For paired 16S | Reverse 16S primer used by the study |
| `--cutadapt_its_forward SEQ` | Yes for ITS | Forward ITS primer used by the study |
| `--cutadapt_its_reverse SEQ` | For paired ITS | Reverse ITS primer used by the study |
| `--outdir DIR` | Recommended | Output directory; default `results` |

Optional and advanced parameters:

| Parameter | Meaning |
| --- | --- |
| `--taxonomy_crosswalk FILE` | Explicit source-to-canonical taxonomy crosswalk |
| `--assay_performance FILE` | Override assay reliability, depth, coverage, and false-positive defaults |
| `--evidence_policy FILE` | Override thresholds, prior, model version, and calibration metadata |
| `--fastp_extra STRING` | Extra arguments passed to fastp |
| `--bracken_read_length N` | Bracken read length; default `150` |

Path behavior:

- relative FASTQ paths in the sample manifest are resolved relative to the
  sample manifest file;
- reference paths are resolved before processes enter isolated Nextflow work
  directories;
- production runs should use stable absolute paths for large references.

Preflight before a large production run:

```bash
bin/triarchy validate \
  --input /data/study/samples.tsv \
  --paired-controls /data/study/control_relationships.tsv \
  --databases conf/databases.example.yml \
  --database-destination /refs/triarchy \
  --report results/preflight.json
```

## Input Data Model

TRIARCHY uses two required tab-separated metadata inputs.

### Sample Manifest

The sample manifest has one row per sequenced library. A physical specimen
assayed by WMS, 16S, and ITS therefore appears as three rows sharing the same
`biological_sample_id`.

Sample manifest columns:

| Field | Required | Meaning |
| --- | --- | --- |
| `library_id` | Yes | Unique sequenced-library identifier |
| `biological_sample_id` | Yes | Identifier shared by matched assay libraries |
| `specimen_id` | Yes | Physical specimen identifier |
| `subject_id` | Yes | Source subject/site identifier; controls should use a documented control value |
| `assay` | Yes | `WMS`, `16S`, or `ITS` |
| `sample_role` | Yes | `biological`, `negative_control`, `positive_control`, or `mock` |
| `sample_material` | Yes | Material or substrate sampled |
| `environment_category` | Yes | Airstream, surface, dust, water, soil, or another prespecified environment label |
| `library_layout` | Yes | `paired` or `single` |
| `fastq_r1` | Yes | Forward read path |
| `fastq_r2` | For paired libraries | Reverse read path |
| `extraction_batch` | Yes | Extraction batch identifier |
| `pcr_plate` | For amplicons | PCR plate identifier |
| `library_batch` | Yes | Library-preparation batch identifier |
| `sequencing_run` | Yes | Sequencing run identifier |
| `technical_replicate` | Yes | Positive integer within specimen and assay |
| `biomass_measure` | No | Quantitative biomass proxy |
| `biomass_unit` | If measured | Unit for `biomass_measure` |

Rules:

- identifiers must be stable strings;
- assay values are case-sensitive;
- paths must exist for production validation;
- paired libraries must provide `fastq_r2`;
- biological samples matched across assays must not conflict in specimen,
  subject, material, environment, or batch metadata that should be shared.

### Sample-Control Edge Table

Controls are represented as graph edges rather than a single generic control
column. One biological library can have multiple linked controls.

Control edge columns:

| Field | Required | Meaning |
| --- | --- | --- |
| `library_id` | Yes | Biological library |
| `control_library_id` | Yes | Linked control library |
| `control_stage` | Yes | `field`, `extraction`, `pcr`, `library`, or `sequencing` |
| `relationship` | Yes | `paired`, `same_batch`, or `run_level` |
| `priority` | Yes | Positive integer; `1` is the closest control |
| `reason` | Yes | Human-readable experimental justification |

Both identifiers must exist in the manifest, must use the same assay, and the
target must have a control sample role.

Validation rejects:

- duplicate library identifiers;
- unknown assay labels or sample roles;
- missing required batch fields;
- control edges whose identifiers do not exist;
- cross-assay control edges;
- biological samples used as negative controls;
- paired assays whose biological sample metadata conflict.

### Example Minimal Metadata Shape

Example sample rows:

```text
library_id	biological_sample_id	specimen_id	subject_id	assay	sample_role	sample_material	environment_category	library_layout	fastq_r1	fastq_r2	extraction_batch	pcr_plate	library_batch	sequencing_run	technical_replicate	biomass_measure	biomass_unit
S01_WMS	S01	SPEC01	SITE01	WMS	biological	air_filter	return_air	paired	reads/S01_WMS_R1.fastq.gz	reads/S01_WMS_R2.fastq.gz	EXT01		LIB01	RUN01	1	12.4	ng
S01_16S	S01	SPEC01	SITE01	16S	biological	air_filter	return_air	paired	reads/S01_16S_R1.fastq.gz	reads/S01_16S_R2.fastq.gz	EXT01	PCR01	LIB01	RUN01	1	12.4	ng
BLANK01_WMS	BLANK01	CTRL01	control	WMS	negative_control	extraction_blank	control	paired	reads/BLANK01_WMS_R1.fastq.gz	reads/BLANK01_WMS_R2.fastq.gz	EXT01		LIB01	RUN01	1	0	ng
```

Example control edge rows:

```text
library_id	control_library_id	control_stage	relationship	priority	reason
S01_WMS	BLANK01_WMS	extraction	same_batch	1	Same extraction batch blank
S01_16S	BLANK01_16S	extraction	same_batch	1	Same extraction batch blank
```

See [docs/data-contracts.md](docs/data-contracts.md) for the full machine-readable
schema and output table contracts.

## Workflow Overview

TRIARCHY is implemented as a Nextflow DSL2 workflow.

```text
metadata + controls
  -> validated sample-control graph
  -> assay-specific processing
       WMS: fastp -> Bowtie2 -> Kraken2/Bracken -> HUMAnN
       16S: cutadapt -> DADA2 -> bacterial/archaeal taxonomy
       ITS: cutadapt -> ITSxpress -> DADA2 -> fungal taxonomy
  -> linked-control annotation
  -> taxonomic harmonization
  -> TRIARCHY authenticity and cross-assay explanation model
  -> assay-specific diversity and cross-assay concordance
  -> HTML report and provenance
```

The process modules live under `modules/local/`:

- `metadata.nf`: manifest validation and assay split;
- `wms.nf`: WMS QC, host removal, taxonomy, and functions;
- `amplicon.nf`: 16S/ITS trimming, ITS extraction, and DADA2 profiles;
- `integration.nf`: QC summaries, contamination annotation, harmonization,
  authenticity/explanation inference, diversity, provenance, and report generation.

The workflow graph is documented in [flowchart.md](flowchart.md). The canonical
Mermaid source is `docs/overall-workflow-dag.mmd`.

## Evidence Model

TRIARCHY's integration layer compares evidence only where assay scope overlaps:

- bacteria and archaea: WMS <-> 16S;
- fungi: WMS <-> ITS;
- functional potential: WMS only.

For each biological sample and harmonized taxon, the evidence model tracks
whether each assay was:

- measured and detected;
- measured but not detected;
- not measured;
- technically underpowered;
- reference-limited;
- contaminated or control-associated.

Detection power is estimated from assay reliability, library depth, quality
state, kingdom, primer compatibility when supplied, and reference coverage. A
measured nondetection from a high-power assay is interpreted differently from a
missing assay or a low-depth library.

The evidence model is two-stage. Stage I emits `biological_authenticity`, a
posterior-style estimate that the taxon is biologically authentic in the
sample, with approximate uncertainty bounds (`authenticity_ci_low`,
`authenticity_ci_high`) and a conservative `authenticity_state`. Stage II emits
`cross_assay_status`, a parent explanation class, and the detailed
`classification` used for audit trails and compatibility with earlier outputs.
Wide uncertainty intervals are reported as `indeterminate` rather than forced
into present/absent calls.

The default model version is serialized in `integrated/evidence_model.json`.
Configuration lives in:

- `conf/assay_performance.tsv`;
- `conf/evidence_policy.json`.

The default probability mapping is labeled
`synthetic_only_not_empirically_calibrated`. It is useful for transparent
ranking and development, but it must not be presented as clinical or empirical
environmental calibration.

### Evidence Signals

For each sample, taxon, and assay, TRIARCHY combines:

- detection support;
- within-assay abundance/rank support;
- taxonomic resolution support;
- mapping confidence;
- assay reliability;
- reference coverage;
- library depth;
- library quality;
- linked-control contamination evidence.

Per-assay positive evidence is represented conceptually as:

```text
e(s,t,m) = detection_support
           * sqrt(abundance_support)
           * resolution_support
           * mapping_support
           * reliability
           * (1 - contamination_probability)
```

Detection power for a measured assay is represented conceptually as:

```text
power(s,t,m) = reliability(m,t)
               * [1 - exp(-library_depth / minimum_depth)]
               * reference_coverage(m,t)
               * quality_factor(s,m)
```

These formulas are used to keep the reasoning explicit. The default parameters
are configuration, not universal biological constants.

### Control Evidence

For each linked control relationship, TRIARCHY evaluates:

1. **Prevalence:** whether a feature appears more often in linked controls than
   in biological samples.
2. **Abundance:** whether a feature is enriched in the closest blank or
   relevant batch controls.
3. **Batch localization:** whether detections cluster in a specific extraction,
   PCR, library, or sequencing batch.

The resulting contamination score is a bounded evidence score, not a calibrated
probability of contamination. The raw profile, annotated profile, and
high-confidence profile are all retained.

### Taxonomic Rank Handling

TRIARCHY does not merge by display name alone. It compares taxa at the lowest
rank valid in both assays. For example, if WMS reports a species and 16S
supports only the genus, the comparison is evaluated at genus and classified as
rank-limited rather than species-level concordance.

### Calibration

The shipped calibration mapping is synthetic-only. To fit a study-specific
calibration against truth data:

```bash
python3 workflow/scripts/triarchy_integrate.py fit-calibration \
  --evidence results/integrated/taxon_evidence.tsv \
  --truth benchmark_truth.tsv \
  --output fitted_calibration.json
```

Use a fitted calibration only within the scope of the truth data and assay
conditions used to create it.

### Authenticity And Explanation Classes

`cross_assay_status` is the parent interpretation used for summaries:

| Status | Interpretation |
| --- | --- |
| `concordant` | Kingdom-appropriate assays support the same biological signal |
| `assay_limited` | Apparent disagreement is explained by assay scope, depth, quality, or missing modality |
| `taxonomy_limited` | Agreement is limited by taxonomic rank, reference representation, or compatibility |
| `contamination_associated` | Linked controls or batch localization explain the signal |
| `genuine_discordance` | Technically adequate assays provide incompatible evidence |
| `indeterminate` | Evidence or uncertainty is insufficient for a stronger explanation |

The detailed `classification` column retains the software-level explanation:

| Class | Interpretation |
| --- | --- |
| `concordant_detection` | All kingdom-appropriate measured assays support the taxon |
| `single_assay_credible` | One appropriate assay supports the taxon and no adequate assay contradicts it |
| `assay_limited_nondetection` | A nondetecting assay lacked enough power for a strong negative call |
| `rank_limited_agreement` | Assays agree only at a shared higher taxonomic rank |
| `reference_limited` | Reference representation is insufficient for a defensible negative call |
| `probable_contaminant` | Linked controls support contamination |
| `cross_contamination_candidate` | Control-associated evidence is localized to a batch or run |
| `genuine_discordance` | Technically adequate assays provide incompatible evidence |
| `indeterminate` | Depth, controls, references, or resolution are insufficient |

These labels are evidence summaries, not declarations of absolute biological
truth.

## Output Files

A completed run writes a structured output directory under `params.outdir`.

Important files:

| Path | Purpose |
| --- | --- |
| `triarchy_report.html` | Main analytical report |
| `qc/library_qc.tsv` | Library-level assay QC and depth summary |
| `profiles/unfiltered/taxa.tsv` | Assay-native taxonomic profiles before filtering |
| `profiles/control_annotated/taxa.tsv` | Profiles with linked-control contamination evidence |
| `profiles/high_confidence/taxa.tsv` | Taxa passing the configured evidence policy |
| `profiles/sample_confidence.tsv` | Biological sample suitability and failure reasons |
| `integrated/harmonized_taxa.tsv` | Harmonized taxon table |
| `integrated/taxonomy_crosswalk.tsv` | Source-to-canonical taxon mapping provenance |
| `integrated/taxon_evidence.tsv` | Sample-taxon authenticity estimates and cross-assay explanations |
| `integrated/evidence_model.json` | Serialized model, policy, and calibration metadata |
| `integrated/cross_assay_concordance.tsv` | Kingdom-aware cross-assay agreement metrics |
| `diversity/alpha_diversity.tsv` | Within-sample diversity by assay |
| `diversity/beta_diversity.tsv` | Between-sample distances by assay |
| `pipeline_info/provenance.json` | Run provenance |
| `pipeline_info/timeline.html` | Nextflow execution timeline |
| `pipeline_info/execution_report.html` | Nextflow resource report |
| `pipeline_info/trace.tsv` | Task-level execution trace |
| `pipeline_info/dag.html` | Workflow DAG for the run |

TRIARCHY retains assay-native counts and normalized values. Cross-assay tables
contain detection evidence, confidence, contamination status, concordance, and
rank-aware summaries, not averaged relative abundance.

See [docs/examining-results.md](docs/examining-results.md) for practical result
inspection commands.

## Output Table Schemas

### Control-Annotated Feature Profile

`profiles/control_annotated/taxa.tsv` and the related function table preserve
assay-native measurements and add linked-control evidence.

| Field | Meaning |
| --- | --- |
| `library_id` | Sequenced library |
| `assay` | `WMS`, `16S`, or `ITS` |
| `triarchy_taxon_id` | Harmonized feature identifier |
| `native_count` | Original assay count |
| `native_normalized_value` | Within-assay normalized value |
| `normalization_method` | Method used within this assay |
| `control_prevalence` | Prevalence among linked negative controls |
| `biological_prevalence` | Prevalence among relevant biological samples |
| `control_abundance_ratio` | Linked-control to biological abundance evidence |
| `batch_localization_score` | Concentration of signal in a laboratory batch |
| `contamination_probability` | Bounded linked-control evidence score from 0 to 1 |
| `associated_control_types` | Semicolon-separated control stages |
| `associated_batches` | Semicolon-separated batch identifiers |
| `decision` | `retain`, `flag`, `exclude`, or `indeterminate` |
| `decision_reason` | Human-readable rule/model explanation |

### Taxonomy Crosswalk

`integrated/taxonomy_crosswalk.tsv` records how source assignments map into a
canonical representation.

| Field | Meaning |
| --- | --- |
| `triarchy_taxon_id` | Canonical TRIARCHY taxon identifier |
| `canonical_name` | Canonical taxon name |
| `canonical_rank` | Canonical rank used by TRIARCHY |
| `canonical_lineage` | Canonical lineage string |
| `source_database` | Source database or classifier |
| `source_taxon_id` | Source taxon identifier |
| `original_name` | Source-reported name |
| `original_rank` | Source-reported rank |
| `mapping_method` | Exact, synonym, lineage, manual, or other mapping method |
| `mapping_confidence` | Numeric confidence from 0 to 1 |
| `database_version` | Source database version |

### Integrated Taxon Evidence

`integrated/taxon_evidence.tsv` has one row per biological sample and
harmonized taxon.

| Field | Meaning |
| --- | --- |
| `biological_sample_id` | Matched biological sample identifier |
| `triarchy_taxon_id` | Harmonized taxon |
| `comparison_rank` | Lowest valid rank shared by compared assays |
| `*_measured` | Whether each assay was performed for the biological sample |
| `*_detected` | Whether each assay produced a compatible taxon call |
| `*_library_depth` | Assay-native total used for power estimation |
| `*_detection_power` | Estimated detection probability under the configured model |
| `*_reference_coverage` | Configured or taxon-specific reference representation |
| `wms_evidence_score` | WMS support from detection, abundance, resolution, mapping, and reliability |
| `16s_evidence_score` | 16S support from detection, abundance, resolution, mapping, and reliability |
| `its_evidence_score` | ITS support from detection, abundance, resolution, mapping, and reliability |
| `abundance_support` | Rank-based within-assay abundance support |
| `taxonomic_resolution_support` | Support implied by the shared comparison rank |
| `mapping_support` | Minimum mapping confidence among supporting calls |
| `assay_reliability` | Kingdom-specific assay weights used |
| `control_evidence` | Summary contamination evidence |
| `field_contamination_score` | False-detection evidence linked to field controls |
| `extraction_contamination_score` | False-detection evidence linked to extraction controls |
| `pcr_contamination_score` | False-detection evidence linked to PCR controls |
| `library_contamination_score` | False-detection evidence linked to library controls |
| `run_contamination_score` | False-detection evidence linked to sequencing/run controls |
| `other_control_contamination_score` | False-detection evidence linked to other declared controls |
| `probable_contamination_source` | Stage or batch with the strongest contamination evidence |
| `taxonomic_compatibility_score` | Compatibility score across kingdom-appropriate assay calls |
| `taxonomic_compatibility_status` | Exact/canonical, hierarchical, single-assay, ambiguous, or incompatible compatibility state |
| `biological_authenticity` | Stage-I posterior-style estimate that the taxon is biologically authentic |
| `authenticity_ci_low` | Lower approximate uncertainty bound for `biological_authenticity` |
| `authenticity_ci_high` | Upper approximate uncertainty bound for `biological_authenticity` |
| `authenticity_uncertainty` | Width of the approximate authenticity interval |
| `authenticity_state` | `biologically_supported`, `not_supported`, or `indeterminate` |
| `cross_assay_status` | Stage-II parent explanation class |
| `presence_score` | Legacy alias mirroring `biological_authenticity` for older consumers |
| `probability_authentic` | Legacy calibrated-score column retained for compatibility |
| `calibration_status` | Calibration scope, including explicit non-empirical status |
| `biological_state` | `present`, `absent`, or `indeterminate` |
| `sample_evidence_rank` | Descending score rank within a biological sample |
| `global_evidence_rank` | Descending score rank across the run |
| `classification` | Detailed TRIARCHY evidence class retained for compatibility and audit trails |
| `discordance_reason` | Resolution, scope, control, direction, or depth explanation |
| `model_version` | Version of the decision policy or model |

### Cross-Assay Concordance

`integrated/cross_assay_concordance.tsv` contains one sample, kingdom, and
biologically valid assay pair per row.

Typical fields include:

- biological sample identifier;
- kingdom;
- assay pair;
- assay-specific taxon counts;
- shared taxon count;
- union taxon count;
- presence Jaccard;
- Spearman rank correlation on shared taxa;
- rank agreement score.

### Sample Confidence

`profiles/sample_confidence.tsv` summarizes whether a biological sample is
suitable for downstream use:

```text
biological_sample_id
wms_qc_status
16s_qc_status
its_qc_status
control_status
cross_assay_status
recommended_use
failure_reasons
```

### Run Provenance

`pipeline_info/provenance.json` records:

- workflow version and Git commit when available;
- command and resolved configuration;
- software versions and container digests when available;
- database names, releases, paths, checksums, and retrieval dates;
- checksums for manifests and control edges;
- run start and completion timestamps;
- task-level CPU, memory, storage, and runtime information from Nextflow.

The provenance record should be sufficient to identify the inputs and
environment. It must not embed credentials or protected sample metadata.

## Inspecting Results

Start with:

```text
<OUTDIR>/triarchy_report.html
<OUTDIR>/pipeline_info/execution_report.html
<OUTDIR>/pipeline_info/dag.html
```

On WSL, open a results folder in Windows Explorer:

```bash
explorer.exe results/lite
```

On Linux:

```bash
xdg-open results/lite/triarchy_report.html
```

On macOS:

```bash
open results/lite/triarchy_report.html
```

Inspect a TSV in the terminal:

```bash
column -t -s $'\t' results/lite/integrated/taxon_evidence.tsv | less -S
```

List failed or incomplete tasks:

```bash
awk -F '\t' 'NR == 1 || $5 != "COMPLETED"' \
  results/lite/pipeline_info/trace.tsv | column -t -s $'\t'
```

Compare a benchmark truth table with observed evidence:

```bash
column -t -s $'\t' benchmarks/TRI-12/truth/expected_taxon_status.tsv | less -S
column -t -s $'\t' results/benchmarks/TRI-12/integrated/taxon_evidence.tsv | less -S
```

## Benchmark Suite

The benchmark suite exercises metadata validation, assay branches,
contamination modeling, harmonization, authenticity/explanation inference, and
ecological recovery.

Dataset tiers include:

| Tier | IDs | Purpose |
| --- | --- | --- |
| Software/CI | TRI-00 to TRI-02 | Contracts, smoke tests, and isolated assay branches |
| Compact bundle | TRI-01, TRI-03, TRI-05 to TRI-07 | Small multi-dataset resource profile |
| Synthetic truth | TRI-10 to TRI-19 | Accuracy, sensitivity, contamination, ecology, harmonization, failures, carryover, repeatability |
| Physical mocks | TRI-20 to TRI-21 | Mock-community design surrogates and laboratory validation targets |
| Environmental | TRI-30 | Matched environmental application design surrogate |
| Showcase | TRI-31 | 30-minute home-computer demo design with expected report assets |

Dataset inventory:

| Dataset | Tier | Libraries | Description | Primary metric |
| --- | --- | ---: | --- | --- |
| TRI-00 | Software/CI | 0 | Metadata contract tests for valid manifests, controls, and deterministic validation errors | Schema validation |
| TRI-01 | Software/CI | 12 | Tiny end-to-end smoke test with biological, negative-control, and positive-control inputs | Workflow completion |
| TRI-02 | Software/CI | 18 | Assay-branch fixture that isolates WMS, 16S, and ITS processing failures | Branch correctness |
| TRI-03 | Synthetic truth | 18 | Balanced compact mixed-kingdom truth community for the CI/resource bundle | Cross-assay status accuracy |
| TRI-05 | Synthetic truth | 24 | Compact low-biomass and batch-specific contamination gradient | Contaminant precision and recall |
| TRI-06 | Synthetic truth | 18 | Compact taxonomy challenge with synonyms, rank differences, historical names, and held-out references | Mapping and comparison-rank accuracy |
| TRI-07 | Synthetic truth | 25 | Compact missing-modality, low-depth, host-dominated, empty, and failed-library challenge | Failure interpretation accuracy |
| TRI-10 | Synthetic truth | 21 | Even bacterial, archaeal, and fungal mock community | Taxonomic F1 |
| TRI-11 | Synthetic truth | 24 | Log-abundance mock spanning six orders of magnitude | Detection calibration |
| TRI-12 | Synthetic truth | 112 | Two-batch, six-level low-biomass contamination ladder for paired-control filtering | Contaminant precision and recall |
| TRI-13 | Synthetic truth | 75 | Structured environmental gradient across four habitat states | Ecological effect recovery |
| TRI-14 | Synthetic truth | 48 | Taxonomic ambiguity challenge with synonyms, rank limits, close relatives, and held-out references | Agreement-rank accuracy |
| TRI-15 | Synthetic truth | 97 | Batch and control-linkage challenge comparing paired controls against global-control assumptions | Batch-aware contamination accuracy |
| TRI-16 | Synthetic truth | 30 | Missing modalities and explicit library failures | Missingness classification |
| TRI-17 | Synthetic truth | 36 | Depth, primer-recovery, copy-number, and held-out-reference factorial | Detection and explanation calibration |
| TRI-18 | Synthetic truth | 54 | Sequencing-run carryover and index-leakage challenge with stage-linked controls | Carryover precision and recall |
| TRI-19 | Synthetic truth | 84 | Technical and extraction repeatability across replicates, runs, and batches | Evidence-status repeatability |
| TRI-20 | Physical mock surrogate | 45 | Commercial physical-standard acquisition and analysis design surrogate | Empirical taxonomic accuracy target |
| TRI-21 | Physical mock surrogate | 114 | Environmental mixed-kingdom mock-community acquisition and analysis design surrogate | Kingdom-specific calibration target |
| TRI-30 | Environmental surrogate | 192 | Balanced 48-sample matched built-environment study design surrogate | Robust ecological conclusions |
| TRI-31 | Showcase | 51 | Thirty-minute showcase with four environmental states, controls, host contamination, reagent/batch contamination, assay bias, taxonomy ambiguity, and expected report assets | Cross-assay explanatory narrative |

The inventory intentionally skips unused IDs such as TRI-04, TRI-08, and
TRI-09; the authoritative machine-readable list is
[benchmarks/registry.tsv](benchmarks/registry.tsv).

Run the benchmark metadata and fixture checks:

```bash
scripts/run_benchmark_suite.sh
```

Exercise generated datasets through the deterministic Nextflow graph:

```bash
scripts/run_benchmark_suite.sh --nextflow
```

Generate one dataset:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-12
```

Generate the showcase dataset:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-31
```

See [benchmarks/README.md](benchmarks/README.md) for dataset designs, truth
scales, publication package targets, and reproducibility requirements.

## Real-World Validation Downloads

Use `scripts/download_real_world_datasets.py` to build ENA-backed manifests
for empirical validation datasets under `data/real_world/`. The script defaults
to metadata and file manifests only; add `--download` when you are ready to
fetch FASTQs.

```bash
python3 scripts/download_real_world_datasets.py --dataset TRI-REAL-01 --max-groups 30
python3 scripts/download_real_world_datasets.py --dataset TRI-REAL-01 --max-groups 30 --download --verify-md5
```

Configured targets include the primary Shaffer mixed bacterial/fungal mock
benchmark (`TRI-REAL-01`, ENA `ERP124610`) plus optional environmental and
stress-test cohorts. Ambiguous marker-gene records are labeled `AMPLICON`
in `selected_runs.tsv` and mapped to `16S` by default in draft `samples.tsv`;
curate source study metadata before making biological claims. See
[docs/real-world-validation.md](docs/real-world-validation.md).

## Repository Layout

```text
.
|-- main.nf                         # Nextflow workflow entrypoint
|-- nextflow.config                 # Default params and execution profiles
|-- run_triarchy.sh                 # Lite/stub/production workflow wrapper
|-- install.sh                      # Local runtime, environment, and reference installer
|-- environment.yml                 # Conda/Mamba production environment
|-- bin/
|   `-- triarchy                    # Preflight and database helper CLI
|-- conf/
|   |-- assay_performance.tsv       # Assay reliability and power defaults
|   |-- evidence_policy.json        # Evidence thresholds and calibration metadata
|   `-- databases.example.yml       # Reference manifest example
|-- modules/local/
|   |-- metadata.nf                 # Metadata validation and manifest splitting
|   |-- wms.nf                      # WMS processing processes
|   |-- amplicon.nf                 # 16S and ITS processes
|   `-- integration.nf              # Evidence, diversity, report, provenance
|-- workflow/scripts/
|   |-- triarchy_integrate.py       # Integration commands used by Nextflow
|   |-- triarchy_evidence.py        # Evidence model implementation
|   `-- dada2_profile.R             # DADA2 profile generation
|-- scripts/
|   |-- validate_benchmark_metadata.py
|   |-- generate_benchmark_datasets.py
|   |-- download_real_world_datasets.py
|   |-- install_databases.sh
|   `-- run_benchmark_suite.sh
|-- benchmarks/                     # Deterministic fixtures and truth tables
|-- data/                           # Local real-world validation downloads
|-- docs/                           # Methodology, running, contracts, results
|-- tests/                          # Unit and integration tests
|-- references/                     # Local reference destination/checkpoints
|-- results/                        # Workflow outputs
`-- work/                           # Nextflow work directory
```

## Command Line Helpers

The primary execution path is `run_triarchy.sh`, but `bin/triarchy` provides
cheap preflight checks before scheduling a large production run.

Validate sample metadata and control edges:

```bash
bin/triarchy validate \
  --input benchmarks/TRI-01/metadata/samples.tsv \
  --paired-controls benchmarks/TRI-01/metadata/control_relationships.tsv
```

Validate FASTQ paths and references before production:

```bash
bin/triarchy validate \
  --input /data/study/samples.tsv \
  --paired-controls /data/study/control_relationships.tsv \
  --databases conf/databases.example.yml \
  --database-destination /refs/triarchy \
  --report results/preflight.json
```

Check a database manifest independently:

```bash
bin/triarchy databases check \
  --manifest conf/databases.example.yml \
  --destination references/data
```

Run the stable database installer helper:

```bash
bin/triarchy databases install \
  --profile minimal \
  --destination references/data \
  --manifest conf/databases.example.yml \
  --yes
```

Helper command reference:

| Command | Option | Meaning |
| --- | --- | --- |
| `validate` | `--input`, `--samples` | Sample manifest TSV |
| `validate` | `--paired-controls`, `--controls` | Sample-control edge TSV |
| `validate` | `--databases` | Optional reference manifest to check |
| `validate` | `--database-destination` | Root directory for relative manifest paths |
| `validate` | `--skip-fastq-check` | Validate metadata without requiring FASTQ paths to exist |
| `validate` | `--report` | Write JSON validation report |
| `databases check` | `--manifest` | Reference manifest to check |
| `databases check` | `--destination` | Root directory for relative reference paths |
| `databases check` | `--report` | Write JSON database-check report |
| `databases install` | `--destination` | Reference install destination |
| `databases install` | `--profile` | `minimal` or `standard` |
| `databases install` | `--manifest` | Optional manifest to verify after installation |
| `databases install` | `--kraken-archive` | Existing Kraken2/Bracken archive |
| `databases install` | `--download-attempts` | Override retry attempts |
| `databases install` | `--yes` | Non-interactive install |
| `databases install` | `--report` | Write JSON install/check report |

## Configuration And References

Key configuration files:

| File | Purpose |
| --- | --- |
| `.nextflow-version` | Pinned Nextflow version |
| `nextflow.config` | Params, resource labels, profiles, reports, trace, DAG |
| `conf/assay_performance.tsv` | Assay reliability, depth, coverage, and power defaults |
| `conf/evidence_policy.json` | Evidence thresholds, priors, model version, calibration scope |
| `conf/databases.example.yml` | Reference manifest schema and example entries |
| `conf/references.config` | Generated production reference paths after install |

Execution profiles include:

- `local`: local executor defaults;
- `lite`: one-task-at-a-time deterministic laptop profile;
- `demo`: TRI-31 showcase defaults and moderate local resource caps;
- `laptop`: constrained local production profile;
- `workstation`: larger local production profile;
- `high_memory`: higher-memory local process labels;
- `docker`, `apptainer`, `singularity`: container-enabled profiles;
- `slurm`: cluster execution profile hook;
- `test`: test fixture defaults.

Production reference installation writes checkpoint receipts under the
configured reference destination and records the last successful destination in
`references/.database-dir`.

## Database Manifest Format

Reference manifests let TRIARCHY check that production references exist and,
for files, optionally verify checksums. The example manifest is
`conf/databases.example.yml`.

Minimal YAML shape:

```yaml
schema_version: 1
profile: minimal
description: Example TRIARCHY reference manifest.
references:
  - id: host_bowtie2
    kind: bowtie2_index
    version: artificial-minimal
    path: host/artificial_host
    required: true
  - id: kraken2_bracken
    kind: kraken2_database
    version: minimal-test
    path: kraken2
    required: true
  - id: silva_16s
    kind: dada2_taxonomy
    version: example
    path: marker/16s_taxonomy.fa.gz
    required: true
    sha256: optional_checksum_for_single_files
```

Reference fields:

| Field | Required | Meaning |
| --- | --- | --- |
| `id` | Yes | Stable reference identifier |
| `kind` | Yes | Reference type, such as `bowtie2_index`, `kraken2_database`, or `dada2_taxonomy` |
| `version` | Yes | Database or reference release label |
| `path` | Yes | Path to the directory, prefix, or file |
| `required` | Yes | Whether a missing reference fails validation |
| `sha256` | For checksummed files | Expected SHA-256 for a single file |

Relative `path` values are resolved beneath the `--destination` passed to
`bin/triarchy databases check`, or beneath the manifest directory when no
destination is supplied.

Check references:

```bash
bin/triarchy databases check \
  --manifest conf/databases.example.yml \
  --destination references/data \
  --report results/database-check.json
```

The command returns nonzero when a required reference is missing, empty, or has
a checksum mismatch.

## Reference Installation Details

The production database installer is called by `install.sh`, and is also
available through:

```bash
bin/triarchy databases install \
  --profile standard \
  --destination /refs/triarchy \
  --manifest conf/databases.example.yml \
  --cpus 8 \
  --yes
```

Reference stages are checkpointed. A checkpoint is trusted only when its
corresponding output also passes validation. Existing valid output without a
receipt can be adopted and checkpointed automatically.

Standard-profile stages include:

- host reference retrieval and Bowtie2 index construction;
- SILVA-compatible 16S taxonomy reference installation;
- ITS fungal taxonomy reference installation;
- prebuilt Kraken2/Bracken PlusPFP archive download, checksum validation, and
  extraction;
- HUMAnN nucleotide database installation;
- HUMAnN protein database installation;
- `conf/references.config` generation.

The standard Kraken2/Bracken archive is large. If you download it separately,
wait for the transfer to finish before passing it with `--kraken-archive`; the
installer checks the exact byte size before calculating the checksum.

Production releases should replace example manifest paths and checksums with
the exact submitted reference locations, release labels, and file hashes.

## Development And Testing

Run unit tests:

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
```

Validate the TRI-00 metadata fixture:

```bash
python3 scripts/validate_benchmark_metadata.py \
  benchmarks/TRI-00/valid/samples.tsv \
  benchmarks/TRI-00/valid/control_relationships.tsv
```

Run a deterministic workflow smoke test:

```bash
./run_triarchy.sh --lite
```

Run benchmarks:

```bash
scripts/run_benchmark_suite.sh
scripts/run_benchmark_suite.sh --nextflow
```

When editing shell scripts on Windows, preserve LF line endings. The repository
includes `.gitattributes` for this because WSL cannot execute scripts whose
shebang has CRLF endings.

## Troubleshooting

### `/usr/bin/env: 'bash\r': No such file or directory`

The script has Windows CRLF line endings. Convert scripts to LF:

```bash
sudo apt install -y dos2unix
dos2unix install.sh run_triarchy.sh scripts/*.sh bin/triarchy
```

The repository includes `.gitattributes` to keep shell scripts LF-normalized on
future checkouts.

### `Nextflow is not installed; run ./install.sh first`

Install the lightweight runtime:

```bash
./install.sh --skip-dependencies --skip-databases --skip-tests --yes
```

Then rerun:

```bash
./run_triarchy.sh --lite
```

### Java Version Errors

Nextflow requires Java 17 or later. Check:

```bash
java -version
```

On WSL Ubuntu:

```bash
sudo apt update
sudo apt install -y openjdk-17-jre
```

### Metadata Validation Fails

Run validation directly to get stable error messages:

```bash
python3 scripts/validate_benchmark_metadata.py \
  path/to/samples.tsv \
  path/to/control_relationships.tsv
```

Common causes:

- duplicate `library_id` values;
- typo in `assay`, `sample_role`, `library_layout`, or `control_stage`;
- paired library missing `fastq_r2`;
- control edge references a missing library;
- control edge links across assays;
- control target is not a control sample;
- matched assay rows disagree on biological sample metadata.

### Production Tool Missing

If a production task fails with a message such as `fastp is required`,
`kraken2 is required`, or `Rscript is required`, the production environment is
not active or was not installed.

Use the wrapper after a full install so `.tools/env/bin` is added to `PATH`:

```bash
./install.sh --database-profile standard --yes
./run_triarchy.sh -profile local ...
```

For an existing environment:

```bash
export PATH="$PWD/.tools/env/bin:$PWD/.tools/bin:$PATH"
export NXF_HOME="$PWD/.nextflow"
export NXF_VER="$(<.nextflow-version)"
```

### Conda Solver Cannot Find `xopen`

If `mamba env create` fails while installing `cutadapt` with a message such as
`xopen >=1.6.0 ... does not exist`, refresh the channel metadata and rerun the
installer from the repository root:

```bash
mamba clean --index-cache
./install.sh --database-profile standard --database-dir /path/to/reference/storage --cpus 8 --yes
```

If the solver still fails and `environment.lock.txt` is present for the same
Linux platform, bypass solving and use the pinned explicit environment:

```bash
./install.sh --locked-dependencies --database-profile standard --database-dir /path/to/reference/storage --cpus 8 --yes
```

### Reference Or Database Missing

Run the database checker before scheduling production:

```bash
bin/triarchy databases check \
  --manifest conf/databases.example.yml \
  --destination /refs/triarchy
```

If `conf/references.config` is missing, rerun the installer without
`--skip-databases` or pass all reference parameters explicitly to
`run_triarchy.sh`.

### A Run Fails Midway

After fixing the input, tool, or reference problem, resume with the same
parameters:

```bash
./run_triarchy.sh -resume [same parameters]
```

Inspect failed tasks:

```bash
nextflow log
nextflow log last -f process,hash,status,exit
```

Nextflow prints the failed work directory. Inspect:

```text
work/<hash>/.command.sh
work/<hash>/.command.err
work/<hash>/.command.log
```

### WSL Performance Is Poor Under `/mnt/c`

For real workflows, copy the repository and data into the Ubuntu filesystem:

```bash
cp -a /mnt/c/Users/agluc/Documents/MASON-LAB/TRIARCHY ~/TRIARCHY
cd ~/TRIARCHY
```

Large Nextflow work directories are usually faster under `~/...` than under
Windows-mounted paths.

### Conda Or Mamba Is Slow

Prefer Mamba or Micromamba for the production environment:

```bash
./install.sh --environment-manager mamba --database-profile minimal --yes
```

Use `--skip-dependencies` only for lite/stub development modes that do not need
production assay tools.

### Apple-Silicon Notes

On Apple-silicon macOS, the installer configures native Homebrew Bowtie2,
Kraken2, and Bracken shims to avoid incompatible Bioconda executables. Do not
use `--skip-dependencies` for a production Apple-silicon setup unless you are
managing the native tool bridge yourself.

## Glossary

| Term | Meaning |
| --- | --- |
| WMS | Whole-metagenome shotgun sequencing |
| 16S | 16S rRNA gene marker sequencing for bacteria and archaea |
| ITS | Internal transcribed spacer marker sequencing for fungi |
| Assay-native | Kept on the measurement scale of the assay that produced it |
| Biological sample | A real specimen or environmental sample being studied |
| Library | One sequenced assay-specific preparation |
| `biological_sample_id` | Identifier used to connect matched WMS, 16S, and ITS libraries |
| Negative control | Blank or control expected not to contain target biological material |
| Positive control | Mock or control expected to contain known material |
| Sample-control graph | Explicit links from biological libraries to relevant controls |
| Control stage | Field, extraction, PCR, library, or sequencing stage where a control is relevant |
| Linked-control evidence | Contamination evidence estimated from controls explicitly connected to a sample |
| Taxonomic harmonization | Mapping source classifier names/IDs to canonical comparable taxon records |
| Detection power | Estimated ability of an assay to detect a taxon if it is present |
| Reference coverage | Whether the relevant database can support detection or nondetection for a taxon |
| Rank-limited agreement | Agreement only at a higher shared taxonomic rank |
| Biological state | TRIARCHY's `present`, `absent`, or `indeterminate` output state |
| Stub run | Nextflow run that executes process stubs rather than production tools |
| Lite mode | A laptop-friendly wrapper mode using stub execution and conservative resources |
| Production mode | Full assay processing with real FASTQs, tools, and references |

## Project Status

This repository contains the executable workflow, local installer, metadata
validator, integration scripts, report builder, deterministic benchmark
fixtures, and tests needed for development and orchestration checks.

The implemented claim is the control-aware evidence integration method, not
the existence of generic WMS, 16S, or ITS processing branches. Biological
validation still depends on running the staged benchmark designs described in
[docs/methodology.md](docs/methodology.md) and [benchmarks/README.md](benchmarks/README.md).

TRIARCHY is being prepared as a Bioinformatics Application Note. A public
publication release should include:

- source code at https://github.com/aglucaci/TRIARCHY;
- an archived source release with DOI;
- archived benchmark data and truth tables;
- container tags and immutable digests;
- a checked database manifest;
- an example `triarchy_report.html`;
- release notes and a maintenance commitment.

## License

The repository is currently licensed as GPL-3.0-only. See [LICENSE](LICENSE).
