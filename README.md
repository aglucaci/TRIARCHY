# TRIARCHY

![TRIARCHY logo](assets/logo/triarchy_logo_transparent.png)

**TRIARCHY is a production workflow for control-aware integration of matched
whole-genome shotgun (WGS) metagenomics, 16S rRNA gene sequencing, and fungal
ITS sequencing in environmental microbiome studies.**

TRIARCHY processes each assay independently, links every production library to
its experimentally relevant controls, harmonizes taxonomy with provenance, and
classifies where assays agree, disagree, or lack enough power for a defensible
conclusion. This repository supports production execution only.

## Requirements

- Linux, macOS, or WSL Ubuntu
- Bash and Java 17 or later
- Mamba, Micromamba, or Conda
- Nextflow (the pinned version is installed locally)
- Production reference storage; the complete reference set can exceed 400 GB
- Real WGS, 16S, and ITS FASTQs with study-specific primer information

Apple-silicon macOS production installs use native tool shims where required.
For large studies, use a workstation, HPC scheduler, or cloud executor rather
than a memory-constrained laptop.

## Install

Create the named environment and install the complete production references:

```bash
mamba env create --file environment.yml
mamba activate TRIARCHY

./install.sh \
  --database-dir /path/to/production/references \
  --cpus 8 \
  --yes
```

The installer downloads resumably, validates completed stages, writes
`conf/references.config`, and records the selected reference root in
`references/.database-dir`.

Useful installer options:

| Option | Meaning |
| --- | --- |
| `--database-dir DIR` | Production reference destination |
| `--silva-archive FILE` | Install a pre-existing, checksum-verified SILVA archive |
| `--kraken-archive FILE` | Reuse an existing full PlusPFP archive |
| `--cpus N` | Installation/indexing threads |
| `--locked-dependencies` | Install the pinned dependency lock |
| `--skip-nextflow` | Keep the existing Nextflow installation |
| `--skip-dependencies` | Keep the existing `TRIARCHY` environment |
| `--skip-databases` | Keep the current production references |
| `--yes` | Accept the large non-interactive installation |

SILVA downloads are automatic and mirror-aware. TRIARCHY tries the official
SILVA release server first and then a checksum-identical Zenodo mirror. Each
source has its own resumable partial download, unreachable sources fail over
automatically, and the completed archive must match MD5
`41f1d231ce5079b88f4c01742c288ba9`. Existing installations of the earlier
full-length DADA2 artifact remain accepted by their pinned checksum. The
optional `--silva-archive` input is only needed for fully air-gapped systems.

## Production inputs

TRIARCHY requires two tab-separated metadata files.

### Sample manifest

Required columns:

```text
library_id
biological_sample_id
specimen_id
subject_id
assay
sample_role
sample_material
environment_category
library_layout
fastq_r1
fastq_r2
extraction_batch
pcr_plate
library_batch
sequencing_run
technical_replicate
```

Valid assays are `WMS`, `16S`, and `ITS`. Valid sample roles are `biological`,
`negative_control`, `positive_control`, and `mock`. Paired libraries require
both FASTQ paths. Relative paths are resolved from the manifest directory.

### Control relationships

Required columns:

```text
library_id
control_library_id
control_stage
relationship
priority
reason
```

Control stages are `field`, `extraction`, `pcr`, `library`, and `sequencing`.
Controls must use the same assay as the production library to which they are
linked.

Validate metadata before launching an expensive run:

```bash
triarchy validate \
  --input /data/study/samples.tsv \
  --paired-controls /data/study/control_relationships.tsv
```

## Run production

Run the production workflow with real FASTQs:

```bash
triarchy run \
  --profile local \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --outdir results/study
```

Assay primers and the ITS target region are maintained centrally in
`conf/primers.config`; users do not pass them on the command line. Update that
controlled configuration when a production protocol changes.

Resume an interrupted run with `--resume`:

```bash
triarchy run \
  --profile local \
  --resume \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --outdir results/study
```

Available production resource profiles:

| Profile | Intended system |
| --- | --- |
| `laptop` | Constrained local production; queue size 2, high tasks capped at 8 GB |
| `local` | Default local executor resources |
| `workstation` | Local workstation with at least 32 GB RAM |
| `high_memory` | High-memory host; high tasks request 64 GB |
| `slurm` | SLURM cluster execution |
| `cloud` | Cloud executor configuration |
| `docker` | Production Docker execution |
| `apptainer` | Production Apptainer execution |

## Production outputs

The output directory contains:

```text
qc/library_qc.tsv
profiles/unfiltered/taxa.tsv
profiles/unfiltered/functions.tsv
profiles/control_annotated/taxa.tsv
profiles/control_annotated/functions.tsv
profiles/high_confidence/taxa.tsv
profiles/sample_confidence.tsv
integrated/harmonized_taxa.tsv
integrated/taxonomy_crosswalk.tsv
integrated/taxon_evidence.tsv
integrated/cross_assay_concordance.tsv
integrated/evidence_model.json
diversity/alpha_diversity.tsv
diversity/beta_diversity.tsv
pipeline_info/provenance.json
pipeline_info/execution_report.html
pipeline_info/timeline.html
triarchy_report.html
```

Open `triarchy_report.html` for the scientific evidence dashboard. It follows
the interpretation sequence:

1. Did the experiment work?
2. Which libraries and samples are trustworthy?
3. Which taxa are supported?
4. Why do assays disagree?
5. Which ecological conclusions survive?
6. Can the evidence and provenance be exported?

Cross-assay classifications summarize evidence. They do not make WGS, 16S,
and ITS relative abundances quantitatively interchangeable.

## Workflow architecture

The production graph is:

```text
metadata validation
  → assay-specific QC
  → WGS host removal, taxonomy, and functions
  → 16S trimming and DADA2
  → ITS trimming, ITSxpress, and DADA2
  → linked-control annotation
  → taxonomy harmonization
  → authenticity and discordance classification
  → assay-native diversity and cross-assay concordance
  → dashboard and provenance
```

The implementation lives in `main.nf`, `modules/local/`, and
`workflow/scripts/`.

## Troubleshooting

Inspect a failed Nextflow task under its `work/<hash>/` directory:

```text
.command.sh
.command.err
.command.log
```

Confirm the named environment and production references before rerunning:

```bash
mamba activate TRIARCHY
nextflow config -flat | grep -E 'host_index|kraken2_db|humann|taxonomy'
```

On laptops, prefer `-profile laptop`, close memory-intensive applications, and
place the Nextflow work directory and full references on fast external storage.

## License

See [LICENSE](LICENSE).
