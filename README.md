# TRIARCHY

![TRIARCHY logo](logo/triarchy_logo.png)

**Control-aware integration of shotgun metagenomic, 16S rRNA gene, and fungal
ITS sequencing for environmental microbiome studies.**

TRIARCHY is a multi-assay microbiome framework for matched environmental
samples. It treats whole-metagenome shotgun sequencing (WMS), 16S, and ITS as
complementary, biased observations of a microbial community—not as
interchangeable abundance measurements.

## Quick start

Requirements:

- a POSIX environment (Linux, macOS, or WSL);
- Bash;
- Java 17 or later;
- Python 3 for metadata validation and repository tests;
- Mamba (default), Micromamba, or Conda;
- `curl` for the local Nextflow installation.

From the TRIARCHY repository, run the installer and then launch the included
deterministic stub workflow:

```bash
cd /path/to/TRIARCHY
./install.sh --database-profile minimal --yes
./run_triarchy.sh --stub
```

The minimal database profile is intended only for installation checks and stub
runs; it is not suitable for biological analysis. Stub results are written to
`results/stub`.

To resume an interrupted stub run:

```bash
./run_triarchy.sh --stub -resume
```

The wrapper automatically configures the repository's local Nextflow and Conda
environment, so no manual `PATH`, `NXF_HOME`, or `NXF_VER` exports are needed
when using `run_triarchy.sh`.

## Installation

For a complete production installation, including the standard reference
databases, run:

```bash
./install.sh --yes
```

The complete installer:

- installs the pinned Nextflow runtime under `.tools/bin`;
- creates `.tools/env` with all WMS, 16S, ITS, R, and Python dependencies;
- downloads the prebuilt Kraken2 PlusPFP/Bracken database and installs the
  host, HUMAnN, SILVA, and ITS references under `references/data`;
- records installed packages and database provenance;
- generates `conf/references.config`;
- validates TRI-00 and runs the test suite.

The standard reference installation can exceed 400 GB and take many hours, so
`--yes` is required for non-interactive installation. For a smaller
installation check, use:

```bash
./install.sh --database-profile minimal --yes
```

Database installation is checkpointed by validated stage. Rerunning the same
command skips completed references and resumes the prebuilt Kraken archive
transfer. Its official MD5 manifest validates both the archive and every
extracted database file. The terminal shows colored `[step/total percent]`
progress; set `NO_COLOR=1` for plain output.

When `--kraken-archive` points to a separately downloaded archive, allow that
transfer to finish before starting the installer. The installer checks its
exact byte size before calculating the checksum.

Minimal databases are not suitable for biological analysis. Run
`./install.sh --help` for component-level skip and destination options.

## Running the workflow

After installation, verify the complete graph with the included metadata
fixture and deterministic stub outputs:

```bash
./run_triarchy.sh --stub
```

The TRI-00 FASTQ paths are intentionally synthetic, so this fixture must use
`-stub-run`. A production run uses the same entrypoint without `-stub-run` and
requires real FASTQs, assay tools, primers, taxonomy references, and WMS
databases.

The wrapper configures `PATH`, `NXF_HOME`, and the pinned `NXF_VER`, then passes
all remaining arguments to Nextflow. Use `./run_triarchy.sh --help` for
production and stub examples.

See [running.md](docs/running.md) for production parameters, software
requirements, execution profiles, outputs, and restart behavior.

Run the complete multi-dataset benchmark suite with:

```bash
scripts/run_benchmark_suite.sh
```

Add `--nextflow` to execute every generated dataset through the full
deterministic Nextflow graph. See
[benchmarks/README.md](benchmarks/README.md) for their truth designs.

> TRIARCHY does not force three assays into one abundance table. It determines
> which biological conclusions survive their different biases, resolutions,
> and contamination patterns.

## Contribution boundary

TRIARCHY does not claim to be the first workflow to process WMS, 16S, and ITS,
the first cross-assay harmonization method, or the first microbiome
contamination filter. Its testable contribution is narrower: a stage-linked,
control-aware evidence layer that explains sample-level taxon agreement and
disagreement using measured assay availability, detection power, taxonomic
rank, reference coverage, and contamination evidence. Existing classifiers,
Greengenes2-style bacterial reconciliation, decontam, Com-2seq, and other
multi-table methods are comparators or upstream evidence sources rather than
claims displaced by TRIARCHY.

## Scientific premise

- **16S** preferentially observes bacterial and archaeal marker genes.
- **ITS** preferentially observes fungal communities.
- **WMS** observes taxonomic and functional potential across kingdoms, with
  sensitivity determined by depth, biomass, extraction, and reference coverage.
- Relative abundances from these assays are on different scales and are never
  averaged directly.
- Negative and positive controls are evidence. They remain linked to the
  biological samples, laboratory stages, and batches for which they are
  relevant.

## How TRIARCHY works

1. **Validate the experimental design.** Input manifests become a sample–control
   graph linking every library to its material, environment, assay, batches,
   and stage-specific controls.
2. **Process each assay independently.**
   - **WMS:** fastp quality control, Bowtie2 host removal, Kraken2/Bracken
     taxonomy, and HUMAnN functional profiling.
   - **16S:** cutadapt trimming, DADA2 ASV inference, and bacterial/archaeal
     taxonomy.
   - **ITS:** cutadapt trimming, ITSxpress extraction, DADA2 ASV inference, and
     fungal taxonomy.
3. **Annotate contamination.** Each feature is evaluated using linked-control
   prevalence, abundance enrichment, and batch localization. TRIARCHY preserves
   unfiltered, control-annotated, and high-confidence profiles.
4. **Harmonize taxonomy with provenance.** Source assignments are mapped to
   canonical taxa without discarding original names, identifiers, database
   versions, or mapping confidence. Assays are compared only at the lowest rank
   valid in both.
5. **Classify cross-assay evidence.** Bacterial and archaeal calls compare 16S
   with WMS; fungal calls compare ITS with WMS. Functional WMS evidence remains
   separate from taxonomic evidence. Taxa receive detection-power estimates,
   explanatory classes, a three-state biological conclusion, an authenticity
   probability with explicit calibration scope, and evidence ranks.
6. **Analyze ecology without collapsing assays.** Alpha and beta diversity are
   calculated within each assay. Cross-assay analysis compares ecological
   conclusions, distance structures, and effect sizes.
7. **Report confidence and disagreement.** The final report explains which
   samples and taxa are trustworthy, where assays disagree, why calls were
   made, and which biological conclusions remain robust.

See [methodology.md](docs/methodology.md) for the evidence model and
[data-contracts.md](docs/data-contracts.md) for the input and output schemas.
The [benchmark suite](benchmarks/README.md) defines truth-based software,
synthetic, physical-mock, and environmental validation datasets.

## Evidence classes

| Status | Interpretation |
| --- | --- |
| Concordant detection | All kingdom-appropriate assays support the taxon |
| Single-assay credible | One appropriate assay supports the taxon and no measured assay contradicts it |
| Assay-limited nondetection | The other assay's estimated detection power is insufficient |
| Rank-limited agreement | Assays agree only at a shared higher taxonomic rank |
| Reference-limited | A database lacks adequate representation for a defensible negative call |
| Probable contaminant | The taxon is enriched in relevant negative controls |
| Cross-contamination candidate | Control-associated evidence is localized to a batch or run |
| Genuine discordance | Technically adequate assays provide incompatible evidence |
| Indeterminate | Depth, controls, or taxonomic resolution are insufficient |

These labels are evidence summaries, not declarations of biological truth.

## Data products

TRIARCHY keeps three distinct levels of output:

```text
profiles/
├── unfiltered/          # assay-native counts and normalized values
├── control_annotated/   # all features with contamination evidence
└── high_confidence/     # decisions passing the configured evidence policy

integrated/
├── sample_confidence.tsv
├── taxon_evidence.tsv
├── evidence_model.json
├── taxonomy_crosswalk.tsv
├── cross_assay_concordance.tsv
└── provenance.json
```

Assay-native counts are retained throughout. Cross-assay products contain
detection evidence, confidence, contamination status, concordance, and
standardized effect summaries—not averaged relative abundance.

`taxon_evidence.tsv` distinguishes measured nondetection from missing-by-design
and reports assay depth, estimated detection power, reference coverage, linked
control evidence, an explanatory class, a three-state biological conclusion,
and `probability_authentic`. The default probability mapping is labeled
`synthetic_only_not_empirically_calibrated`; it is not empirical clinical or
environmental calibration.

## Modality summary

| Step | WMS | 16S | ITS |
| --- | --- | --- | --- |
| Read QC | fastp | cutadapt | cutadapt |
| Host removal | Bowtie2 | — | — |
| ASV inference | — | DADA2 | DADA2 |
| ITS extraction | — | — | ITSxpress |
| Taxonomy | Kraken2/Bracken | SILVA | EUKARYOME |
| Functional profiling | HUMAnN | — | — |
| Linked-control assessment | Yes | Yes | Yes |
| Within-assay diversity | Yes | Yes | Yes |
| Cross-assay comparison | WMS ↔ 16S bacteria/archaea | 16S ↔ WMS | ITS ↔ WMS fungi |

## Project status

This repository contains the executable Nextflow workflow, a pinned local
software environment, automated reference installation, metadata validation,
unit and integration tests, and the TRI-00 orchestration fixture. Biological
validation still depends on running the staged read-based benchmark datasets
described in [methodology.md](docs/methodology.md) and
[benchmarks/README.md](benchmarks/README.md).

## Workflow

See [flowchart.md](flowchart.md) for the complete control-aware data flow.
