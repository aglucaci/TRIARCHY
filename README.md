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
- Production reference storage; the complete reference set can exceed 500 GB
- Real WGS, 16S, and ITS FASTQs with study-specific primer information

Apple-silicon macOS production installs use native tool shims where required.
For large studies, use a workstation, HPC scheduler, or cloud executor rather
than a memory-constrained laptop.

## Minimum production specifications

For a complete production run using the full Kraken2 PlusPFP, MetaPhlAn,
HUMAnN, host, SILVA, and ITS reference set, provision at least:

| Resource | Minimum production specification |
| --- | --- |
| Largest compute-job slot | 12 CPU cores |
| Memory per high-resource job | 96 GB RAM |
| Maximum job wall time | 48 hours |
| Free space for reference installation | 600 GB |
| Recommended shared storage | 1 TB or more |
| GPU | Not required |

The built-in SLURM profile requests the following resources per task:

| Task class | CPUs | Memory | Wall time |
| --- | ---: | ---: | ---: |
| Low | 1 | 2 GB | 1 hour |
| Medium | 4 | 16 GB | 8 hours |
| High | 12 | 96 GB | 48 hours |

For HPC execution, use Linux (preferably x86-64), SLURM with `sbatch`, Bash,
Java 17 or later, and the named `TRIARCHY` Mamba/Conda environment. The
repository, reference databases, input FASTQs, result directory, and Nextflow
work directory must be on shared POSIX storage visible at identical paths from
the login and compute nodes. Internet access is needed during database
installation unless the validated archives are supplied locally.

The database installer enforces 600 GB of free space for a new reference
installation. Production work files and results require additional capacity;
1 TB of shared storage is a practical starting point, with more space needed
for large input studies. Nodes with 128 GB RAM are recommended to leave
operating-system and database-loading headroom above the 96 GB process request.
The `laptop` profile's 8 GB high-task cap is intended for the constrained demo
and is not a dependable baseline for full production datasets.

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

TRIARCHY uses Kraken2/Bracken for WMS taxonomy and installs the complete HUMAnN
ChocoPhlAn and UniRef90 databases for functional profiling. It also installs
the checksum-pinned `mpa_vJun23_CHOCOPhlAnSGB_202307` MetaPhlAn marker database
required by HUMAnN 3.9. HUMAnN prescreens against this shared index before
constructing a small sample-specific ChocoPhlAn search database. The database
path and release are passed explicitly with offline lookup enabled, so compute
nodes do not need internet access and do not build the full ChocoPhlAn index.

The full Kraken2 PlusPFP database is larger than the 96 GB production job
request. TRIARCHY therefore enables Kraken2's supported `--memory-mapping` mode
and limits WMS taxonomy to two concurrent jobs on SLURM. Place the database on
fast parallel or node-local storage when possible; memory-mapped classification
trades lower process memory for additional storage I/O.

## Production demo

TRIARCHY includes a deterministic production demonstration built from
version-pinned public NCBI nucleotide sequences. It contains 21 WMS, 16S, and ITS
libraries spanning three matched biological samples, stage-specific negative
controls, a mixed positive control, human mitochondrial reads, and designed
cross-assay edge cases. It uses the normal production workflow and installed
reference databases—there are no stub tools or precomputed taxonomic calls.

The demo is a synthetic benchmark, not a biological cohort. It generates reads
from the following public reference records and records their exact sequence
checksums in `data/production_demo/metadata/references.tsv`:

| Target | Assay role | Version-pinned NCBI records |
| --- | --- | --- |
| *Escherichia coli* | WMS and 16S | [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3), [NR_114042.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_114042.1) |
| *Pseudomonas aeruginosa* | WMS and 16S | [NC_002516.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002516.2), [NR_026078.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_026078.1) |
| *Staphylococcus epidermidis* | WMS and 16S | [NC_004461.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_004461.1), [NR_113957.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_113957.1) |
| *Saccharomyces cerevisiae* | WMS and ITS | [NC_001133.9](https://www.ncbi.nlm.nih.gov/nuccore/NC_001133.9), [NC_001144.5](https://www.ncbi.nlm.nih.gov/nuccore/NC_001144.5) |
| *Candida albicans* | ITS | [FJ662402.1](https://www.ncbi.nlm.nih.gov/nuccore/FJ662402.1) |
| *Homo sapiens* mitochondrion | host-removal challenge | [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) |

The generated mixture contains 317,000 paired reads at the default scale. Its
truth table intentionally covers concordant bacterial and fungal detections,
assay-limited nondetection, low-biomass contamination, positive controls, and
approximately 15% mitochondrial host reads in one WMS library. Because the
benchmark is constructed from isolated reference sequences, it tests workflow
execution and evidence rules; it must not be used to estimate performance on
natural communities.

Build the dataset and run the complete production pipeline:

```bash
mamba activate TRIARCHY
triarchy demo --profile laptop --resume
```

The command downloads the small source references from NCBI, generates exactly
317,000 deterministic paired reads, validates the manifests and FASTQs, runs
TRIARCHY, and evaluates the resulting evidence and host-removal metrics against
the designed truth table. Outputs are written to `results/production_demo`.

To generate and validate the inputs without starting Nextflow:

```bash
triarchy demo --build-only
```

Use `--offline` after the NCBI references have been cached, or `--force` to
regenerate all reads with the selected `--seed`, `--scale`, and `--error-rate`.

On a SLURM HPC, keep the repository, production references, inputs, results, and
Nextflow work directory on shared storage visible to compute nodes:

```bash
mamba activate TRIARCHY
triarchy demo \
  --profile slurm \
  --site-config /shared/project/triarchy-site.config \
  --work-dir /shared/scratch/$USER/triarchy-work
```

The built-in SLURM profile requests 1 CPU/2 GB for low, 4 CPU/16 GB for medium,
and 12 CPU/96 GB for high-resource tasks. The optional site configuration is
where a cluster-specific partition, allocation, QoS, module, or scratch policy
should be declared. The activated `TRIARCHY` environment and database paths
must also be available at the same paths from submitted jobs.

## Public empirical cohort

`TRI-REAL-01` is a small empirical validation subset of the UC San Diego DNA
extraction comparison study [ENA ERP124610](https://www.ebi.ac.uk/ena/browser/view/ERP124610),
also registered as [BioProject PRJEB40903](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB40903)
and [Qiita study 12201](https://qiita.ucsd.edu/study/description/12201).
The source experiment compared extraction methods for microbial-community DNA
using 16S, ITS, and shotgun metagenomic sequencing [1].

The default TRIARCHY selection is deliberately small enough for an integration
run while retaining matched assay structure:

| Property | Default `TRI-REAL-01` selection |
| --- | --- |
| Biological material | Dust microbial-community samples from the extraction comparison |
| Selected groups | Three biological groups and one blank group |
| Libraries | 11 ENA runs represented by 15 FASTQ files |
| Assays | Three matched biological WMS/16S/ITS triplets; blank 16S and WMS libraries |
| Layout | Single-end amplicon reads and paired-end WMS reads |
| Download size | Approximately 2.96 MB |
| Integrity | ENA-provided MD5 verification plus gzip and manifest validation |
| Local provenance | `metadata/ena_runs.tsv`, `selected_runs.tsv`, `download_manifest.tsv`, and `dataset.json` |

This is a selected subset, not a redistribution or complete mirror of the
7,296-run source study. The downloader queries ENA at execution time, records
the returned accession metadata and URLs, and retrieves the chosen FASTQs
directly from ENA. `--max-groups` changes the subset and may materially change
its size and assay balance.

Download and run it with:

```bash
mamba activate TRIARCHY
triarchy real-world \
  --dataset TRI-REAL-01 \
  --profile slurm \
  --site-config /shared/project/triarchy-site.config \
  --work-dir /shared/scratch/$USER/triarchy-real-01-work \
  --outdir /shared/project/results/TRI-REAL-01
```

Run this command from a login or transfer node with internet access. FASTQs are
downloaded into `data/real_world/TRI-REAL-01`, verified against ENA MD5 values,
and then submitted through the normal production workflow. To stage the data
without launching Nextflow, add `--download-only`; later use `--skip-download`
to run the staged copy on an offline HPC.

The public metadata does not currently provide curated pairings between the
selected blank and biological libraries. The generated control relationship
manifest is therefore empty, and TRIARCHY reports this limitation rather than
inventing control links. Curate that manifest before interpreting contamination
classifications.

### Other catalogued public studies

The metadata downloader also catalogues larger studies for future curation.
They are not exposed through `triarchy real-world` because their ENA metadata
does not currently yield a validated, automatically paired WMS/16S/ITS subset:

| Dataset ID | Public study | Intended use | Current limitation |
| --- | --- | --- | --- |
| `TRI-REAL-02` | Arabidopsis leaf microbiomes, [ENA PRJEB31530](https://www.ebi.ac.uk/ena/browser/view/PRJEB31530) | Plant phyllosphere WMS with associated 16S and ITS data [2] | Assay runs require study-specific sample pairing rather than automatic ENA alias grouping. |
| `TRI-REAL-03` | Earth Microbiome Project Multi-omics (EMP500), [ENA PRJEB42019](https://www.ebi.ac.uk/ena/browser/view/PRJEB42019), [Qiita 13114](https://qiita.ucsd.edu/study/description/13114) | Broad environmental multi-omics benchmark [3] | Amplicon targets and matched sample identities require manual curation. |
| `TRI-REAL-04` | HCHS/SOL GOLD gut microbiome study, [ENA ERP117287](https://www.ebi.ac.uk/ena/browser/view/ERP117287), [Qiita 11666](https://qiita.ucsd.edu/study/description/11666) | Supplemental human-associated amplicon stress test | It is not a complete matched WMS/16S/ITS cohort and is unsuitable as the primary cross-assay benchmark. |

### Dataset references

1. Shaffer JP, Carpenter CS, Martino C, et al. A comparison of six DNA
   extraction protocols for 16S, ITS and shotgun metagenomic sequencing of
   microbial communities. *BioTechniques*. 2022;73(1):34–46.
   [doi:10.2144/btn-2022-0032](https://doi.org/10.2144/btn-2022-0032),
   [PubMed 35713407](https://pubmed.ncbi.nlm.nih.gov/35713407/).
2. Regalado J, Lundberg DS, Deusch O, et al. Combining whole-genome shotgun
   sequencing and rRNA gene amplicon analyses to improve detection of
   microbe-microbe interaction networks in plant leaves. *The ISME Journal*.
   2020;14:2116–2130.
   [doi:10.1038/s41396-020-0665-8](https://doi.org/10.1038/s41396-020-0665-8),
   [PubMed 32405027](https://pubmed.ncbi.nlm.nih.gov/32405027/).
3. Shaffer JP, Nothias LF, Thompson LR, et al.; Earth Microbiome Project 500
   Consortium. Standardized multi-omics of Earth's microbiomes reveals
   microbial and metabolite diversity. *Nature Microbiology*.
   2022;7:2128–2150.
   [doi:10.1038/s41564-022-01266-x](https://doi.org/10.1038/s41564-022-01266-x),
   [PubMed 36443458](https://pubmed.ncbi.nlm.nih.gov/36443458/).

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

## Contributing

Bug reports, feature requests, and pull requests are welcome. See
[CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, checks, and
review process, and the [Code of Conduct](CODE_OF_CONDUCT.md). Report security
issues privately as described in [SECURITY.md](SECURITY.md). Notable changes
are recorded in [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use TRIARCHY, please cite it using the metadata in
[CITATION.cff](CITATION.cff) (GitHub's "Cite this repository" button reads it).

## License

TRIARCHY is released under the GNU General Public License v3.0. See
[LICENSE](LICENSE).
