# Changelog

All notable changes to TRIARCHY are recorded here.

## Unreleased

- Added `bin/triarchy validate` for metadata, FASTQ and database preflight
  checks.
- Added `bin/triarchy databases check` and `bin/triarchy databases install`
  wrappers for reference-manifest validation and checkpointed installation.
- Added a database-manifest example at `conf/databases.example.yml`.
- Added publication-facing release, governance, contribution and citation
  files.
- Added Application Note manuscript source, supplementary outline and Figure 1
  plan.
- Extended benchmark comparison metrics with false-concordance rate,
  false-discordance rate, comparison-rank accuracy and low-abundance authentic
  retention.
- Added Nextflow aliases for `--input` and `--paired_controls`.
- Added test, laptop, workstation, high-memory, Docker, Apptainer, Slurm and
  cloud-template profiles.

## 0.1.0

- Initial executable TRIARCHY workflow skeleton.
- Metadata validation for sample manifests and sample-control relationships.
- WMS, 16S and ITS assay branches.
- Control-aware feature annotation.
- Taxonomic harmonization.
- Cross-assay evidence classification.
- Deterministic benchmark fixtures and repository tests.
