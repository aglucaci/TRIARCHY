# TRIARCHY benchmark suite

TRIARCHY uses multiple datasets because metadata validation, assay processing,
contamination modeling, taxonomic harmonization, integration, and ecological
recovery require different truth designs.

## Tiers

| Tier | IDs | Purpose |
| --- | --- | --- |
| Software/CI | TRI-00–TRI-02 | Contracts, smoke tests, and isolated assay branches |
| Compact bundle | TRI-01, TRI-03, TRI-05–TRI-07 | Five-dataset, 97-library CI and resource-profile bundle |
| Synthetic truth | TRI-10–TRI-19 | Accuracy, sensitivity, contamination, ecology, harmonization, failures, carryover, and repeatability |
| Physical mocks | TRI-20–TRI-21 | Extraction, amplification, library, and laboratory performance |
| Environmental | TRI-30 | Publication-grade matched environmental application |

The machine-readable inventory is [registry.tsv](registry.tsv). `status` has
one of four values:

- `implemented`: data and executable checks are present;
- `implemented_design_surrogate`: the complete executable metadata and truth
  contract is present, but empirical laboratory or field data are still
  required before biological claims can be made;
- `scaffolded`: the contract exists but data generation or acquisition remains;
- `planned`: design only.

## Implemented datasets

### TRI-00 — Metadata contract tests

[TRI-00](TRI-00/README.md) contains a valid multi-assay sample–control graph and
separate malformed cases. The validator returns stable `TRI-META-*` codes so CI
can test both acceptance and failure behavior.

Run it with:

```bash
python3 scripts/validate_benchmark_metadata.py \
  benchmarks/TRI-00/valid/samples.tsv \
  benchmarks/TRI-00/valid/control_relationships.tsv

python3 -m unittest discover -s tests -p 'test_*.py'
```

No external Python packages are required.

### Runnable CI, synthetic, and design-surrogate fixtures

Twelve deterministic datasets implement the inventory while remaining small
enough to keep in the repository:

| Dataset | Libraries | Scenario |
| --- | ---: | --- |
| TRI-01 | 12 | End-to-end smoke test with biological, negative, and positive inputs |
| TRI-02 | 18 | Independent WMS, 16S, and ITS branch stressors |
| TRI-03 | 18 | Compact balanced bacterial, archaeal, and fungal truth community |
| TRI-05 | 24 | Compact biomass and contamination gradient |
| TRI-06 | 18 | Compact synonym, rank, historical-name, and held-out-reference challenge |
| TRI-07 | 25 | Compact missing-modality, low-depth, host, empty, and failure challenge |
| TRI-10 | 21 | Even bacterial, archaeal, and fungal community |
| TRI-11 | 24 | Six-order log-abundance sensitivity series |
| TRI-12 | 112 | Two-batch, six-level low-biomass contamination ladder |
| TRI-13 | 75 | Four-habitat structured environmental gradient |
| TRI-14 | 48 | Synonym, rank, close-relative, and held-out harmonization |
| TRI-15 | 97 | Balanced and deliberately confounded batch layouts |
| TRI-16 | 30 | Missing modalities and explicit library failures |
| TRI-17 | 36 | Depth, primer recovery, copy-number, and held-out-reference factorial |
| TRI-18 | 54 | Sequencing-run carryover and index leakage with stage-linked controls |
| TRI-19 | 84 | Technical and extraction repeatability across three batches |
| TRI-20 | 45 | Commercial-standard executable design surrogate |
| TRI-21 | 114 | Custom mixed-kingdom mock design surrogate |
| TRI-30 | 192 | Balanced 48-sample built-environment design surrogate |

Together they contain 1,047 libraries. Every dataset includes actual
gzip-compressed FASTQs, valid metadata, sample–control edges, the complete truth
package, checksums, and generation provenance. TRI-20, TRI-21, and TRI-30 are
explicitly labeled design surrogates: they make acquisition plans and analysis
contracts executable but do not replace physical mocks or field data.

Regenerate and validate all fixtures:

```bash
scripts/run_benchmark_suite.sh
```

Exercise every dataset through the complete deterministic Nextflow graph:

```bash
scripts/run_benchmark_suite.sh --nextflow
```

Generate one dataset only:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-12
```

### Compact five-dataset resource bundle

The versioned bundle definition at
`benchmarks/bundles/TRIARCHY_benchmark_bundle_v0.1.json` selects TRI-01,
TRI-03, TRI-05, TRI-06, and TRI-07: five datasets and 97 modality-specific
libraries. Generate them with:

```bash
python3 scripts/generate_benchmark_datasets.py \
  --dataset TRI-01 --dataset TRI-03 --dataset TRI-05 \
  --dataset TRI-06 --dataset TRI-07
```

Measure orchestration wall time, process-tree peak memory, task metrics, and a
rolling ETA with:

```bash
.tools/env/bin/python scripts/benchmark_resources.py --mode stub
```

Stub measurements do not size production bioinformatics work. After installing
production references and supplying the required parameters, repeat with
`--mode production`. Resource reports are written as CSV and JSON under
`output/benchmarks/`; individual logs and workflow results are retained under
`results/resource-benchmarks/`.

## Dataset contract

Every dataset directory must eventually contain:

```text
README.md
metadata/
├── samples.tsv
└── control_relationships.tsv
truth/
├── community_members.tsv
├── sample_composition_genomic_dna.tsv
├── sample_composition_genome_copy.tsv
├── sample_composition_16s.tsv
├── sample_composition_its.tsv
├── read_origin_wms.tsv.gz
├── true_16s_asvs.fasta
├── true_its_asvs.fasta
├── taxonomic_crosswalk.tsv
├── contaminants.tsv
├── control_relationships.tsv
├── expected_sample_qc.tsv
├── expected_taxon_status.tsv
├── expected_pathways.tsv
└── generation_manifest.json
```

Physical and environmental datasets may omit inapplicable synthetic read-origin
files, but must document why. Restricted raw reads may be referenced by stable
accession rather than stored in Git.

## Truth scales

Truth is assay appropriate. Genomic-DNA abundance, genome-copy abundance, 16S
template abundance, ITS template abundance, and observed read abundance are
distinct quantities. A benchmark must never reuse one expected percentage
vector for all three assays.

## Minimum publication package

The initial paper target is TRI-01, TRI-10, TRI-11, TRI-12, TRI-14, TRI-20, and
TRI-30. TRI-12 is the primary methodological benchmark because it tests whether
paired-control integration preserves authentic low-abundance taxa while
removing known contaminants.

## Versioning and reproducibility

Synthetic datasets must record genome accessions and checksums, simulator
versions, random seeds, primers, read and insert lengths, error models,
copy-number assumptions, contaminant sources, and frozen reference-database
versions in `generation_manifest.json`.

The tiny deterministic CI FASTQs are committed directly. Publication-scale
FASTQs and reference genomes belong in a versioned data release or Git LFS, not
ordinary Git history.
