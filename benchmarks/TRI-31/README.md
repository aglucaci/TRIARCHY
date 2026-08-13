# TRI-31 - Thirty-minute TRIARCHY showcase dataset

Home-computer showcase with 12 biological samples, seven controls, host DNA, reagent contamination, batch contamination, primer bias, taxonomic ambiguity, fungal reference limitation, pathways, expected figures, and expected tables.

This is the repository-sized implementation of the TRIARCHY 30-minute showcase design. It includes the full 19-sample, 51-library metadata design and expected-output bundle; the committed FASTQs are scaled down for Git and CI while the full-scale target is 2,526,000 paired-end reads. The deterministic fixture contains 51 libraries, 108 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Showcase expected outputs live under `expected_outputs/`; a precomputed narrative HTML report lives under `report/`; mini functional demo references live under `references/mini_function_catalog.fasta` and `references/mini_function_mapping.tsv`.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-31
```
