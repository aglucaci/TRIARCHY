# TRI-16 - Missing-modality and failure challenge

Graceful integration of missing, failed, low-depth, and unreliable assays.

The deterministic fixture contains 30 libraries, 21 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Showcase expected outputs live under `expected_outputs/`; a precomputed narrative HTML report lives under `report/`; mini functional demo references live under `references/mini_function_catalog.fasta` and `references/mini_function_mapping.tsv`.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-16
```
