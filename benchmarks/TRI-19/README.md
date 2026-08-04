# TRI-19 — Technical and extraction repeatability

Measure within-run, between-run, and extraction-batch stability of evidence classifications.

The deterministic fixture contains 84 libraries, 72 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-19
```
