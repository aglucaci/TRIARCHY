# TRI-17 — Depth, primer, and reference factorial

Jointly test depth-limited detection, marker recovery bias, copy number, and held-out references.

The deterministic fixture contains 36 libraries, 30 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-17
```
