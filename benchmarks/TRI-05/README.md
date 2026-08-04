# TRI-05 — Compact contamination gradient

Five-bundle low-biomass and batch-specific contamination challenge with assay-specific truth.

The deterministic fixture contains 24 libraries, 15 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-05
```
