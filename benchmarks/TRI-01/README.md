# TRI-01 — Tiny end-to-end smoke test

Complete three-assay execution with biological samples and controls.

The deterministic fixture contains 12 libraries, 6 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-01
```
