# TRI-03 — Balanced mixed-kingdom truth community

Replicated bacterial, archaeal, and fungal truth with negative and positive controls.

The deterministic fixture contains 18 libraries, 9 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-03
```
