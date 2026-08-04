# TRI-20 — Commercial mixed-community standards

Executable design surrogate for even DNA, whole-cell, and log DNA standards.

This is an executable design surrogate, not empirical laboratory or field data. The deterministic fixture contains 45 libraries, 27 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-20
```
