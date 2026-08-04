# TRI-21 — Custom environmental mixed-kingdom mock

Executable design surrogate spanning profile, biomass, and kingdom diversity.

This is an executable design surrogate, not empirical laboratory or field data. The deterministic fixture contains 114 libraries, 108 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-21
```
