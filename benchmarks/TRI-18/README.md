# TRI-18 — Sequencing carryover and index leakage

Separate run-local carryover from extraction and assay-specific contamination using linked controls.

The deterministic fixture contains 54 libraries, 72 batch-local control links, paired or explicitly single-end FASTQs, assay-specific abundance scales, read origins, ASV truth, artificial references, precomputed profiles, machine-readable assertions, taxonomy mappings, contamination truth, QC expectations, pathways, checksums, and provenance.

Regenerate from the repository root:

```bash
python3 scripts/generate_benchmark_datasets.py --dataset TRI-18
```
