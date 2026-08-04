# TRIARCHY compact benchmark bundle v0.1

This bundle selects TRI-01, TRI-03, TRI-05, TRI-06, and TRI-07: five datasets,
34 underlying specimens, and 97 modality-specific libraries.

Each dataset contains deterministic gzip-compressed FASTQs, metadata and
paired-control manifests, artificial bacterial/archaeal/fungal/host reference
FASTAs, per-read origins, assay-specific abundance truth, a taxonomic
crosswalk, precomputed profiles, expected assertions, SHA-256 checksums, and
generation provenance. These repository-scale data validate software and CI
behavior; they do not replace publication benchmarks built from real genomes,
frozen production databases, physical mocks, and field samples.

Regenerate the bundle datasets:

```bash
python3 scripts/generate_benchmark_datasets.py \
  --dataset TRI-01 --dataset TRI-03 --dataset TRI-05 \
  --dataset TRI-06 --dataset TRI-07
```

Profile wall time, memory-time, throughput, and rolling ETA:

```bash
python3 scripts/benchmark_resources.py --mode stub
```

Use `--mode production` only after the production tools, references, primers,
and database parameters are configured. Stub results measure orchestration and
report generation, not biological processing cost.
