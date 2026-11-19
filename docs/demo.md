# Production demo

TRIARCHY includes a deterministic production demonstration built from
version-pinned public NCBI nucleotide sequences. It contains 21 WMS, 16S, and ITS
libraries spanning three matched biological samples, stage-specific negative
controls, a mixed positive control, human mitochondrial reads, and designed
cross-assay edge cases. It uses the normal production workflow and installed
reference databases—there are no stub tools or precomputed taxonomic calls.

The demo is a synthetic benchmark, not a biological cohort. It generates reads
from the following public reference records and records their exact sequence
checksums in `data/production_demo/metadata/references.tsv`:

| Target | Assay role | Version-pinned NCBI records |
| --- | --- | --- |
| *Escherichia coli* | WMS and 16S | [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3), [NR_114042.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_114042.1) |
| *Pseudomonas aeruginosa* | WMS and 16S | [NC_002516.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002516.2), [NR_026078.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_026078.1) |
| *Staphylococcus epidermidis* | WMS and 16S | [NC_004461.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_004461.1), [NR_113957.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_113957.1) |
| *Saccharomyces cerevisiae* | WMS and ITS | [NC_001133.9](https://www.ncbi.nlm.nih.gov/nuccore/NC_001133.9), [NC_001144.5](https://www.ncbi.nlm.nih.gov/nuccore/NC_001144.5) |
| *Candida albicans* | ITS | [FJ662402.1](https://www.ncbi.nlm.nih.gov/nuccore/FJ662402.1) |
| *Homo sapiens* mitochondrion | host-removal challenge | [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) |

The generated mixture contains 317,000 paired reads at the default scale. Its
truth table intentionally covers concordant bacterial and fungal detections,
assay-limited nondetection, low-biomass contamination, positive controls, and
approximately 15% mitochondrial host reads in one WMS library. Because the
benchmark is constructed from isolated reference sequences, it tests workflow
execution and evidence rules; it must not be used to estimate performance on
natural communities.

Build the dataset and run the complete production pipeline:

```bash
mamba activate TRIARCHY
triarchy demo --profile laptop --resume
```

The command downloads the small source references from NCBI, generates exactly
317,000 deterministic paired reads, validates the manifests and FASTQs, runs
TRIARCHY, and evaluates the resulting evidence and host-removal metrics against
the designed truth table. Outputs are written to `results/production_demo`.

To generate and validate the inputs without starting Nextflow:

```bash
triarchy demo --build-only
```

Use `--offline` after the NCBI references have been cached, or `--force` to
regenerate all reads with the selected `--seed`, `--scale`, and `--error-rate`.

On a SLURM HPC, keep the repository, production references, inputs, results, and
Nextflow work directory on shared storage visible to compute nodes:

```bash
mamba activate TRIARCHY
triarchy demo \
  --profile slurm \
  --site-config /shared/project/triarchy-site.config \
  --work-dir /shared/scratch/$USER/triarchy-work
```

The built-in SLURM profile requests 1 CPU/2 GB for low, 4 CPU/16 GB for medium,
and 12 CPU/96 GB for high-resource tasks. The optional site configuration is
where a cluster-specific partition, allocation, QoS, module, or scratch policy
should be declared. The activated `TRIARCHY` environment and database paths
must also be available at the same paths from submitted jobs.
