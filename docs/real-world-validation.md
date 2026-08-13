# Real-World Validation Datasets

TRIARCHY includes a downloader for public empirical datasets that can be used
after the deterministic benchmark suite. The downloader uses ENA run metadata
to build reproducible file manifests under `data/real_world/`; it downloads
FASTQs only when `--download` is provided.

Recommended datasets:

| ID | Source | Use |
| --- | --- | --- |
| `TRI-REAL-01` | ENA `ERP124610`, Qiita `12201` | Primary mixed bacterial/fungal mock and extraction benchmark with blanks and dilution structure |
| `TRI-REAL-02` | ENA `PRJEB31530` | Secondary matched environmental WMS/16S/ITS complete-case validation |
| `TRI-REAL-03` | ENA `PRJEB42019`, Qiita `13114` | Optional multi-environment breadth check |
| `TRI-REAL-04` | ENA `ERP117287`, Qiita `11666` | Optional supplemental stress test; not environmental |

Plan a subset and write manifests:

```bash
python scripts/download_real_world_datasets.py --dataset TRI-REAL-01 --max-groups 30
```

Download the selected FASTQs after inspecting the manifest:

```bash
python scripts/download_real_world_datasets.py --dataset TRI-REAL-01 --max-groups 30 --download --verify-md5
```

Downloads print one progress line per file with transferred bytes, percent,
speed and ETA, plus an overall summary. For more frequent updates, use
`--progress-interval 0.5`; for batch logs, use `--quiet`.

Plan all configured datasets:

```bash
python scripts/download_real_world_datasets.py --dataset all
```

Each dataset directory contains:

```text
data/real_world/<dataset_id>/
|-- metadata/
|   |-- ena_runs.tsv
|   |-- selected_runs.tsv
|   |-- download_manifest.tsv
|   |-- samples.tsv
|   |-- control_relationships.tsv
|   `-- dataset.json
`-- reads/
```

`ena_runs.tsv` is the full ENA read-run result. `selected_runs.tsv` adds
TRIARCHY helper columns for inferred assay, sample group, and sample role.
When ENA does not expose enough indexed metadata to distinguish 16S from ITS,
the helper assay is recorded as `AMPLICON` rather than guessed.
`download_manifest.tsv` records one FASTQ per row with source URL, local path,
MD5, and byte count. `dataset.json` records the query URL, selection summary,
and download summary.

The generated `samples.tsv` and `control_relationships.tsv` are draft TRIARCHY
contracts for workflow testing. When ENA does not expose enough indexed
metadata to distinguish 16S from ITS, unresolved `AMPLICON` rows are mapped to
`16S` by default; use `--amplicon-assay ITS` for an ITS-oriented draft. The
control relationship file is intentionally empty because ENA run metadata alone
does not prove stage-specific sample-control edges. Before biological
validation, curate both files from the original study metadata so controls,
batches, primers, sample roles, and matched biological sample IDs are explicit.

Run the complete TRIARCHY graph in stub mode with the draft metadata:

```bash
./run_triarchy.sh --stub \
  --samples data/real_world/TRI-REAL-01/metadata/samples.tsv \
  --controls data/real_world/TRI-REAL-01/metadata/control_relationships.tsv \
  --outdir results/TRI-REAL-01-stub
```

FASTQ downloads can be large. Keep `data/` outside Git and use a filesystem
with enough space for reads, Nextflow `work/`, and final results.
