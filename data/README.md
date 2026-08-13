# TRIARCHY Data Folder

Downloaded real-world validation data lives here during local development.
Large FASTQ files, ENA metadata dumps, and generated manifests are intentionally
not committed to Git.

Use:

```bash
python scripts/download_real_world_datasets.py --dataset TRI-REAL-01
python scripts/download_real_world_datasets.py --dataset TRI-REAL-01 --download --verify-md5
```

Download mode prints per-file progress, speed and ETA by default.

Downloader output is written under `data/real_world/<dataset_id>/`.
