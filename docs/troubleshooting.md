# Troubleshooting

Inspect a failed Nextflow task under its `work/<hash>/` directory:

```text
.command.sh
.command.err
.command.log
```

Confirm the named environment and production references before rerunning:

```bash
mamba activate TRIARCHY
nextflow config -flat | grep -E 'host_index|kraken2_db|humann|taxonomy'
```

On laptops, prefer `-profile laptop`, close memory-intensive applications, and
place the Nextflow work directory and full references on fast external storage.
