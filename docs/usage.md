# Run production

Run the production workflow with real FASTQs:

```bash
triarchy run \
  --profile local \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --outdir results/study
```

Assay primers and the ITS target region are maintained centrally in
`conf/primers.config`; users do not pass them on the command line. Update that
controlled configuration when a production protocol changes.

Resume an interrupted run with `--resume`:

```bash
triarchy run \
  --profile local \
  --resume \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --outdir results/study
```

Available production resource profiles:

| Profile | Intended system |
| --- | --- |
| `laptop` | Constrained local production; queue size 2, high tasks capped at 8 GB |
| `local` | Default local executor resources |
| `workstation` | Local workstation with at least 32 GB RAM |
| `high_memory` | High-memory host; high tasks request 64 GB |
| `slurm` | SLURM cluster execution |
| `cloud` | Cloud executor configuration |
| `docker` | Production Docker execution |
| `apptainer` | Production Apptainer execution |

## Container profiles

With `--profile docker`, `apptainer`, or `singularity`, `triarchy run` reads
the reference root recorded by the installer (`references/.database-dir`) and
passes it as `--reference_root`, so the databases are mounted into each task
container. When calling `nextflow run` directly, pass `--reference_root`
yourself:

```bash
nextflow run main.nf -profile docker -c conf/references.config \
  --reference_root /path/to/production/references \
  --samples samples.tsv --controls control_relationships.tsv --outdir results/study
```

Container profiles do not set an executor, so they combine with scheduler
profiles. For containers on SLURM, call Nextflow directly with
`-profile slurm,apptainer` (`triarchy run` accepts a single profile).

## Parameters

Every pipeline parameter is described in [parameters.md](parameters.md),
generated from `nextflow_schema.json`.

## Version

```bash
triarchy --version
```

Each run records the TRIARCHY version, Nextflow version, and git commit in
`pipeline_info/provenance.json`.
