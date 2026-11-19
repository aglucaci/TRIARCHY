# Production inputs

TRIARCHY requires two tab-separated metadata files.

## Sample manifest

Required columns:

```text
library_id
biological_sample_id
specimen_id
subject_id
assay
sample_role
sample_material
environment_category
library_layout
fastq_r1
fastq_r2
extraction_batch
pcr_plate
library_batch
sequencing_run
technical_replicate
```

Valid assays are `WMS`, `16S`, and `ITS`. Valid sample roles are `biological`,
`negative_control`, `positive_control`, and `mock`. Paired libraries require
both FASTQ paths. Relative paths are resolved from the manifest directory.

## Control relationships

Required columns:

```text
library_id
control_library_id
control_stage
relationship
priority
reason
```

Control stages are `field`, `extraction`, `pcr`, `library`, and `sequencing`.
Controls must use the same assay as the production library to which they are
linked.

Validate metadata before launching an expensive run:

```bash
triarchy validate \
  --input /data/study/samples.tsv \
  --paired-controls /data/study/control_relationships.tsv
```
