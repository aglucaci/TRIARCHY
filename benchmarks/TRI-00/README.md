# TRI-00 - Metadata contract tests

TRI-00 validates the sample manifest and sample-control graph without requiring
FASTQs. Paths are deliberately synthetic; file existence belongs to the
execution preflight rather than this structural contract test.

## Valid design

The valid fixture includes:

- one specimen with WMS, 16S, and ITS;
- one specimen with WMS and 16S only;
- stage-specific extraction, PCR, library, and sequencing controls;
- a positive control;
- controls shared across biological samples;
- technical replicates;
- two extraction batches and two sequencing runs.

## Invalid cases

Each malformed case has a separate TSV under `invalid/`. The machine-readable
`expected_errors.tsv` maps each fixture to its required diagnostic:

| Case | Expected code |
| --- | --- |
| `duplicate_library_id` | `TRI-META-001` |
| `missing_biological_sample_id` | `TRI-META-002` |
| `invalid_assay` | `TRI-META-003` |
| `paired_missing_r2` | `TRI-META-004` |
| `duplicate_fastq_path` | `TRI-META-005` |
| `unknown_control_id` | `TRI-META-010` |
| `pcr_control_on_wms` | `TRI-META-011` |
| `cross_batch_control` | `TRI-META-014` |
| `empty_required_field` | `TRI-META-008` |
| `inconsistent_subject_mapping` | `TRI-META-009` |
| `conflicting_sample_role` | `TRI-META-016` |
| `cyclic_control_relationship` | `TRI-META-015` |

The automated tests also mutate the full valid graph to confirm each failure in
context. The small invalid fixtures make the contract failures easy to inspect
without duplicating the complete valid manifest.

## Schema notes

`library_layout` is `paired` or `single`. The `fastq_r2` field is required only
for paired libraries. `biological_sample_id` links matched assays;
`technical_replicate` distinguishes repeated libraries without changing the
underlying specimen.
