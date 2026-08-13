# TRIARCHY data contracts

This document defines the minimum machine-readable interfaces for the
control-aware integration layer. TSV files use UTF-8 encoding, a header row,
tab separators, and one value per cell. Missing values are empty, not sentinel
strings such as `NA`.

## Sample manifest

One row represents one sequenced library. A biological specimen assayed by WMS,
16S, and ITS therefore has three rows sharing `biological_sample_id`.

| Field | Required | Meaning |
| --- | --- | --- |
| `library_id` | Yes | Unique sequenced-library identifier |
| `biological_sample_id` | Yes | Identifier shared by matched assays |
| `specimen_id` | Yes | Physical specimen identifier |
| `subject_id` | Yes | Source subject/site identifier; use a documented control value for controls |
| `assay` | Yes | `WMS`, `16S`, or `ITS` |
| `sample_role` | Yes | `biological`, `negative_control`, `positive_control`, or `mock` |
| `sample_material` | Yes | Material or substrate sampled |
| `environment_category` | Yes | Airstream or other prespecified environment |
| `library_layout` | Yes | `paired` or `single` |
| `fastq_r1` | Yes | Forward read path |
| `fastq_r2` | For paired libraries | Reverse read path |
| `extraction_batch` | Yes | Extraction batch identifier |
| `pcr_plate` | For amplicons | PCR plate identifier |
| `library_batch` | Yes | Library-preparation batch identifier |
| `sequencing_run` | Yes | Sequencing run identifier |
| `technical_replicate` | Yes | Positive integer within specimen and assay |
| `biomass_measure` | No | Quantitative biomass proxy |
| `biomass_unit` | If measured | Unit for `biomass_measure` |

Identifiers must be stable strings. Paths must exist and assay values are
case-sensitive.

## Sample-control edges

Controls are represented as a graph rather than a single generic `control_id`.
One row links one biological library to one experimentally relevant control.
A sample can have multiple controls at the same or different stages.

| Field | Required | Meaning |
| --- | --- | --- |
| `library_id` | Yes | Biological library |
| `control_library_id` | Yes | Linked control library |
| `control_stage` | Yes | `field`, `extraction`, `pcr`, `library`, or `sequencing` |
| `relationship` | Yes | `paired`, `same_batch`, or `run_level` |
| `priority` | Yes | Positive integer; 1 is the closest control |
| `reason` | Yes | Human-readable experimental justification |

Both identifiers must exist in the sample manifest, must use the same assay,
and the target must have a control `sample_role`.

## Taxonomy crosswalk

Every source assignment remains traceable after harmonization.

| Field | Required |
| --- | --- |
| `triarchy_taxon_id` | Yes |
| `canonical_name` | Yes |
| `canonical_rank` | Yes |
| `canonical_lineage` | Yes |
| `source_database` | Yes |
| `source_taxon_id` | Yes |
| `original_name` | Yes |
| `original_rank` | Yes |
| `mapping_method` | Yes |
| `mapping_confidence` | Yes |
| `database_version` | Yes |

`mapping_confidence` is numeric from 0 to 1. Comparisons use the lowest
canonical rank supported by both assays; a lower-resolution match is recorded
as taxonomic ambiguity rather than species-level disagreement.

## Control-annotated feature profile

| Field | Meaning |
| --- | --- |
| `library_id` | Sequenced library |
| `assay` | WMS, 16S, or ITS |
| `triarchy_taxon_id` | Harmonized feature identifier |
| `native_count` | Original assay count |
| `native_normalized_value` | Within-assay normalized value |
| `normalization_method` | Method used within this assay |
| `control_prevalence` | Prevalence among linked negative controls |
| `biological_prevalence` | Prevalence among relevant biological samples |
| `control_abundance_ratio` | Linked-control to biological abundance evidence |
| `batch_localization_score` | Concentration of signal in a laboratory batch |
| `contamination_probability` | Bounded linked-control evidence score from 0 to 1 |
| `associated_control_types` | Semicolon-separated control stages |
| `associated_batches` | Semicolon-separated batch identifiers |
| `decision` | `retain`, `flag`, `exclude`, or `indeterminate` |
| `decision_reason` | Human-readable rule/model explanation |

Unfiltered values are never overwritten by a contamination decision.

## Integrated taxon evidence

One row represents a biological sample and harmonized taxon.

| Field | Meaning |
| --- | --- |
| `biological_sample_id` | Matched specimen identifier |
| `triarchy_taxon_id` | Harmonized taxon |
| `comparison_rank` | Lowest valid rank shared by compared assays |
| `*_measured` | Whether that assay was performed for this biological sample |
| `*_detected` | Whether that assay produced a compatible taxon call |
| `*_library_depth` | Assay-native total used for power estimation |
| `*_detection_power` | Estimated probability of detection under the configured model |
| `*_reference_coverage` | Configured or taxon-specific reference representation |
| `wms_evidence_score` | WMS detection, abundance-rank, resolution, and mapping support |
| `16s_evidence_score` | 16S detection, abundance-rank, resolution, and mapping support |
| `its_evidence_score` | ITS detection, abundance-rank, resolution, and mapping support |
| `abundance_support` | Rank-based within-assay abundance support |
| `taxonomic_resolution_support` | Support implied by the shared comparison rank |
| `mapping_support` | Minimum mapping confidence among supporting calls |
| `assay_reliability` | Kingdom-specific assay weights used |
| `control_evidence` | Summary contamination evidence |
| `field_contamination_score` | False-detection evidence linked to field controls |
| `extraction_contamination_score` | False-detection evidence linked to extraction controls |
| `pcr_contamination_score` | False-detection evidence linked to PCR controls |
| `library_contamination_score` | False-detection evidence linked to library controls |
| `run_contamination_score` | False-detection evidence linked to sequencing/run controls |
| `other_control_contamination_score` | False-detection evidence linked to other declared controls |
| `probable_contamination_source` | Stage or batch with the strongest contamination evidence |
| `taxonomic_compatibility_score` | Compatibility score across kingdom-appropriate assay calls |
| `taxonomic_compatibility_status` | Exact/canonical, hierarchical, single-assay, ambiguous, or incompatible compatibility state |
| `biological_authenticity` | Stage-I posterior-style estimate that the taxon is biologically authentic |
| `authenticity_ci_low` | Lower approximate uncertainty bound for `biological_authenticity` |
| `authenticity_ci_high` | Upper approximate uncertainty bound for `biological_authenticity` |
| `authenticity_uncertainty` | Width of the approximate authenticity interval |
| `authenticity_state` | `biologically_supported`, `not_supported`, or `indeterminate` |
| `cross_assay_status` | Stage-II parent explanation class |
| `presence_score` | Legacy alias mirroring `biological_authenticity` for older consumers |
| `probability_authentic` | Legacy calibrated-score column retained for compatibility |
| `calibration_status` | Scope of calibration, including explicit non-empirical status |
| `biological_state` | `present`, `absent`, or `indeterminate` |
| `sample_evidence_rank` | Descending score rank within a biological sample |
| `global_evidence_rank` | Descending score rank across the run |
| `classification` | Detailed evidence class retained for compatibility and audit trails |
| `discordance_reason` | Resolution, scope, control, direction, or depth explanation |
| `model_version` | Version of the decision policy or model |

An absent assay is different from a measured non-detection and must remain
distinguishable in serialized output.

The valid `cross_assay_status` values are `concordant`, `assay_limited`,
`taxonomy_limited`, `contamination_associated`, `genuine_discordance`, and
`indeterminate`.

The valid detailed `classification` values are `concordant_detection`,
`single_assay_credible`, `assay_limited_nondetection`,
`rank_limited_agreement`, `reference_limited`, `probable_contaminant`,
`cross_contamination_candidate`, `genuine_discordance`, and `indeterminate`.

`integrated/evidence_model.json` serializes the exact policy, assay-performance
parameters, calibration status, and model version used for the run.

## Cross-assay concordance

`integrated/cross_assay_concordance.tsv` contains one sample, kingdom, and
biologically valid assay pair per row. It reports assay-specific taxon counts,
shared and union counts, presence Jaccard, Spearman rank correlation, and a
rank-agreement score.

## Sample confidence

```text
biological_sample_id
wms_qc_status
16s_qc_status
its_qc_status
control_status
cross_assay_status
recommended_use
failure_reasons
```

## Run provenance

Each run records a JSON document containing:

- workflow version and Git commit;
- command and complete resolved configuration;
- container digests and software versions;
- database names, releases, checksums, and retrieval dates;
- checksums for manifests and sample-control edges;
- start and completion timestamps;
- task-level CPU, memory, storage, and runtime usage.

The provenance record must be sufficient to identify inputs and environments;
it must not embed credentials or protected sample metadata.
