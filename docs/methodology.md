# TRIARCHY methodology

## Scope

TRIARCHY jointly interprets matched WMS, 16S, and ITS data while preserving the
measurement scale and biological scope of each assay. It integrates evidence,
not raw relative abundance.

## 1. Sample-control graph

The experimental design is a bipartite graph connecting biological libraries
to field, extraction, PCR, library, and sequencing controls. Edges carry the
laboratory stage, batch relationship, and priority. Contamination evidence for
a sample is estimated from its linked controls; study-wide controls are used
only when the design explicitly declares a run-level relationship.

Validation rejects:

- duplicate library identifiers;
- unknown assays or sample roles;
- missing required batch fields;
- control edges whose identifiers do not exist;
- cross-assay control edges;
- biological samples incorrectly used as negative controls;
- paired assays whose biological sample metadata conflict.

## 2. Assay-native profiles

Counts and normalization remain assay specific. Raw and normalized profiles are
immutable inputs to the integration layer.

- 16S describes bacterial and archaeal marker-gene evidence.
- ITS describes fungal marker-gene evidence.
- WMS describes taxonomic evidence and functional potential.

TRIARCHY does not average relative abundances across assays. Cross-assay
visualization may use prevalence, rank abundance, within-assay CLR values,
detection probabilities, or standardized effect sizes with the transformation
recorded.

## 3. Control-aware feature annotation

For each sample-taxon-assay combination, the implementation calculates three
transparent signals:

1. **Prevalence:** whether the feature occurs more often in linked negative
   controls than in biological samples.
2. **Abundance:** whether the feature is unexpectedly enriched in the closest
   paired blank or relevant batch controls.
3. **Batch localization:** whether detections concentrate in one extraction,
   PCR, library, or sequencing batch.

These signals produce a bounded contamination-evidence score and a versioned
policy decision. The score is not described as a calibrated probability.
Every decision includes its inputs, policy version, and explanation.
Profiles are emitted in unfiltered, control-annotated, and high-confidence
forms so filtering remains auditable.

## 4. Taxonomic harmonization

Kraken2/Bracken, marker-gene, and UNITE assignments are not merged by display
name. Each source identifier maps into a canonical representation with database
version, method, and confidence.

Assay comparisons roll both calls up to the lowest rank valid in both. For
example, a species-level WMS assignment and genus-level 16S assignment are
compared at genus. Agreement at genus is classified as rank-limited rather than
as species-level concordance or discordance.

## 5. Formal evidence inference

Cross-assay comparisons respect biological scope:

- bacteria and archaea: WMS <-> 16S;
- fungi: WMS <-> ITS;
- functional potential: WMS only.

The TRIARCHY evidence model is two-stage. Stage I estimates
`biological_authenticity`, a posterior-style evidence estimate that a taxon is
biologically authentic in a sample. Stage II assigns `cross_assay_status`, the
explanation for agreement or disagreement among the kingdom-appropriate assays.

The model treats every sample-taxon-assay combination as detected,
measured-but-not-detected, or not measured. It estimates assay-specific
detection power from library depth, assay reliability, quality state, kingdom,
primer compatibility when supplied, and reference coverage. This prevents a
missing assay or underpowered fungal WMS library from being interpreted as
biological absence.

The model emits a parent `cross_assay_status` for summaries:

- `concordant`;
- `assay_limited`;
- `taxonomy_limited`;
- `contamination_associated`;
- `genuine_discordance`;
- `indeterminate`.

It also retains nine detailed explanatory classes in `classification`:

- `concordant_detection`;
- `single_assay_credible`;
- `assay_limited_nondetection`;
- `rank_limited_agreement`;
- `reference_limited`;
- `cross_contamination_candidate`;
- `genuine_discordance`;
- `probable_contaminant`;
- `indeterminate`.

The Stage-I `biological_authenticity` estimate starts from a versioned prior
and updates log odds using positive detections, technically adequate
nondetections, mapping support, taxonomic compatibility, reference coverage,
and linked-control contamination evidence. It is a posterior-style evidence
estimate under the configured policy, not an empirically calibrated probability
of presence. Approximate uncertainty bounds are serialized as
`authenticity_ci_low`, `authenticity_ci_high`, and `authenticity_uncertainty`;
wide intervals are assigned `authenticity_state = indeterminate`.
Per-assay positive evidence is:

```text
e(s,t,m) = detection_support
           * sqrt(abundance_support)
           * resolution_support
           * mapping_support
           * reliability
           * primer_compatibility
           * (1 - contamination_probability)
```

Detection power for a measured assay is:

```text
power(s,t,m) = reliability(m,t)
               * [1 - exp(-library_depth / minimum_depth)]
               * reference_coverage(m,t)
               * quality_factor(s,m)
```

False-detection evidence is tracked by experimental stage when linked controls
are present: field, extraction, PCR, library, sequencing/run, other controls,
and batch localization. The strongest source is serialized as
`probable_contamination_source`.

Defaults live in `conf/assay_performance.tsv`; thresholds, priors, model
version, and calibration metadata live in `conf/evidence_policy.json`. Both
files are inputs to every run and are copied into `evidence_model.json`.

Evidence rows receive `sample_evidence_rank` and `global_evidence_rank`.
Comparisons use the least-resolved supported rank and also report
`taxonomic_compatibility_score` and `taxonomic_compatibility_status`. Calls
sharing a canonical taxon but differing in resolution are
`rank_limited_agreement`; explicitly incompatible assignments are
`genuine_discordance` when both assay contexts are technically adequate.

The repository includes an isotonic calibration command and a truth-based
evaluation command. The shipped mapping is explicitly
`synthetic_only_not_empirically_calibrated`; it must not be presented as
empirically calibrated confidence. Independence is a simplifying assumption. A
latent-community model remains a future method, not current behavior.

Benchmark comparisons are serialized by
`scripts/compare_evidence_methods.py`. It accepts named outputs from independent
assay analysis, naive merging, shared-taxonomy analysis, Greengenes2,
per-assay decontam, global controls, Com-2seq, and full TRIARCHY, and reports
authentic-taxon precision/recall/F1, contaminant precision/recall,
explanatory-status accuracy, and Brier score against the same truth table.

## 6. Diversity and functional interpretation

Diversity is calculated separately for:

- 16S bacterial/archaeal profiles;
- ITS fungal profiles;
- WMS taxonomic profiles;
- WMS functional profiles.

Cross-assay concordance compares distance matrices, sample clustering, and
biological effect sizes. A joint diversity value requires a separately
validated common or latent representation.

The executable concordance table reports kingdom-aware presence Jaccard,
Spearman rank correlation on shared taxa, and rank agreement for WMS-16S
bacterial/archaeal and WMS-ITS fungal comparisons. These compare detection and
ordering, never raw abundance across assays.

HUMAnN pathways, coverage, contributing taxa, and unclassified fractions remain
distinct from taxonomic confidence. 16S and ITS may support the presence of
organisms associated with a pathway; they do not independently validate that
functional potential.

## 7. Analytical report

The report must answer:

1. Did sequencing and assay processing succeed?
2. Which samples are suitable for downstream use, and why?
3. Which taxa have reliable evidence?
4. Where do assays disagree, and can scope, depth, taxonomy, or controls explain
   the disagreement?
5. Which ecological conclusions are robust before and after control annotation
   and across assay-specific analyses?

## 8. Validation plan

### Minimum benchmark

- in silico mixed bacterial, archaeal, and fungal communities;
- matched tri-assay mixed-kingdom mock community;
- biomass dilution series and negative controls;
- real environmental case study;
- ablations for naive merging, harmonization, global versus paired controls,
  reliability weighting, and modality removal.

### Metrics

- taxonomic precision, recall, F1, false-discovery rate, and resolved rank;
- Spearman correlation, Aitchison distance, and compositional error;
- contaminant precision/recall and biological taxa incorrectly removed;
- rank-aware Jaccard overlap and concordance of ecological effect sizes;
- replicate consistency, PERMANOVA effect size, and downsampling stability;
- CPU hours, peak memory, storage, runtime, and scaling by sample count.

### Release gate

TRIARCHY should not claim calibrated integrated confidence until:

1. expected outputs are tested on a small CI dataset;
2. truth-known mock results establish score calibration;
3. paired-control performance is compared with global and per-assay filtering;
4. thresholds and database releases are versioned;
5. local and containerized runs reproduce checksummed outputs.

## 9. Implementation roadmap

1. Restore or add executable workflow source and tests to version control.
2. Implement strict validation for the two input contracts.
3. Emit immutable assay-native and control-annotated profiles.
4. Add versioned taxonomy crosswalk generation.
5. Implement two-stage authenticity/explanation inference with complete decision reasons.
6. Add assay-specific diversity and ecological concordance.
7. Build the five-question analytical report.
8. Add containers, reference manifests, CI data, provenance, and benchmarks.
9. Calibrate the evidence score; evaluate a latent-community model separately.
