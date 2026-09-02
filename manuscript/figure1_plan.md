# Figure 1 plan

Use a single four-panel figure for the Application Note.

Generated draft asset:

```text
manuscript/figures/figure1_composite.svg
```

Convert the final SVG to `manuscript/figures/figure1_composite.pdf` before
compiling the submission manuscript.

## Panel A: Workflow architecture (about 20% of figure)

Manifests plus paired controls -> sample-control graph -> WMS, 16S and ITS
branches -> taxonomic harmonization -> two-stage authenticity/explanation
model -> authenticity estimates, cross-assay explanations, ecology and report.

## Panel B: Evidence logic (about 25% of figure)

Show Stage I and Stage II outputs. Stage I reports `biological_authenticity`,
approximate uncertainty bounds, and `authenticity_state`. Stage II reports six
parent `cross_assay_status` explanations, with detailed statuses in the
supplement:

- Concordant: concordant detection.
- Assay-limited: single-assay credible, assay-limited nondetection.
- Taxonomy-limited: rank-limited agreement, reference-limited.
- Contamination-associated: probable contaminant, cross-contamination candidate.
- Genuine discordance: technically adequate incompatible evidence.
- Indeterminate: unresolved evidence or wide uncertainty.

## Panel C: Benchmark comparison (about 40% of figure)

Plot authentic-taxon F1 and contaminant F1 for:

- independent assay analysis;
- naive genus concatenation;
- shared taxonomy only;
- per-assay contamination filtering;
- global negative-control filtering;
- TRIARCHY without detection-power modelling;
- full TRIARCHY.

Panel C should be the dominant panel. It should report measured authentic-taxon
F1, contaminant F1, status accuracy, and false-concordance or
false-discordance rates for the frozen release. The generated panel is a
template until `scripts/compare_evidence_methods.py` has been run on frozen
release outputs.

## Panel D: Runtime and portability (about 15% of figure)

Report wall time, peak memory and tested profile for the frozen release across
the compact test bundle and the publication benchmark dataset.

Current status: the generated panel uses the repository's stub resource profile
only. Replace or supplement it with production measurements before submission.
