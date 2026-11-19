# Production outputs

The output directory contains:

```text
qc/library_qc.tsv
profiles/unfiltered/taxa.tsv
profiles/unfiltered/functions.tsv
profiles/control_annotated/taxa.tsv
profiles/control_annotated/functions.tsv
profiles/high_confidence/taxa.tsv
profiles/sample_confidence.tsv
integrated/harmonized_taxa.tsv
integrated/taxonomy_crosswalk.tsv
integrated/taxon_evidence.tsv
integrated/cross_assay_concordance.tsv
integrated/evidence_model.json
diversity/alpha_diversity.tsv
diversity/beta_diversity.tsv
pipeline_info/provenance.json
pipeline_info/execution_report.html
pipeline_info/timeline.html
triarchy_report.html
```

Open `triarchy_report.html` for the scientific evidence dashboard. It follows
the interpretation sequence:

1. Did the experiment work?
2. Which libraries and samples are trustworthy?
3. Which taxa are supported?
4. Why do assays disagree?
5. Which ecological conclusions survive?
6. Can the evidence and provenance be exported?

Cross-assay classifications summarize evidence. They do not make WGS, 16S,
and ITS relative abundances quantitatively interchangeable.
