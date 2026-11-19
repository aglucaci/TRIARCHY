# Workflow architecture

The production graph is:

```text
metadata validation
  → assay-specific QC
  → WGS host removal, taxonomy, and functions
  → 16S trimming and DADA2
  → ITS trimming, ITSxpress, and DADA2
  → linked-control annotation
  → taxonomy harmonization
  → authenticity and discordance classification
  → assay-native diversity and cross-assay concordance
  → dashboard and provenance
```

The implementation lives in `main.nf`, `modules/local/`, and
`workflow/scripts/`.
