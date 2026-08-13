# Release readiness for a Bioinformatics Application Note

This checklist separates work that is already executable in the repository from
release artifacts that must exist at submission time.

## Implemented contribution

- Stage-linked sample-control metadata validation.
- Assay-native WMS, 16S, and ITS processing branches.
- Control annotation using linked-control prevalence, paired-control abundance,
  and batch localization.
- Taxonomic harmonization with source identifiers, canonical identifiers,
  mapping method, mapping confidence, and comparison rank.
- Two-stage authenticity/explanation inference in `taxon_evidence.tsv`, with
  `biological_authenticity`, uncertainty bounds, parent `cross_assay_status`,
  detailed `classification`, and legacy score columns retained for older
  consumers.
- Stage-specific contamination evidence and probable contamination-source
  reporting across field, extraction, PCR, library, run, and other controls.
- Taxonomic compatibility scoring and status reporting for rank-aware
  cross-assay interpretation.
- Versioned evidence policy in `conf/evidence_policy.json`.
- Assay performance parameters in `conf/assay_performance.tsv`.
- Truth-based comparison utility in `scripts/compare_evidence_methods.py`.
- Resource profiling utility in `scripts/benchmark_resources.py`.

## Submission gates

- Public source repository with issue tracker enabled.
- Tagged `v1.0.0` release.
- Zenodo or Software Heritage archive for the submitted source version.
- Public benchmark bundle with checksums and a persistent identifier.
- Container image tags and immutable digests for every production profile.
- Documentation site or rendered documentation with stable URLs.
- Example HTML report generated from the test profile.
- At least one empirical matched WMS/16S/ITS dataset or mock-community run.
- Baseline comparison table for independent assays, naive merging, shared
  taxonomy, per-assay contamination filtering, global controls, ablated
  TRIARCHY, and full TRIARCHY.
- Cover letter explaining novelty, maintenance, availability, and any use of
  AI-assisted writing or code generation required by the journal.

## Maintainer commitment

For the submitted release, maintainers should commit to:

- keeping the `v1.0.0` source archive, containers, benchmark data, and
  documentation available for at least two years after publication;
- preserving submitted artifacts even if later releases change interfaces;
- assigning at least two people release permissions;
- responding to security or correctness reports through the public issue
  tracker;
- documenting supported Nextflow, container, operating-system, and reference
  database versions.

## Recommended one-command checks

```bash
bin/triarchy validate \
  --input benchmarks/TRI-01/metadata/samples.tsv \
  --paired-controls benchmarks/TRI-01/metadata/control_relationships.tsv
```

```bash
bin/triarchy databases check \
  --manifest conf/databases.example.yml \
  --destination references/data \
  --report output/database-validation.json
```

```bash
./run_triarchy.sh --stub \
  --samples benchmarks/TRI-01/metadata/samples.tsv \
  --controls benchmarks/TRI-01/metadata/control_relationships.tsv \
  --outdir results/TRI-01-stub
```
