# Contributing to TRIARCHY

Contributions should preserve the central design rule: WMS, 16S, and ITS are
not interchangeable abundance scales. New features should either improve an
assay-native branch or make the cross-assay evidence layer more transparent.

## Development checks

Run focused tests before opening a pull request:

```bash
python -m unittest discover tests
```

For workflow changes, also run a deterministic stub dataset:

```bash
./run_triarchy.sh --stub \
  --samples benchmarks/TRI-01/metadata/samples.tsv \
  --controls benchmarks/TRI-01/metadata/control_relationships.tsv \
  --outdir results/TRI-01-stub
```

## Pull request expectations

- Describe the behavior change and the reason for it.
- Add or update tests for metadata validation, harmonization, evidence
  classification, reporting, or benchmark comparisons as appropriate.
- Do not replace assay-native values with cross-assay averaged abundances.
- Record new output columns in `docs/data-contracts.md`.
- Update `docs/methodology.md` when evidence logic changes.

## Data and references

Do not commit large production databases. Commit manifests, checksums, small
fixtures, and reproducible acquisition instructions instead.
