## Summary

<!-- What does this change and why? Link the issue it addresses, e.g. "Closes #12". -->

## Type of change

- [ ] Bug fix
- [ ] New feature
- [ ] Scientific change (alters calls, thresholds, evidence rules, or reference versions)
- [ ] Documentation
- [ ] Maintenance / CI

## Testing

<!-- Describe what you ran. For scientific changes, compare the production demo evaluation before and after. -->

- [ ] `pytest`
- [ ] `ruff check .`
- [ ] `shellcheck -S warning -s bash scripts/*.sh hpc/*.sbatch`
- [ ] `triarchy demo --build-only`
- [ ] `triarchy demo --profile <profile> --resume` (for workflow changes)

## Checklist

- [ ] I have read [CONTRIBUTING.md](../CONTRIBUTING.md)
- [ ] I updated `docs/`, `nextflow_schema.json`, and configuration where needed
- [ ] I added or updated tests in `tests/`
- [ ] I added an entry under `## [Unreleased]` in `CHANGELOG.md`
- [ ] No protected, identifiable, or large data files are included
