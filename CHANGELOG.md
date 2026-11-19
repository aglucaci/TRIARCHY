# Changelog

All notable changes to TRIARCHY are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `TRI-REAL-05`, the new default public cohort for `triarchy real-world`:
  whole extraction plates from the Shaffer et al. (2022) extraction study
  (Qiita 12201), built from the public Qiita sample information. Every
  extraction blank is linked to the samples on its plate, and the
  ZymoBIOMICS D6310 dilution series is scored against a truth table
  transcribed from the Zymo instruction manual. Options select extraction
  kits, plates and a per-plate sample cap.
- `conf/cohorts/shaffer_emp.config` with the study's EMP primers and ITS1
  extraction, applied automatically to that cohort.
- `scripts/evaluate_mock_truth.py`, run automatically after a real-world run,
  reporting per-assay recall by abundance and dilution and the evidence model's
  calls on expected and unexpected taxa.

## [1.0.0] - 2026-11-19

### Added

- `VERSION` file, `manifest.version` in `nextflow.config`, and `triarchy --version`.
  `pipeline_info/provenance.json` now records `triarchy_version`, shown in the report.
- `Dockerfile` for the task image and a GitHub Actions workflow that publishes
  `ghcr.io/aglucaci/triarchy:<version>` on version tags.
- Committed per-platform Conda lock files in `locks/` (linux-64, osx-arm64);
  `./install.sh --locked-dependencies` installs the lock for the current platform.
- `nextflow_schema.json` describing every parameter, and generated
  `docs/parameters.md`.
- pytest suite in `tests/` covering metadata validation, the evidence model, and
  version/schema consistency.
- CI: ruff, ShellCheck, pytest, `nextflow config` for every profile, and
  `nextflow lint`.
- User documentation in `docs/`; the README is now an overview.
- Contributor and maintenance documentation: `CONTRIBUTING.md`,
  `CODE_OF_CONDUCT.md`, `SECURITY.md`, `CITATION.cff`, and this changelog.
- GitHub issue forms, pull request template, CODEOWNERS, and Dependabot.
- `.editorconfig` for consistent formatting across editors.

### Changed

- The `docker`, `apptainer`, and `singularity` profiles now load
  `conf/*.config`, use the published image, mount the project and reference
  directories read-only, and no longer force the local executor (so they
  combine with `slurm`). Previously they used no image.
- The `manuscript/` directory is no longer tracked in the repository.

### Fixed

- Removed an unused import and variable reported by ruff, and ShellCheck
  warnings in the installer scripts.

[Unreleased]: https://github.com/aglucaci/TRIARCHY/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/aglucaci/TRIARCHY/releases/tag/v1.0.0
