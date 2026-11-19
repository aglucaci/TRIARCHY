# Contributing to TRIARCHY

Thank you for helping improve TRIARCHY. This guide explains how to report
problems, propose changes, and get a pull request merged.

By participating you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Ways to contribute

- **Report a bug** with the [bug report form](https://github.com/aglucaci/TRIARCHY/issues/new?template=bug_report.yml).
- **Request a feature** or a new assay, database, or executor profile with the
  [feature request form](https://github.com/aglucaci/TRIARCHY/issues/new?template=feature_request.yml).
- **Ask a usage question** in [GitHub Discussions](https://github.com/aglucaci/TRIARCHY/discussions)
  or as an issue if Discussions is not enabled.
- **Report a security problem** privately as described in [SECURITY.md](SECURITY.md).
  Do not open a public issue for it.
- **Improve documentation**, examples, or error messages.
- **Fix bugs or add features** through a pull request.

For anything larger than a small fix, please open an issue first so the design
can be agreed before you invest time in code.

## Development setup

```bash
git clone https://github.com/aglucaci/TRIARCHY.git
cd TRIARCHY
mamba env create --file environment.yml
mamba activate TRIARCHY
```

Most changes to the CLI, metadata validation, and evidence model can be
developed without the full production reference set. Changes to Nextflow
processes that call Kraken2, HUMAnN, MetaPhlAn, DADA2, or ITSxpress need the
references installed with `./install.sh` (see the README).

## Repository layout

| Path | Contents |
| --- | --- |
| `triarchy` | Python command-line entry point (`run`, `demo`, `real-world`, `validate`, `databases`) |
| `main.nf`, `nextflow.config` | Nextflow entry workflow and global configuration |
| `modules/local/` | Nextflow subworkflows: `metadata`, `wms`, `amplicon`, `integration` |
| `workflow/scripts/` | Evidence model, integration, and DADA2 scripts run inside tasks |
| `workflow/bin/` | Executables placed on the task `PATH` |
| `scripts/` | Installer, dataset download, validation, and demo build/evaluation helpers |
| `conf/` | Executor profiles, primers, assay performance, and evidence policy |
| `data/` | Benchmark and demo inputs (generated data is not committed) |
| `results/production_demo/` | Reference outputs of the production demo |

## Making a change

1. Fork the repository and create a branch from `main`
   (for example `fix/kraken-memory-mapping` or `feat/nanopore-16s`).
2. Keep each pull request focused on one change.
3. Follow the style guidelines below.
4. Run the checks under [Testing](#testing).
5. Add an entry under `## [Unreleased]` in [CHANGELOG.md](CHANGELOG.md) for any
   user-visible change.
6. Update the README and any affected configuration documentation.
7. Open a pull request and complete the template.

### Commit messages

Write commit messages in the imperative mood with a short summary line
(about 72 characters or fewer), for example `Add Apptainer cache directory option`.
Explain *why* in the body when it is not obvious.

## Style guidelines

- **Python**: target Python 3.12, standard library only unless a dependency is
  added to `environment.yml`. Follow PEP 8, use type hints, and keep
  `from __future__ import annotations` at the top of modules.
- **Nextflow**: DSL2. Keep process resource labels (`low`, `medium`, `high`)
  consistent with `nextflow.config`. Every process should emit versions and
  provenance consistent with the existing modules.
- **Shell**: Bash with `set -euo pipefail`. Scripts must pass `bash -n` and
  should pass [ShellCheck](https://www.shellcheck.net/).
- **Configuration**: primers belong in `conf/primers.config` and evidence
  thresholds in `conf/evidence_policy.json`. Do not hard-code them in processes.
- **Line endings**: LF only (enforced by `.gitattributes` and `.editorconfig`).

## Scientific changes

Changes that alter taxonomic calls, evidence classification, control handling,
or default thresholds affect published results. For these:

- Explain the scientific rationale in the pull request and cite sources.
- Show the effect on the production demo (`triarchy demo`) by comparing
  `results/production_demo/production_demo_evaluation.json` before and after.
- Update the pinned database versions and checksums together, never separately.

## Testing

Run the fast checks before every pull request. Continuous integration runs the
same checks:

```bash
python -m py_compile triarchy scripts/*.py workflow/scripts/*.py
bash -n install.sh setup_and_run_triarchy.sh scripts/*.sh
./triarchy --help
```

Validate study metadata and the demo inputs without running Nextflow:

```bash
triarchy demo --build-only
```

When you change workflow logic and have the references installed, run the full
demo and confirm the evaluation still passes:

```bash
triarchy demo --profile laptop --resume
```

## Review and merging

A maintainer reviews every pull request. Expect questions about scientific
rationale and reproducibility as well as code. Pull requests are merged into
`main` once CI passes and review comments are resolved.

## Releases

Maintainers cut releases by:

1. Moving `## [Unreleased]` entries in `CHANGELOG.md` under a new version heading.
2. Updating `version` and `date-released` in `CITATION.cff`.
3. Tagging the commit as `vX.Y.Z` and publishing a GitHub release.

TRIARCHY follows [Semantic Versioning](https://semver.org/). Changes that alter
outputs for the same inputs and references are treated as at least a minor
version change.

## License

By contributing, you agree that your contributions are licensed under the
[GNU General Public License v3.0](LICENSE), the same license as the project.
