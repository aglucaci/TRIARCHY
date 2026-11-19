# TRIARCHY

![TRIARCHY logo](assets/logo/triarchy_logo_transparent.png)

**TRIARCHY is a production workflow for control-aware integration of matched
whole-genome shotgun (WGS) metagenomics, 16S rRNA gene sequencing, and fungal
ITS sequencing in environmental microbiome studies.**

TRIARCHY processes each assay independently, links every production library to
its experimentally relevant controls, harmonizes taxonomy with provenance, and
classifies where assays agree, disagree, or lack enough power for a defensible
conclusion. This repository supports production execution only.


## Quick start

```bash
git clone https://github.com/aglucaci/TRIARCHY.git
cd TRIARCHY
# Software: the TRIARCHY conda environment, Nextflow, and the triarchy command
scripts/install_environment.sh --yes
mamba activate TRIARCHY

# Reference databases (> 500 GB; rerun the same command to resume)
triarchy databases install --destination /path/to/production/references --cpus 8 --yes

# Verify the installation end to end on the deterministic demo
triarchy demo --profile laptop --resume

# Run a study
triarchy run \
  --profile local \
  --samples /data/study/samples.tsv \
  --controls /data/study/control_relationships.tsv \
  --outdir results/study
```

Production runs need Linux, macOS, or WSL with Java 17+, and a workstation, HPC
cluster, or cloud executor sized for the full reference set (12 CPUs and 96 GB
RAM per high-resource task). See [installation](docs/installation.md) for the
full requirements.

## Documentation

| Topic | Contents |
| --- | --- |
| [Installation](docs/installation.md) | Requirements, sizing, installer options, locked installs, containers |
| [Inputs](docs/inputs.md) | Sample manifest and control relationship formats |
| [Running TRIARCHY](docs/usage.md) | `triarchy run`, resuming, execution profiles, containers |
| [Parameters](docs/parameters.md) | Every pipeline parameter (generated from `nextflow_schema.json`) |
| [Outputs](docs/outputs.md) | Result files and the evidence dashboard |
| [Production demo](docs/demo.md) | Deterministic synthetic benchmark and its truth table |
| [HPC guide](hpc/README.md) | Step-by-step SLURM setup and job scripts for `TRI-REAL-05` |
| [Public datasets](docs/datasets.md) | `TRI-REAL-05` empirical cohort with linked blanks and ZymoBIOMICS truth; HPC instructions |
| [Architecture](docs/architecture.md) | Workflow graph and source layout |
| [Troubleshooting](docs/troubleshooting.md) | Debugging failed tasks |

## Contributing

Bug reports, feature requests, and pull requests are welcome. See
[CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, checks, and
review process, and the [Code of Conduct](CODE_OF_CONDUCT.md). Report security
issues privately as described in [SECURITY.md](SECURITY.md). Notable changes
are recorded in [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use TRIARCHY, please cite it using the metadata in
[CITATION.cff](CITATION.cff) (GitHub's "Cite this repository" button reads it).

## License

TRIARCHY is released under the GNU General Public License v3.0. See
[LICENSE](LICENSE).
