# Installation

## Requirements

- Linux, macOS, or WSL Ubuntu
- Bash and Java 17 or later
- Mamba, Micromamba, or Conda
- Nextflow (the pinned version is installed locally)
- Production reference storage; the complete reference set can exceed 500 GB
- Real WGS, 16S, and ITS FASTQs with study-specific primer information

Apple-silicon macOS production installs use native tool shims where required.
For large studies, use a workstation, HPC scheduler, or cloud executor rather
than a memory-constrained laptop.

## Minimum production specifications

For a complete production run using the full Kraken2 PlusPFP, MetaPhlAn,
HUMAnN, host, SILVA, and ITS reference set, provision at least:

| Resource | Minimum production specification |
| --- | --- |
| Largest compute-job slot | 12 CPU cores |
| Memory per high-resource job | 96 GB RAM |
| Maximum job wall time | 48 hours |
| Free space for reference installation | 600 GB |
| Recommended shared storage | 1 TB or more |
| GPU | Not required |

The built-in SLURM profile requests the following resources per task:

| Task class | CPUs | Memory | Wall time |
| --- | ---: | ---: | ---: |
| Low | 1 | 2 GB | 1 hour |
| Medium | 4 | 16 GB | 8 hours |
| High | 12 | 96 GB | 48 hours |

For HPC execution, use Linux (preferably x86-64), SLURM with `sbatch`, Bash,
Java 17 or later, and the named `TRIARCHY` Mamba/Conda environment. The
repository, reference databases, input FASTQs, result directory, and Nextflow
work directory must be on shared POSIX storage visible at identical paths from
the login and compute nodes. Internet access is needed during database
installation unless the validated archives are supplied locally.

The database installer enforces 600 GB of free space for a new reference
installation. Production work files and results require additional capacity;
1 TB of shared storage is a practical starting point, with more space needed
for large input studies. Nodes with 128 GB RAM are recommended to leave
operating-system and database-loading headroom above the 96 GB process request.
The `laptop` profile's 8 GB high-task cap is intended for the constrained demo
and is not a dependable baseline for full production datasets.

## Install

Create the named environment and install the complete production references:

```bash
mamba env create --file environment.yml
mamba activate TRIARCHY

./install.sh \
  --database-dir /path/to/production/references \
  --cpus 8 \
  --yes
```

The installer downloads resumably, validates completed stages, writes
`conf/references.config`, and records the selected reference root in
`references/.database-dir`.

Useful installer options:

| Option | Meaning |
| --- | --- |
| `--database-dir DIR` | Production reference destination |
| `--silva-archive FILE` | Install a pre-existing, checksum-verified SILVA archive |
| `--kraken-archive FILE` | Reuse an existing full PlusPFP archive |
| `--cpus N` | Installation/indexing threads |
| `--locked-dependencies` | Install the pinned dependency lock |
| `--skip-nextflow` | Keep the existing Nextflow installation |
| `--skip-dependencies` | Keep the existing `TRIARCHY` environment |
| `--skip-databases` | Keep the current production references |
| `--yes` | Accept the large non-interactive installation |

SILVA downloads are automatic and mirror-aware. TRIARCHY tries the official
SILVA release server first and then a checksum-identical Zenodo mirror. Each
source has its own resumable partial download, unreachable sources fail over
automatically, and the completed archive must match MD5
`41f1d231ce5079b88f4c01742c288ba9`. Existing installations of the earlier
full-length DADA2 artifact remain accepted by their pinned checksum. The
optional `--silva-archive` input is only needed for fully air-gapped systems.

TRIARCHY uses Kraken2/Bracken for WMS taxonomy and installs the complete HUMAnN
ChocoPhlAn and UniRef90 databases for functional profiling. It also installs
the checksum-pinned `mpa_vJun23_CHOCOPhlAnSGB_202307` MetaPhlAn marker database
required by HUMAnN 3.9. HUMAnN prescreens against this shared index before
constructing a small sample-specific ChocoPhlAn search database. The database
path and release are passed explicitly with offline lookup enabled, so compute
nodes do not need internet access and do not build the full ChocoPhlAn index.

The full Kraken2 PlusPFP database is larger than the 96 GB production job
request. TRIARCHY therefore enables Kraken2's supported `--memory-mapping` mode
and limits WMS taxonomy to two concurrent jobs on SLURM. Place the database on
fast parallel or node-local storage when possible; memory-mapped classification
trades lower process memory for additional storage I/O.

## Reproducible (locked) installs

`environment.yml` lists TRIARCHY's dependencies with loose version
constraints, so a fresh solve can pick up newer releases. For reproducible
installs, use the committed explicit lock for your platform:

```bash
./install.sh --locked-dependencies --database-dir /path/to/production/references --yes
```

| Platform | Lock file |
| --- | --- |
| Linux x86-64 | `locks/environment.linux-64.lock.txt` |
| macOS Apple silicon | `locks/environment.osx-arm64.lock.txt` |

On Apple silicon, Bioconda does not provide current builds of every tool (for
example Bracken); the installer adds native Homebrew shims for those, as in an
unlocked install.

Maintainers regenerate the locks after changing `environment.yml`:

```bash
conda-lock lock -f environment.yml -p linux-64 -p osx-arm64 \
  --kind explicit --filename-template 'locks/environment.{platform}.lock.txt'
```

Each install also records the exact packages it installed in
`environment.lock.txt` (machine-specific and not committed).

## Containers

The `docker`, `apptainer`, and `singularity` profiles run every task in the
versioned image `ghcr.io/aglucaci/triarchy:<version>`, built from the
repository `Dockerfile` with the linux-64 lock. The image contains the tools
only; reference databases stay on the host and are mounted read-only. Nextflow
itself and the `triarchy` CLI still run on the host, so you need Java and
Nextflow there, but not the Conda environment's bioinformatics tools.

The image is published by the **Container image** GitHub Actions workflow when
a `vX.Y.Z` tag is pushed. To build it locally:

```bash
docker build --build-arg TRIARCHY_VERSION="$(cat VERSION)" \
  -t "ghcr.io/aglucaci/triarchy:$(cat VERSION)" .
```

For publication runs, pin the image by digest
(`ghcr.io/aglucaci/triarchy@sha256:...`) in a site config passed with
`--site-config`.
