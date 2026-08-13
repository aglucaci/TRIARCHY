# TRIARCHY references

Run `./install.sh --yes` to create `references/data/`, download the configured
reference datasets, build indexes, and generate:

- `references/manifest.tsv` with versions, sources, timestamps, and file
  checksums;
- `conf/references.config` with absolute Nextflow parameter paths.

These generated files and the databases are intentionally ignored by Git.

The standard database profile installs:

- NCBI GRCh38.p14 (`GCF_000001405.40`) and a Bowtie2 host index;
- the June 2026 Kraken2 PlusPFP database and matching Bracken distributions;
- full HUMAnN ChocoPhlAn and UniRef90 Diamond databases;
- SILVA 138.2 formatted for DADA2;
- EUKARYOME 2.0 99%-clustered ITS taxonomy formatted for DADA2.

The ITS default is EUKARYOME because it provides a stable direct, versioned
DADA2 archive. To install a different DADA2-compatible ITS reference, set
`TRIARCHY_ITS_URL` before running the installer.

Database licenses and citation requirements remain those of their upstream
providers. Review them before redistributing any downloaded content.
