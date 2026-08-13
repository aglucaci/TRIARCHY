# TRIARCHY governance and maintenance

TRIARCHY is maintained as a research software project for environmental
microbiome evidence integration.

## Release responsibility

Each public release must identify:

- the release manager;
- at least one additional maintainer with permission to publish emergency
  fixes;
- the Git tag, source archive, container digests, benchmark-data archive, and
  documentation version;
- the supported Nextflow and Python versions;
- the supported execution profiles.

The submitted publication release must remain accessible for at least two years
after publication, even if later releases change the workflow interface.

## Issue handling

Use the public issue tracker for reproducible bug reports, documentation
problems, benchmark-data problems, and feature requests. Security-sensitive or
data-sensitive reports should be routed to the contact address listed in the
publication and release notes.

## Scientific claims

TRIARCHY evidence classes summarize assay support, linked controls, taxonomic
resolution, and detection power. They are not declarations of biological truth.
Any manuscript, release note, or documentation page that reports accuracy must
state the benchmark or empirical data on which the claim is based.

## Archival policy

For each publication-associated release, archive:

- source code;
- test data;
- benchmark truth tables;
- reference manifests and checksums;
- container image digests;
- generated example report;
- scripts required to reproduce tables and figures.
