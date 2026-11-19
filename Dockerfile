# TRIARCHY production task image.
#
# Contains every tool the Nextflow processes call (fastp, Bowtie2, Kraken2,
# Bracken, HUMAnN, MetaPhlAn, DIAMOND, cutadapt, ITSxpress, VSEARCH, R/DADA2)
# installed from the pinned linux-64 lock file. Reference databases are NOT
# included; they are mounted from the host at run time by the container
# profiles in conf/docker.config and conf/apptainer.config.
#
# Build:  docker build -t ghcr.io/aglucaci/triarchy:$(cat VERSION) .
FROM mambaorg/micromamba:1.5.10-bookworm-slim

ARG TRIARCHY_VERSION=dev
LABEL org.opencontainers.image.title="TRIARCHY" \
      org.opencontainers.image.description="Control-aware integration of WMS, 16S, and ITS microbiome assays" \
      org.opencontainers.image.source="https://github.com/aglucaci/TRIARCHY" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.version="${TRIARCHY_VERSION}"

# Nextflow needs `ps` inside task containers to collect resource metrics.
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER locks/environment.linux-64.lock.txt /tmp/environment.lock.txt
RUN micromamba install --yes --name base --file /tmp/environment.lock.txt \
    && micromamba clean --all --yes \
    && rm /tmp/environment.lock.txt

COPY --chown=$MAMBA_USER:$MAMBA_USER VERSION triarchy /opt/triarchy/
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ /opt/triarchy/scripts/
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/ /opt/triarchy/workflow/

ENV PATH="/opt/conda/bin:/opt/triarchy/workflow/bin:/opt/triarchy:${PATH}" \
    TRIARCHY_HOME=/opt/triarchy
