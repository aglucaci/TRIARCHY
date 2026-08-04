# TRIARCHY workflow

TRIARCHY preserves assay-native measurements while integrating taxonomic,
control, and ecological evidence.

The canonical standalone Mermaid source is
[`docs/overall-workflow-dag.mmd`](docs/overall-workflow-dag.mmd). Generate the
publication-ready vector PDF with:

```bash
.tools/env/bin/python scripts/render_workflow_dag_pdf.py
```

The default output is `output/pdf/triarchy-overall-workflow-dag.pdf`. Override
it with `--output PATH`. The renderer uses `reportlab`, which is declared in
`environment.yml` and installed by `install.sh`.

```mermaid
flowchart TD
    META["Manifests + metadata + stage-specific controls"] --> GRAPH["Validated sample–control graph"]

    GRAPH --> WMS["WMS reads"]
    GRAPH --> S16["16S reads"]
    GRAPH --> ITS["ITS reads"]

    WMS --> WQC["fastp QC"]
    WQC --> HOST["Bowtie2 host removal"]
    HOST --> WTAX["Kraken2/Bracken taxonomy"]
    HOST --> FUNC["HUMAnN functional potential"]

    S16 --> SQC["cutadapt QC"]
    SQC --> SDADA["DADA2 ASVs"]
    SDADA --> STAX["Bacterial/archaeal taxonomy"]

    ITS --> IQC["cutadapt QC"]
    IQC --> IX["ITSxpress"]
    IX --> IDADA["DADA2 ASVs"]
    IDADA --> ITAX["UNITE fungal taxonomy"]

    GRAPH --> CONTROL["Linked-control model"]
    WTAX --> CONTROL
    STAX --> CONTROL
    ITAX --> CONTROL

    WQC --> QMET["Measured library QC + depth"]
    HOST --> QMET
    SQC --> QMET
    IQC --> QMET

    CONTROL --> RAW["Unfiltered profiles"]
    CONTROL --> ANNOTATED["Control-annotated profiles"]
    CONTROL --> HIGH["High-confidence profiles"]

    WTAX --> HARM["Taxonomic harmonization + provenance"]
    STAX --> HARM
    ITAX --> HARM
    ANNOTATED --> EVIDENCE["Evidence-v3 inference"]
    HARM --> EVIDENCE
    QMET --> EVIDENCE
    PERF["Assay performance + reference coverage"] --> EVIDENCE
    POLICY["Versioned thresholds + calibration scope"] --> EVIDENCE
    EVIDENCE --> MATRIX["Taxon evidence matrix + biological state"]
    EVIDENCE --> MODEL["Serialized evidence_model.json"]

    WMSDIV["WMS taxonomic diversity"]
    S16DIV["16S diversity"]
    ITSDIV["ITS diversity"]
    FUNCDIV["WMS functional diversity"]
    HIGH --> WMSDIV
    HIGH --> S16DIV
    HIGH --> ITSDIV
    FUNC --> FUNCDIV

    WMSDIV --> ECO["Cross-assay ecological concordance"]
    S16DIV --> ECO
    ITSDIV --> ECO
    EVIDENCE --> ECO

    MATRIX --> REPORT["TRIARCHY analytical report"]
    MODEL --> REPORT
    ECO --> REPORT
    FUNC --> REPORT
    GRAPH --> REPORT
```

The integration layer compares evidence only where assay scope overlaps:

- WMS and 16S for bacteria and archaea;
- WMS and ITS for fungi;
- WMS alone for functional potential.

Joint diversity is not calculated unless a separately validated multi-table or
latent-community model supplies a justified common representation.
