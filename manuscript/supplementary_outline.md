# TRIARCHY supplementary material outline

This file is the source checklist for the single Supplementary Material PDF
that should accompany a Bioinformatics Application Note submission.

## Supplementary Methods

S1. Workflow implementation and command-line entry points

S2. Input schemas

S3. Sample-control graph

S4. Taxonomic harmonization

S5. Authenticity and cross-assay explanation algorithm

- Stage-I `biological_authenticity` definition, uncertainty bounds and calibration status.
- Stage-II mapping from six parent `cross_assay_status` explanations to detailed software statuses.
- Stage-specific contamination scores and taxonomic compatibility states.

S6. Benchmark-data generation

S7. Baseline comparison and statistical evaluation

S8. Software, container, database-location and database-version records

S9. Empirical mock-community or environmental validation dataset

## Supplementary Tables

S1. Supported inputs

S2. Parameters and defaults

S3. Authenticity and explanation decision rules

S4. Benchmark compositions

S5. Baseline methods

S6. Performance metrics

S7. Runtime and resource use

S8. Tested computing environments

S9. Generated artifacts excluded from version control

## Supplementary Figures

S1. Complete Nextflow DAG: `figures/supplementary_figure_s1_dag.svg`

S2. Control-link examples: `figures/supplementary_figure_s2_control_graph.svg`

S3. Taxonomic-crosswalk validation: `figures/supplementary_figure_s3_taxonomy_crosswalk.svg`

S4. Detection limit curves: `figures/supplementary_figure_s4_detection_power.svg`

S5. Contaminant performance/design: `figures/supplementary_figure_s5_contamination_design.svg`

S6. Environmental case-study design: `figures/supplementary_figure_s6_environmental_design.svg`

S7. Installation and execution tests: `figures/supplementary_figure_s7_installation_tests.svg`

## Required frozen-release values

- Public source URL.
- Documentation URL.
- Version tag.
- Source archive DOI.
- Benchmark-data DOI.
- Container image digests.
- Reference-manifest checksums.
- External database destination and `references/.database-dir` record.
- `setup_and_run_triarchy.sh` install-and-run command examples.
- Generated-output cleanup policy for `work/`, `results/`, `.nextflow/`,
  `output/`, `tmp/`, Python bytecode and manuscript build products.
- Contact email.
- Baseline comparison results.
- Empirical mock-community or environmental validation results.
- Runtime, memory and tested-platform table, including Linux/Docker,
  macOS/Docker, WSL2 and Slurm/Apptainer status when available.
