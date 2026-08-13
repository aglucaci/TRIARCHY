# Examining TRIARCHY benchmark results

Each completed benchmark is stored under `results/benchmarks/<DATASET_ID>`.
Start with the two HTML reports:

- `triarchy_report.html` summarizes quality, assay completeness, taxonomic
  evidence, contamination, diversity, and cross-assay concordance.
- `pipeline_info/execution_report.html` summarizes Nextflow task execution and
  resource use.

The workflow graph and execution timing are available at:

- `pipeline_info/dag.html` - interactive workflow DAG;
- `pipeline_info/timeline.html` - process timing;
- `pipeline_info/trace.tsv` - task-level status and resource data;
- `pipeline_info/provenance.json` - run provenance.

For example, to open the TRI-12 low-biomass contamination benchmark on macOS:

```bash
open results/benchmarks/TRI-12/triarchy_report.html
open results/benchmarks/TRI-12/pipeline_info/dag.html
open results/benchmarks/TRI-12/pipeline_info/execution_report.html
```

## Important tables

| File | Question it answers |
| --- | --- |
| `qc/library_qc.tsv` | Which libraries passed, failed, or had limited evidence? |
| `profiles/sample_confidence.tsv` | Which biological samples are suitable for analysis? |
| `integrated/taxon_evidence.tsv` | What is the authenticity estimate and cross-assay explanation for each sample-taxon call? |
| `integrated/cross_assay_concordance.tsv` | How strongly do the appropriate assays agree? |
| `integrated/taxonomy_crosswalk.tsv` | How were source classifications harmonized? |
| `diversity/alpha_diversity.tsv` | What is within-sample diversity by assay? |
| `diversity/beta_diversity.tsv` | What are between-sample distances by assay? |

Inspect a table in the terminal with:

```bash
column -t -s $'\t' results/benchmarks/TRI-12/integrated/taxon_evidence.tsv | less -S
```

In `integrated/taxon_evidence.tsv`, start with these columns:

- `biological_authenticity`: Stage-I estimate that the taxon is biologically authentic in the sample.
- `authenticity_ci_low`, `authenticity_ci_high`, `authenticity_uncertainty`: approximate uncertainty interval.
- `authenticity_state`: conservative support state; wide intervals remain `indeterminate`.
- `cross_assay_status`: parent explanation such as `concordant`, `assay_limited`, `taxonomy_limited`, `contamination_associated`, `genuine_discordance`, or `indeterminate`.
- `classification`: detailed audit/status label such as `rank_limited_agreement` or `assay_limited_nondetection`.
- `probable_contamination_source`: strongest linked-control stage when contamination evidence is present.
- `taxonomic_compatibility_status`: whether compatible calls are exact/canonical, hierarchical, single-assay, ambiguous, or incompatible.

List task failures or non-completed tasks with:

```bash
awk -F '\t' 'NR == 1 || $5 != "COMPLETED"' \
  results/benchmarks/TRI-12/pipeline_info/trace.tsv | column -t -s $'\t'
```

Compare an observed table with the corresponding truth table. For example:

```bash
column -t -s $'\t' benchmarks/TRI-12/truth/expected_taxon_status.tsv | less -S
column -t -s $'\t' results/benchmarks/TRI-12/integrated/taxon_evidence.tsv | less -S
```

## Creating a DAG on a new run

The `dag` block in `nextflow.config` enables DAG generation automatically. Any
run through `run_triarchy.sh` writes the graph to
`<OUTDIR>/pipeline_info/dag.html`:

```bash
./run_triarchy.sh --stub \
  --samples benchmarks/TRI-12/metadata/samples.tsv \
  --controls benchmarks/TRI-12/metadata/control_relationships.tsv \
  --outdir results/benchmarks/TRI-12
```

Nextflow regenerates the DAG at the end of the run. Use `-resume` to reuse
cached processes while rebuilding reports:

```bash
./run_triarchy.sh --stub -resume \
  --samples benchmarks/TRI-12/metadata/samples.tsv \
  --controls benchmarks/TRI-12/metadata/control_relationships.tsv \
  --outdir results/benchmarks/TRI-12
```
