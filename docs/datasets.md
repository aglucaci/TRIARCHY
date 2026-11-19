# Public empirical cohort

`TRI-REAL-01` is a small empirical validation subset of the UC San Diego DNA
extraction comparison study [ENA ERP124610](https://www.ebi.ac.uk/ena/browser/view/ERP124610),
also registered as [BioProject PRJEB40903](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB40903)
and [Qiita study 12201](https://qiita.ucsd.edu/study/description/12201).
The source experiment compared extraction methods for microbial-community DNA
using 16S, ITS, and shotgun metagenomic sequencing [1].

The default TRIARCHY selection is deliberately small enough for an integration
run while retaining matched assay structure:

| Property | Default `TRI-REAL-01` selection |
| --- | --- |
| Biological material | Dust microbial-community samples from the extraction comparison |
| Selected groups | Three biological groups and one blank group |
| Libraries | 11 ENA runs represented by 15 FASTQ files |
| Assays | Three matched biological WMS/16S/ITS triplets; blank 16S and WMS libraries |
| Layout | Single-end amplicon reads and paired-end WMS reads |
| Download size | Approximately 2.96 MB |
| Integrity | ENA-provided MD5 verification plus gzip and manifest validation |
| Local provenance | `metadata/ena_runs.tsv`, `selected_runs.tsv`, `download_manifest.tsv`, and `dataset.json` |

This is a selected subset, not a redistribution or complete mirror of the
7,296-run source study. The downloader queries ENA at execution time, records
the returned accession metadata and URLs, and retrieves the chosen FASTQs
directly from ENA. `--max-groups` changes the subset and may materially change
its size and assay balance.

Download and run it with:

```bash
mamba activate TRIARCHY
triarchy real-world \
  --dataset TRI-REAL-01 \
  --profile slurm \
  --site-config /shared/project/triarchy-site.config \
  --work-dir /shared/scratch/$USER/triarchy-real-01-work \
  --outdir /shared/project/results/TRI-REAL-01
```

Run this command from a login or transfer node with internet access. FASTQs are
downloaded into `data/real_world/TRI-REAL-01`, verified against ENA MD5 values,
and then submitted through the normal production workflow. To stage the data
without launching Nextflow, add `--download-only`; later use `--skip-download`
to run the staged copy on an offline HPC.

The public metadata does not currently provide curated pairings between the
selected blank and biological libraries. The generated control relationship
manifest is therefore empty, and TRIARCHY reports this limitation rather than
inventing control links. Curate that manifest before interpreting contamination
classifications.

## Other catalogued public studies

The metadata downloader also catalogues larger studies for future curation.
They are not exposed through `triarchy real-world` because their ENA metadata
does not currently yield a validated, automatically paired WMS/16S/ITS subset:

| Dataset ID | Public study | Intended use | Current limitation |
| --- | --- | --- | --- |
| `TRI-REAL-02` | Arabidopsis leaf microbiomes, [ENA PRJEB31530](https://www.ebi.ac.uk/ena/browser/view/PRJEB31530) | Plant phyllosphere WMS with associated 16S and ITS data [2] | Assay runs require study-specific sample pairing rather than automatic ENA alias grouping. |
| `TRI-REAL-03` | Earth Microbiome Project Multi-omics (EMP500), [ENA PRJEB42019](https://www.ebi.ac.uk/ena/browser/view/PRJEB42019), [Qiita 13114](https://qiita.ucsd.edu/study/description/13114) | Broad environmental multi-omics benchmark [3] | Amplicon targets and matched sample identities require manual curation. |
| `TRI-REAL-04` | HCHS/SOL GOLD gut microbiome study, [ENA ERP117287](https://www.ebi.ac.uk/ena/browser/view/ERP117287), [Qiita 11666](https://qiita.ucsd.edu/study/description/11666) | Supplemental human-associated amplicon stress test | It is not a complete matched WMS/16S/ITS cohort and is unsuitable as the primary cross-assay benchmark. |

## Dataset references

1. Shaffer JP, Carpenter CS, Martino C, et al. A comparison of six DNA
   extraction protocols for 16S, ITS and shotgun metagenomic sequencing of
   microbial communities. *BioTechniques*. 2022;73(1):34–46.
   [doi:10.2144/btn-2022-0032](https://doi.org/10.2144/btn-2022-0032),
   [PubMed 35713407](https://pubmed.ncbi.nlm.nih.gov/35713407/).
2. Regalado J, Lundberg DS, Deusch O, et al. Combining whole-genome shotgun
   sequencing and rRNA gene amplicon analyses to improve detection of
   microbe-microbe interaction networks in plant leaves. *The ISME Journal*.
   2020;14:2116–2130.
   [doi:10.1038/s41396-020-0665-8](https://doi.org/10.1038/s41396-020-0665-8),
   [PubMed 32405027](https://pubmed.ncbi.nlm.nih.gov/32405027/).
3. Shaffer JP, Nothias LF, Thompson LR, et al.; Earth Microbiome Project 500
   Consortium. Standardized multi-omics of Earth's microbiomes reveals
   microbial and metabolite diversity. *Nature Microbiology*.
   2022;7:2128–2150.
   [doi:10.1038/s41564-022-01266-x](https://doi.org/10.1038/s41564-022-01266-x),
   [PubMed 36443458](https://pubmed.ncbi.nlm.nih.gov/36443458/).
