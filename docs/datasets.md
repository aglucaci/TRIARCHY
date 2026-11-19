# Public empirical cohorts

Both cohorts come from the UC San Diego DNA extraction comparison study
[ENA ERP124610](https://www.ebi.ac.uk/ena/browser/view/ERP124610)
([BioProject PRJEB40903](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB40903),
[Qiita study 12201](https://qiita.ucsd.edu/study/description/12201)), which
compared DNA extraction kits using 16S, ITS, and shotgun metagenomic
sequencing of the same extractions [1].

| Cohort | Use | Samples | Controls | Ground truth |
| --- | --- | --- | --- | --- |
| `TRI-REAL-05` (default) | Empirical validation on HPC | Whole extraction plates (190 samples by default) | Every extraction blank linked to the samples on its plate | ZymoBIOMICS D6310 dilution series |
| `TRI-REAL-01` | Quick integration smoke test | 4 groups | None (empty control table) | None |

## TRI-REAL-05: extraction plates with linked blanks and mock truth

### Design

The Qiita sample information for study 12201 records, for every extraction,
the kit (`extraction_kit_round`), extraction round, extraction plate
(`biomass_plate`: high or low biomass) and plate well. `TRI-REAL-05` uses it
to build a cohort that the ENA metadata alone cannot support:

- **Whole plates.** Every sample on the selected extraction plates that has
  WMS, 16S, and ITS runs is included. By default these are the two round-1
  PowerSoil plates (the Earth Microbiome Project reference kit).
- **Linked extraction blanks.** Each extraction blank is linked, per assay, to
  every sample extracted on the same plate (`control_stage=extraction`,
  `relationship=same_batch`). Blanks are never linked across plates.
- **Physical ground truth.** Each plate carries a dilution series (`d1`–`d8`,
  with three replicates on the low-biomass plate) of the ZymoBIOMICS
  Microbial Community Standard II, Log Distribution (D6310): eight bacteria and
  two yeasts spanning 89.1% to 0.000089% of genomic DNA. These samples are
  processed as biological samples, so TRIARCHY scores them like any other, and
  `truth/mock_truth.tsv` lists their expected composition, transcribed from
  Table 1 of the Zymo instruction manual (genomic DNA, 16S, 16S+18S, genome
  copy, and cell number percentages).
- **Diverse materials.** The default plates include human saliva, feces,
  skin, urine and milk; cat and mouse feces; mouse jejunum; bare and
  rhizosphere soil; kimchi; yogurt; floor dust; keyboard swabs; freshwater; and
  seawater, mostly in technical triplicate.

| Property | Default selection (`PowerSoil r1`, both plates) |
| --- | --- |
| Samples | 190: 148 biological, 32 ZymoBIOMICS D6310, 10 extraction blanks |
| Libraries | 570 (190 each of WMS, 16S, ITS) |
| Control relationships | 2,694 blank-to-library links |
| Download | About 6.1 GB (WMS 5.7 GB, 16S 0.2 GB, ITS 0.2 GB) |
| Layout | Paired-end WMS (HiSeq 4000); single-end 16S and ITS (MiSeq) |
| Amplicons | EMP 515F/806R (16S V4) and ITS1f/ITS2 (ITS1) |

The study used Earth Microbiome Project primers, not TRIARCHY's default
V3–V4 and ITS2 configuration. `triarchy real-world` therefore applies
[`conf/cohorts/shaffer_emp.config`](../conf/cohorts/shaffer_emp.config), which
sets the EMP primers and makes ITSxpress extract ITS1. Without it, ITSxpress
would find no ITS2 region in these amplicons.

### Scaling the cohort

| Option | Samples | Libraries | Download |
| --- | --- | --- | --- |
| Default (`PowerSoil r1`) | 190 | 570 | 6.1 GB |
| `--max-biological-per-plate 20` | 82 | 246 | 2.3 GB |
| `--extraction-kit "PowerSoil r1" --extraction-kit "PowerSoil Pro"` | 377 | 1,131 | 14 GB |
| All six round-1 kits | 1,282 | 3,846 | 97 GB |

`--max-biological-per-plate` caps the non-mock biological samples on each
plate, rotating through sample types so every material stays represented;
blanks and mocks are always kept. Adding kits turns the cohort into an
extraction-method comparison. Round-1 kits are `PowerSoil r1`,
`PowerSoil Pro`, `MagMax Beta`, `Norgen`, `Zymo MagBead`, and `Homebrew`;
round-2 kits (`PowerSoil r2`, `MagMAX Microbiome`, `NucleoMag Food`) carry
the even D6300 standard instead. `Zymo MagBead` was used in both rounds and
selecting it includes both. `dataset.json` lists the complete samples available
on every plate.

### Running on an HPC cluster

1. **Download on a node with internet access** (login or transfer node). This
   fetches the ENA run table and the public Qiita sample information, writes
   the manifests and truth table, downloads the FASTQs with MD5 verification,
   and validates the cohort:

   ```bash
   mamba activate TRIARCHY
   triarchy real-world --dataset TRI-REAL-05 --download-only
   ```

   Files land in `data/real_world/TRI-REAL-05/` (`reads/`, `metadata/`,
   `truth/`). Use `--data-dir` to place them on project storage.

2. **Write a site config** for your scheduler, for example
   `/shared/project/triarchy-site.config`:

   ```groovy
   process {
       clusterOptions = '--account=YOUR_ACCOUNT'
       queue = 'YOUR_PARTITION'
       // WMS taxonomy dominates runtime: each task scans the full Kraken 2
       // PlusPFP database. Allow more concurrent tasks if your storage can
       // sustain it, and place the database on fast storage.
       withName: WGS_TAXONOMY {
           maxForks = 8
       }
   }
   ```

3. **Submit the Nextflow head job** so it survives logout. It needs little
   CPU or memory but must run for the whole workflow:

   ```bash
   #!/bin/bash
   #SBATCH --job-name=triarchy-tri-real-05
   #SBATCH --time=7-00:00:00
   #SBATCH --cpus-per-task=2
   #SBATCH --mem=8G
   #SBATCH --output=triarchy-tri-real-05-%j.log

   source "$(conda info --base)/etc/profile.d/conda.sh"
   conda activate TRIARCHY
   cd /path/to/TRIARCHY
   triarchy real-world --dataset TRI-REAL-05 --skip-download \
     --profile slurm \
     --site-config /shared/project/triarchy-site.config \
     --work-dir /shared/scratch/$USER/triarchy-tri-real-05-work \
     --outdir /shared/project/results/TRI-REAL-05
   ```

   Pass the same `--extraction-kit`, `--biomass-plate`, and
   `--max-biological-per-plate` options that you used for the download.
   `--resume` is on by default, so resubmitting the same script continues an
   interrupted run.

4. **Read the evaluation.** After a successful run, `triarchy real-world`
   scores the ZymoBIOMICS samples against `truth/mock_truth.tsv` and writes
   `mock_evaluation/` inside the results directory:

   | File | Contents |
   | --- | --- |
   | `mock_evaluation.json` | Per-assay recall overall, by expected genomic-DNA abundance (≥1%, 0.1–1%, 0.01–0.1%, <0.01%) and by dilution; evidence states of expected taxa; how often taxa absent from the standard were called biologically supported; contamination-associated calls in the other samples |
   | `mock_taxon_calls.tsv` | One row per mock sample and expected genus: detection and relative abundance in each assay, and TRIARCHY's classification |
   | `mock_unexpected_taxa.tsv` | Evidence rows in mock samples for genera not in the standard |

   To rerun the evaluation by hand:

   ```bash
   python3 scripts/evaluate_mock_truth.py \
     --results /shared/project/results/TRI-REAL-05 \
     --truth data/real_world/TRI-REAL-05/truth/mock_truth.tsv \
     --samples data/real_world/TRI-REAL-05/metadata/samples.tsv
   ```

   **Download problems.** Some ENA runs have a broken processed FASTQ (for
   example `ERR8577066_2.fastq.gz` is an empty directory on ENA's server).
   When that happens the downloader fetches that run's original submitted
   FASTQs instead, both mates together, and verifies them against their own
   ENA checksums. Runs that still fail are listed in
   `metadata/download_failures.tsv` and the command exits with an error after
   trying every other file; rerun the same command to retry them, since
   completed files are kept.

### Caveats

- Comparisons are at genus level. 16S V4 cannot reliably separate
  *Escherichia* from *Shigella*, and reference databases use different names
  for some organisms (for example *Limosilactobacillus* for *Lactobacillus
  fermentum*); the evaluator accepts these synonyms.
- 16S does not target the two yeasts and ITS does not target the bacteria, so
  each assay is scored only on the organisms it can detect.
- The truth table gives the composition of the standard as manufactured. Lot
  certificates of analysis report small deviations, and extraction bias is
  part of what this study measured, so expect differences in abundance even
  when detection is correct.
- Dilution labels (`d1`–`d8`) are taken from the Qiita sample names. DNA
  yields fall with increasing `d`, but the dilution factors are not recorded
  in the public metadata.
- The WMS libraries are shallow (on average about 30 MB of compressed reads
  per sample), which limits detection of the rarest members of D6310.

## TRI-REAL-01: small integration subset

`TRI-REAL-01` is a small subset of the same study, kept for quick integration
runs. It has no control relationships and no ground truth.

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
