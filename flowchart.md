```mermaid
flowchart TD

    A[Input Manifests<br>Paired Controls TSV] --> B[Build Unified Samples.tsv]

    B --> WGSFASTQ[WGS FASTQs]
    B --> S16FASTQ[16S FASTQs]
    B --> ITSFASTQ[ITS FASTQs]

    %% WGS branch
    WGSFASTQ --> WGSQC[WGS QC<br>(fastp)]
    WGSQC --> HOSTREM[Host Removal<br>(Bowtie2)]
    HOSTREM --> HUMANN[HUMAnN Functional Profiling]
    HOSTREM --> KRBRACK[Kraken2 + Bracken<br>WGS Taxonomy]
    KRBRACK --> WGSnorm[WGS Taxonomy Normalization]

    %% 16S branch
    S16FASTQ --> S16QC[16S QC<br>(cutadapt)]
    S16QC --> S16DADA[DADA2 Denoising]
    S16DADA --> S16TAX[16S Taxonomy]
    S16TAX --> S16TABLE[16S ASV Table]

    %% ITS branch
    ITSFASTQ --> ITSQC[ITS QC<br>(cutadapt)]
    ITSQC --> ITSx[ITS Extraction<br>(ITSxpress)]
    ITSx --> ITSDADA[ITS DADA2]
    ITSDADA --> ITSTAX[ITS Taxonomy]
    ITSTAX --> ITSTABLE[ITS ASV Table]

    %% Controls + Diversity
    WGSnorm --> CTRL[Controls + Negatives<br>Pair QC]
    S16TABLE --> CTRL
    ITSTABLE --> CTRL
    HUMANN --> CTRL

    CTRL --> DIV[Alpha/Beta Diversity]

    DIV --> REPORT[Final TRINITY Report<br>HTML + Figures]
    WGSnorm --> REPORT
    S16TABLE --> REPORT
    ITSTABLE --> REPORT
    HUMANN --> REPORT
```
