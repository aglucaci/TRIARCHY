# TRINITY
Three sequencing modalities (WGS, 16S, ITS), one unified workflow for environmental microbiome analysis.




| Step                 | WGS             | 16S    | ITS    |
| -------------------- | --------------- | ------ | ------ |
| QC                   | yes             | yes    | yes    |
| Host Removal         | **yes**         | **NO** | **NO** |
| DADA2                | no              | yes    | yes    |
| ITSxpress            | no              | no     | yes    |
| Taxonomy             | Kraken2/Bracken | SILVA  | UNITE  |
| Functional Profiling | HUMAnN          | no     | no     |
| Diversity            | yes             | yes    | yes    |
