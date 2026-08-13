# TRIARCHY manuscript

The main manuscript is formatted as a Bioinformatics Application Note:

- short four-heading abstract;
- approximately 2,000-word budget when the final Figure 1 is included;
- title includes the application name;
- availability, contact and supplementary-information fields are explicit.

Compile the manuscript from this directory with:

```bash
latexmk -pdf triarchy.tex
```

Clean generated LaTeX files with:

```bash
latexmk -C triarchy.tex
```

Generated manuscript build products, including PDF, DOCX and LaTeX auxiliary
files, are ignored by Git. Keep the manuscript source in `triarchy.tex`,
`references.bib`, figure sources and supplementary source files.

Files requiring author input or frozen-release data:

- `triarchy.tex`: red `\todo{...}` placeholders for public URLs, DOI, contact,
  funding and co-authors;
- `figure1_plan.md`: panel specification for the single composite figure;
- `supplementary_outline.md`: required contents for the single supplementary
  PDF, including frozen-release benchmark, empirical validation, runtime and
  portability values.

Do not submit until all placeholders are replaced with public, stable resources
and the Results section includes measured accuracy from the frozen release plus
at least one real matched mock-community or environmental dataset.
