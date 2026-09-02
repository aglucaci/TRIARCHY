# TRIARCHY manuscript

The main manuscript is formatted as a Bioinformatics Application Note:

- short four-heading abstract;
- approximately 2,000-word budget when the final Figure 1 is included;
- title includes the application name;
- availability, contact and supplementary-information fields are explicit.

## Overleaf upload

Upload the contents of this `manuscript` directory as a new Overleaf project.
Overleaf will automatically use `main.tex` as the main document and run BibTeX
for `references.bib`. The project is compatible with the default pdfLaTeX
compiler; no custom fonts, shell escape or external build steps are required.

The files required to compile the manuscript are:

- `main.tex` (the Overleaf entry point);
- `triarchy.tex` (the manuscript source);
- `references.bib`;
- `figures/figure1_workflow_architecture.pdf`.

The Markdown planning files, SVG figure sources and alternate composite figure
are optional working files and may also be uploaded safely.

Compile the manuscript locally from this directory with:

```bash
latexmk -pdf main.tex
```

Clean generated LaTeX files with:

```bash
latexmk -C main.tex
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
