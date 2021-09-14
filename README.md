# R-laslo
Data cleaning and analysis of LASLO (aka LoopMatcher) results.

## Content
- data:
  *.lst files: Gene lists used to download Genbank sequences from NCBI. 
  *.txt files: Original gene lists described in papers.
  stem15: Loopmatcher results in CSV format.

- src:
  main.R : Data analysis of LoopMatcher results.
  configuration.R : Setting for functions and libraries required in the main script.
  genbank_accessions.R : Script to obtain the accession IDs in GenBank given a list of genes.
  plot_functions.R : Functions to generate plots in main script.
  count.R : Counts of genes and transcripts in several results. Unused.
  - misc:
    A collection of scripts of differents analysis, not used in the final version of the thesis.

- markdown:
  Markdown document generated with R. Unfinished.

