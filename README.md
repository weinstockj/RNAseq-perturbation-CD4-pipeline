This repo contains the bioninformatics processing steps (primarily alignment, QC, gene quantification) and most downstream analyses 
included in [Gene regulatory network inference from CRISPR perturbations in primary CD4+ T cells elucidates the genomic basis
of immune disease](https://www.biorxiv.org/content/10.1101/2023.09.17.557749v2). 

## Requirements

The data processing pipeline uses Snakemake. You can recreate the mamba/conda environment here using
`mamba env create --file environment.yaml` . Please note that the Snakefile does contain a few 
hardcoded paths, which are written before the first rule. 

This repo also contains R code to perform several analyses. You can recreate the R envinroment
using [renv](https://rstudio.github.io/renv/index.html), i.e., 
`install.packages("renv"); renv::restore()` . The R analyses are implemented using a single large
[targets](https://books.ropensci.org/targets/) workflow file. 

## Contact

For questions about the code, please file an issue here or contact Josh Weinstock <joshua.s.weinstock@gmail.com>. 
