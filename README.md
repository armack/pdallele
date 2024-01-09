# pdallele

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/642998583.svg)](https://zenodo.org/doi/10.5281/zenodo.10476134)
<!-- badges: end -->

The goal of pdallele is to simplify the downloading, cleanup, and analysis of
[NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) data,
primarailly for exploring the diversity and distribution antimicrobial resistance
alleles across large numbers of bacterial isolates.

## Installation

You can install the development version of pdallele from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("armack/pdallele")
```

No CRAN release in planned.

## Example

This contains all the code the code used to generate the tables and figures in
"β-Lactamase Diversity in *Acinetobacter baumannii* and *Pseudomonas
aeruginosa*: A Bioinformatic Analysis" by Mack et al.

Data downloading, processing, and cleanup was performed as outlined in
`vignette("downloading-data", package = "pdallele")` and
`vignette("processing-data", package = "pdallele")`.

The code used for analysis and generation of tables and figures is included in
`vignette("sample-analysis"", package = "pdallele")` and
`vignette("additional-analysis", package = "pdallele")`.

## Citation
This code is available on Zenodo at https://zenodo.org/doi/10.5281/zenodo.10476134 and by DOI 10.5281/zenodo.10476135

Citation information for the preprint and final versions of the manuscript will added when available
