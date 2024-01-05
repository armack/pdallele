# pdallele

<!-- badges: start -->
<!-- badges: end -->

The goal of pdallele is to simplify the downloading, cleanup, and analysis of
[NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) data,
primarailly for exploring antimicrobial resistance.

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
