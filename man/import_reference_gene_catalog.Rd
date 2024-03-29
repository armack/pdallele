% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_process.R
\name{import_reference_gene_catalog}
\alias{import_reference_gene_catalog}
\title{Import Reference Gene Catalog to tibble}
\usage{
import_reference_gene_catalog(path)
}
\arguments{
\item{path}{Path to an NCBI ReferenceGeneCatalog.txt file}
}
\value{
A \code{\link[=tibble]{tibble()}}.  Values are described in Details.
}
\description{
Imports Reference Gene Catalog, renames columns, and drops data not needed
for downstream processing
}
\details{
The following columns are included in the output:\tabular{ll}{
   \strong{Source File} \tab \strong{Output Tibble} \cr
   allele \tab allele \cr
   gene_family \tab gene \cr
   product_name \tab name \cr
   type \tab type \cr
   subtype \tab subtype \cr
   class \tab class \cr
   subclass \tab subclass \cr
   refseq_protein_accession \tab protein \cr
   blacklisted_taxa \tab blacklisted_taxa \cr
   whitelisted_taxa \tab whitelisted_taxa \cr
}
}
