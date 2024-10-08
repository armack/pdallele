% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{download_identical_protein_groups}
\alias{download_identical_protein_groups}
\title{Download Identical Protein Groups data using NCBI Entrez}
\usage{
download_identical_protein_groups(accessions, n = 10, path)
}
\arguments{
\item{accessions}{A character vector of protein accession numbers}

\item{n}{Number of accessions to include in each request}

\item{path}{Destination directory}
}
\description{
Download Identical Protein Group (IPG) details for a list of protein
accessions using the NCBI Entrez API via the rentrez package. The information
for \code{n} accessions is downloaded in each batch and the remaining accessions
filtered to reduce duplicate downloads.

The downloaded data is saved as 'ipg.tsv' in \code{path}.

A lower value for \code{n} will require more queries to NCBI servers but will
reduce redundant downloads if an IPG is associated with a large number of
accession numbers in the dataset.
}
\seealso{
Other download data: 
\code{\link{download_metadata_clusters}()},
\code{\link{download_microbigge_bq}()},
\code{\link{download_reference_gene_catalog}()}
}
\concept{download data}
