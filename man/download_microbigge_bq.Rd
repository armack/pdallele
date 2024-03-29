% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{download_microbigge_bq}
\alias{download_microbigge_bq}
\title{Download MicroBIGG-E Data from the Google Cloud Platform}
\usage{
download_microbigge_bq(
  billing,
  path,
  taxgroup,
  sciname,
  element,
  debug = FALSE
)
}
\arguments{
\item{billing}{GCP project identifier for billing/usage tracking}

\item{path}{Directory to save downloaded files to}

\item{taxgroup}{An NCBI Pathogen Detection organism group to filter by}

\item{sciname}{Scientific name of the organism of interest to filter by}

\item{element}{Name of a gene/allele to filter by}

\item{debug}{Should the query be previewed rather than executed?}
}
\description{
Download the latest MicroBIGG-E data from the Google Cloud
Platform using \code{bigrquery} with optional filtering by \code{taxgroup}, \code{sciname},
and \code{element}.

The downloaded data is saved as 'microbigge.tsv' in \code{path}.

Note that \code{taxgroup}, \code{sciname}, and \code{element} all use SQL syntax for
queries, meaning '\%' will be treated as a wildcard rather than '*'.

If you are not familiar with \code{bigrquery}, please see
\link{https://bigrquery.r-dbi.org/index.html#important-details} for more
information on authentication and generating a billing project code.
}
\seealso{
Other download data: 
\code{\link{download_identical_protein_groups}()},
\code{\link{download_metadata_clusters}()},
\code{\link{download_reference_gene_catalog}()}
}
\concept{download data}
