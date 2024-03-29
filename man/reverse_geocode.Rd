% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_process.R
\name{reverse_geocode}
\alias{reverse_geocode}
\title{Parse \code{lat_lon} data and reverse geocode to \code{location}}
\usage{
reverse_geocode(data)
}
\arguments{
\item{data}{A data frame or tibble to containing a \code{lat_lon} column}
}
\value{
\code{data} with column \code{lat_lon} replaced by columns \code{lat} and \code{lon} and
missing values for \code{location} added where possible
}
\description{
Separates \code{lat_lon} and convert 'N/S/E/W' labels to '+/-'
format. Uses these values to reverse geocode and determine country when
\code{location} is not provided in the original metadata.

Also removes (and reverse geocodes if possible) any \code{location} values
containing a forward slash. These appear to sometimes be used in metadata
to indicate samples tied to multiple locations (e.g. a hospital in one
country treating a patient who acquired an infection in another country),
but is not part of the INSDC standard and is not consistently reliable.
}
