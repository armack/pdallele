#' @importFrom dplyr %>%
NULL

#' Clean and filter allele listings
#'
#' @description The column `allele_call` contains codes representing to the
#'   calling method used for a given allele. Some of these codes correspond to
#'   "less than perfect" matches and may warrant removal depending on the scope
#'   and goals of the analysis. Detailed explanations are available at:
#'   <https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories>
#'
#'   This function removes rows with `allele_call` matching a value in `remove`.
#'
#'   This function also (optionally) keeps only rows with `allele` matching a
#'   value in `filter`, such as when only a specific subset of `allele` is
#'   needed in the dataset.
#'
#' Optionally, can deduplicate `allele`/`biosample` combinations. This can occur
#' when an isolate contains a large number of contigs resulting in a given
#' allele being listed several times for the isolate. Runs
#' `dplyr::distinct(biosample, allele, .keep_all = TRUE)`.
#'
#' @param data A dataframe or tibble with an `allele` column
#' @param filter A character vector of `allele` names to keep. Uses grepl().
#'   Partial matches are kept. Vector is collapsed with "|" (or operator).
#' @param remove A character vector of allele_call values to remove
#' @param deduplicate Should duplicate `allele`/`biosample` combinations be
#'   removed?
#' @returns `data` with rows removed according to the selected parameters
#' @export

clean_filter_alleles <- function(data, filter, remove, deduplicate = TRUE) {
  missing_filter <- missing(filter)
  missing_remove <- missing(remove)

  .complete <- data %>%
    {if(!missing_filter) dplyr::filter(., grepl(paste(filter, collapse = "|"), .data$allele)) else . } %>%
    {if(!missing_remove) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$allele_call))  else . }%>%
    {if(deduplicate) dplyr::distinct(., .data$biosample, .data$allele, .keep_all = TRUE) else . }

  return(.complete)
}


#' Filter imported MicroBigg-E data
#'
#' @description Removes allele calls that do not meet the minimum `coverage` and
#'   `identity` thresholds. These values may vary based on the needs of a
#'   particular analysis. Set both to `100` to allow only exact allele calls.
#'
#'   Specific allele calling methods can be excluded from the dataset by name in
#'   `remove`. See <<>> for more information about allele calling methods used
#'   in MicroBIGG-E data.
#'
#' @details Isolates Browser (and thus amr.metadata.tsv files on the NCBI FTP
#'   server) uses a 90% `identity` threshold:
#'   https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories
#'
#' @param data A dataframe or tibble
#' @param coverage Minimum percent coverage to keep
#' @param identity Minimum percent identity to keep
#' @param remove Character vector of `method` values to remove alleles
#' @returns `data` with rows removed according to the selected parameters
#' @export
filter_microbigge <- function(data, coverage = 100L, identity = 90L, remove) {
  . <- NULL # Workaround to suppress `no visible binding for global variable`
  remove_missing <- missing(remove)
  .complete <- data %>%
    tidyr::drop_na("protein") %>%
    {if(!remove_missing) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$method)) else .} %>%
    dplyr::filter("coverage" >= coverage) %>%
    dplyr::filter("identity" >= identity)

  return(.complete)
}
