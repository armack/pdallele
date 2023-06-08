#' @importFrom dplyr %>%
NULL

#' @name filter-alleles
#' @title Filter alleles to prepare data set for analysis
#'
#' @description Both functions serve three major roles in preparing data sets
#'   for analysis, all of which are optional:
#'    * keep rows with `allele` matching a value in `filter`
#'    * drop rows with `method` matching a value in `remove`
#'    * drop rows with duplicate `allele`/`biosample` combinations based on a logical value `deduplicate`
#'
#'   Additionally, `filter_microbigge()` will drop rows not meeting minimum
#'   thresholds of `coverage` and `identity`.
#'
#'   ## Allele calling method The `method` column contains codes representing
#'   the calling method used for a given allele, some of which correspond to
#'   "less than perfect" matches and may warrant removal depending on the scope
#'   and goals of the analysis. Detailed explanations are available at:
#'   [https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories]
#'
#'   ## Coverage and identity These values are produced by AMRFinderPlus as part
#'   of the analysis pipeline:
#'   * `coverage` is the percentage of the reference sequence covered by the alignment between the target and the reference sequence
#'   * `identity` is the percentage of identical amino acids or base pairs in the alignment between the target and reference sequence
#'
#'   The default thresholds correspond to those used in Isolates Browser (and
#'   thus `amr.metadata.tsv` files on the NCBI FTP server). See more info at:
#'   [https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories]
#'
#'   ## Deduplication
#'   Removes duplicate `allele`/`biosample` combinations, which can be included
#'   for a variety of reasons. Runs `dplyr::distinct(biosample, allele,
#'   .keep_all = TRUE)`.
#'
#' @param data A dataframe or tibble with an `allele` column
#' @param filter A character vector of `allele` names to keep. Uses regex,
#'   vector is collapsed with "|".
#' @param remove A character vector of `method` values to remove. Uses regex,
#'   vector is collapsed with "|".
#' @param deduplicate Should duplicate `allele`/`biosample` combinations be
#'   removed?
#' @returns `data` with rows removed according to the selected parameters
#' @export
filter_isolates_browser <- function(data, filter, remove, deduplicate = TRUE) {
  has_filter <- missing(filter)
  has_remove <- missing(remove)

  .complete <- data %>%
    {if(has_filter) dplyr::filter(., grepl(paste(filter, collapse = "|"), .data$allele)) else . } %>%
    {if(has_remove) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$method))  else . }%>%
    {if(deduplicate) dplyr::distinct(., .data$biosample, .data$allele, .keep_all = TRUE) else . }

  return(.complete)
}

#' @rdname filter-alleles
#' @param coverage Minimum percent coverage to keep
#' @param identity Minimum percent identity to keep
#' @export
filter_microbigge <- function(data, coverage = 100L, identity = 90L, filter, remove, deduplicate = TRUE) {
  . <- NULL # Workaround to suppress `no visible binding for global variable`

  has_filter <- missing(filter)
  has_remove <- missing(remove)

  .complete <- data %>%
    dplyr::filter("coverage" >= coverage) %>%
    dplyr::filter("identity" >= identity)
    {if(has_filter) dplyr::filter(., grepl(paste(filter, collapse = "|"), .data$allele)) else . } %>%
    {if(has_remove) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$method)) else .} %>%
    {if(deduplicate) dplyr::arrange(., .data$protein) %>% dplyr::distinct(.data$biosample, .data$allele, .keep_all = TRUE) else . }

  return(.complete)
}
