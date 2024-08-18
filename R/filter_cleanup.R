#' @importFrom dplyr %>%
NULL

#' @name filter-alleles
#' @title Filter alleles to prepare data set for analysis
#'
#' @description Both `filter_isolates_browser()` and `filter_microbigge()` serve
#'   four major roles in preparing data sets for analysis, all of which are
#'   optional:
#'    * keep rows with `species` matching a value in `species`
#'    * keep rows with `allele` matching a value in `filter`
#'    * drop rows with `method` matching a value in `remove`
#'    * drop rows with duplicate `allele`/`biosample` combinations based on a logical value `deduplicate`
#'
#'   Note that filtering requiring a specific blaOXA Family in an isolates is
#'   handled separetely using `filter_oxa_family()`.
#'
#'   Additionally, `filter_microbigge()` will drop rows not meeting minimum
#'   thresholds of `coverage` and `identity`.
#'
#'   ## Allele calling method
#'   The `method` column contains codes representing the calling method used for
#'   a given allele, some of which correspond to "less than perfect" matches and
#'   may warrant removal depending on the scope and goals of the analysis.
#'
#'   More detailed explanations of these calling methods and the specific
#'   parameters used to generate them is available from NCBI at
#'   [https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories]
#'   and [https://github.com/ncbi/amr/wiki/Methods]
#'
#'   ## Coverage and identity
#'   These values are produced by AMRFinderPlus as part of the analysis pipeline
#'   and percentages use the closest sequence in the the Reference Gene Catalog
#'   used as a reference.
#'
#'   * `coverage` is the percentage of the reference sequence covered by the alignment between the target and the reference sequence
#'   * `identity` is the percentage of identical amino acids or base pairs in the alignment between the target and reference sequence
#'
#'   More information is available from NCBI at
#'   [https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories]
#'   and [https://github.com/ncbi/amr/wiki/Methods]
#'
#'   ## Deduplication
#'   Removes duplicate `allele`/`biosample` combinations, which can be included
#'   for a variety of reasons. Runs `dplyr::distinct(biosample, allele,
#'   .keep_all = TRUE)`.
#'
#' @param data A data frame or tibble with an `allele` column
#' @param species A character vector of `species` names to keep. Uses regex,
#'   vector is collapsed with "|".
#' @param filter A character vector of `allele` names to keep. Uses regex,
#'   vector is collapsed with "|".
#' @param remove A character vector of `method` values to remove. Uses regex,
#'   vector is collapsed with "|".
#' @param deduplicate Should duplicate `allele`/`biosample` combinations be
#'   removed?
#' @returns `data` with rows removed according to the selected parameters
#' @export
filter_isolates_browser <- function(data, species, filter, remove, deduplicate = TRUE) {
  has_species <- !missing(species)
  has_filter <- !missing(filter)
  has_remove <- !missing(remove)

  .complete <- data %>%
    {if(has_species) dplyr::filter(., grepl(paste(.env$species, collapse = "|"), .data$species)) else . } %>%
    {if(has_filter) dplyr::filter(., grepl(paste(filter, collapse = "|"), .data$allele)) else . } %>%
    {if(has_remove) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$method))  else . }%>%
    {if(deduplicate) dplyr::distinct(., .data$biosample, .data$allele, .keep_all = TRUE) else . }

  return(.complete)
}

#' @rdname filter-alleles
#' @param coverage Minimum percent coverage to keep
#' @param identity Minimum percent identity to keep
#' @export
filter_microbigge <- function(data, coverage = 100L, identity = 100L, species, filter, remove, deduplicate = TRUE) {
  . <- NULL # Workaround to suppress `no visible binding for global variable`

  has_species <- !missing(species)
  has_filter <- !missing(filter)
  has_remove <- !missing(remove)

  .complete <- data %>%
    dplyr::filter("coverage" >= coverage) %>%
    dplyr::filter("identity" >= identity) %>%
    {if(has_species) dplyr::filter(., grepl(paste(.env$species, collapse = "|"), .data$species)) else . } %>%
    {if(has_filter) dplyr::filter(., grepl(paste(filter, collapse = "|"), .data$allele)) else . } %>%
    {if(has_remove) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$method)) else .} %>%
    {if(deduplicate) dplyr::distinct(dplyr::arrange(., .data$protein), .data$biosample, .data$allele, .data$ipg, .keep_all = TRUE) else . }

  return(.complete)
}

#' Filter isolates by presence of blaOXA Families
#'
#' Keep only isolates containing the blaOXA family (or families) listed in the
#' character fector `family`. Filters based on BioSamples.
#'
#' @param data A data frame or tibble
#' @param family A character vector of (potentially partial) blaOXA family names
#'   to keep. Collapsed with '|' and filtered with [grepl()].
#' @export
filter_oxa_family <- function(data, family){
  has_family <- !missing(family)

  if(has_family){
    keep_biosamples <- data %>%
      filter(grepl(paste(family, collapse = "|"), oxa_family)) %>%
      distinct(biosample) %>%
      pull(biosample)

    .complete <- data %>%
      filter(biosample %in% keep_biosamples)

    return(.complete)
  } else {
    return(data)
  }
}
