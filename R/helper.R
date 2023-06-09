#' @importFrom dplyr %>%
#' @importFrom rlang :=
NULL

########## Flow and Piping ##########

#' Pass data along the pipeline and inform
#'
#' Pass data along the pipeline and invoke rlang::inform() with the desired message
#'
#' @param data a dataframe or tibble
#' @param ... parameters to pass to inform()
#' @returns unaltered `data`
#' @seealso [pipe_warn()]
#' @export

pipe_inform <- function(data, ...){
  rlang::inform(...)
  return(data)
}

#' Pass data along the pipeline and warn
#'
#' Pass data along the pipeline and invoke rlang::warn() with the desired message
#'
#' @param data a dataframe or tibble
#' @param ... parameters to pass to warn()
#' @returns unaltered `data`
#' @seealso [pipe_inform()]
#' @export

pipe_warn <- function(data, ...){
  rlang::warn(...)
  return(data)
}

########## Counting ##########

#' Count by columns
#'
#' @description
#' Wraps [dplyr::count()] and adds additional options for ordering of `NA`
#' values and counting distinct BioSamples
#'
#' Two convenience functions are also available:
#' * `count_alleles()` wraps `count_by_column(name = "alleles", isolates = FALSE)`
#' * `count_isolates()` wraps `count_by_column(name = "isolates", isolates = TRUE)`
#'
#' @param data A data frame or tibble.
#' @param ... <data masking> columns to group by.
#' @param sort If TRUE, will show the largest groups at the top.
#' @param name The name of the new column in the output.
#' @param na_last Where to sort `NA`. If `TRUE`, `NA` is put last; if `FALSE`,
#'   `NA` is put first; if `NA`, `NA` is dropped.
#' @param isolates Should isolates (BioSamples) be counted instead of alleles?
#' @export
count_by_column <- function(data, ..., isolates = FALSE, sort = TRUE,
                            name = "n", na_last = TRUE){
  . <- NULL # Workaround to suppress `no visible binding for global variable`

  .complete <- data %>%
    {if (isolates) dplyr::distinct(., .data$biosample, .keep_all = TRUE) else . } %>%
    dplyr::count(..., name = name) %>%
    {if(sort) dplyr::mutate(., dplyr::across(c(...), ~forcats::fct_reorder(., !!rlang::sym(name), .desc = TRUE))) %>%
        dplyr::arrange(., dplyr::across(c(...))) else . } %>%
    {if(identical(na_last, TRUE)) dplyr::arrange(., dplyr::across(c(...), ~ is.na(.)))
      else if(identical(na_last, FALSE)) dplyr::arrange(., dplyr::across(c(...), ~ !is.na(.)))
      else if(is.na(na_last)) tidyr::drop_na(., ...)
      else dplyr::arrange(., dplyr::across(c(...), ~ is.na(.))) %>%
        pipe_inform(., "Invalid `na_last` provided. Defaulted to TRUE.") }

  return(.complete)
}

#' Count alleles by columns
#' @rdname count_by_column
#' @export
count_alleles <- function(data, ..., sort = TRUE, name = "alleles", na_last = TRUE){
  .complete <- count_by_column(data = data, ..., sort = sort,
                               isolates = FALSE, name = name,
                               na_last = na_last)
    return(.complete)
}

#' Count isolates by columns
#' @rdname count_by_column
#' @export
count_isolates <- function(data, ..., sort = TRUE, na_last = TRUE){
  .complete <- count_by_column(data = data, ..., sort = sort, isolates = TRUE,
                  name = "isolates", na_last = na_last)
    return(.complete)
}

########## Releveling ##########

#' Relevel and sort factors in human language order
#'
#' Relevel factors in human language alphanumeric order (i.e. digits are sorted
#' numerically rather than as strings) and sort columns in ascending order.
#'
#' @details
#' Uses [stringr::str_sort()] with `numeric = TRUE` for sorting.
#'
#' @param data A dataframe or tibble
#' @param ... <data-masking> Columns to relevel alphanumerically
#' @family factor releveling
#' @export
relevel_numeric <- function(data, ...) {

  .complete <- data %>%
    dplyr::mutate(dplyr::across(c(...), ~forcats::fct_relevel(., ~ stringr::str_sort(., numeric = TRUE)))) %>%
    dplyr::arrange(...)

  return(.complete)
}

#' Relevel strings to first or last position and sort
#'
#' Relevel values in `...` matching values in `first` or `last` to the beginning
#' or end (respectively) of the factor levels for `...`.
#'
#' If `first` or `last` contain multiple values, they will remain in the order
#' provided as input.
#'
#' @param data A data frame or tibble
#' @param ... <data masking> Columns to reorder
#' @param first Character vector of values to move to the first position
#' @param last Character vector of values to move to the last position
#' @param na_last Where should `NA` go? `TRUE` at the end, `FALSE` at the
#'   beginning, `NA` dropped.
#' @family factor releveling
#' @export
relevel_first_last <- function(data, ..., first = NULL, last = NULL, na_last = TRUE){

  .complete <- data %>%
    dplyr::mutate(dplyr::across(c(...), ~forcats::fct_expand(., cur_column(), first, last))) %>%
    dplyr::mutate(dplyr::across(c(...), ~forcats::fct_relevel(., first, after = 0))) %>%
    dplyr::mutate(dplyr::across(c(...), ~forcats::fct_relevel(., last, after = Inf))) %>%
    dplyr::arrange(...) %>%
    .resolve_na_last(..., na_last = na_last)

  return(.complete)
}

########## Filtering ##########

#' Filter *bla* alleles by assignment status
#'
#' @description These functions serve to filter beta-lactamase (*bla*) alleles
#'   based on their assignment status. Four variants are available:
#'   * `filter_assigned_bla()` keeps only assigned *bla* alleles from `data`
#'   * `filter_unassigned_bla()` keeps only unassigned *bla* alleles from `data`
#'   * `remove_assigned_bla()` removes assigned *bla* alleles from `data` while leaving other alleles untouched
#'   * `remove_unassigned_bla()` removes unassigned *bla* alleles from `data` while leacing other alleles untouched
#'
#' @details In NCBI Pathogen Detection Project datasets, *bla* alleles with
#'   formal designations are listed by their number (e.g. *bla*PDC-3) while
#'   those without formal designations are collectively listed by their gene
#'   name (e.g. *bla*PDC).
#'
#'   Every *bla*PDC-3 allele will have an identical protein sequence and resolve
#'   to the same Identical Protein Group accession number. The same is NOT true
#'   for *bla*PDC, as this is a generic representation of all unassigned
#' *bla*PDC alleles and can theoretical be associated with hundreds or even
#'   thousands of unique protein sequences.
#'
#'   Note that nomenclature varies and 'unassigned' alleles are often referred
#'   to as 'novel' alleles.
#'
#' @param data A data frame or tibble
#' @returns A data frame or tibble with either assigned (`remove_assigned_bla`)
#'   or unassigned (`remove_unassigned_bla`) alleles removed
#' @export
remove_assigned_bla <- function(data){
  .complete <- data %>%
    dplyr::filter(!(grepl("bla", .data$allele, ignore.case = TRUE) & grepl("-", .data$allele)))

  return(.complete)
}

#' @rdname remove_assigned_bla
#' @export
remove_unassigned_bla <- function(data){
  .complete <- data %>%
    dplyr::filter(!(grepl("bla", .data$allele, ignore.case = TRUE) & !grepl("-", .data$allele)))

  return(.complete)
}

#' @rdname remove_assigned_bla
#' @export
filter_assigned_bla <- function(data){
  .complete <- data %>%
    dplyr::filter((grepl("bla", .data$allele, ignore.case = TRUE) & grepl("-", .data$allele)))

  return(.complete)
}

#' @rdname remove_assigned_bla
#' @export
filter_unassigned_bla <- function(data){
  .complete <- data %>%
    dplyr::filter((grepl("bla", .data$allele, ignore.case = TRUE) & !grepl("-", .data$allele)))

  return(.complete)
}

#' Filter attribute column by string
#'
#' A convenience function wrapping `filter(grepl(string, attribute, ignore.case =
#' TRUE))` which looks cleaner in pipes and follows tidyverse argument order
#'
#' @param data A data frame or tibble
#' @param attribute <data-masking> Column to filter
#' @param string String to look for
#' @export
filter_attribute <- function(data, attribute, string){

  .complete <- data %>%
    dplyr::filter(grepl(string, {{attribute}}, ignore.case = TRUE ))

  return(.complete)
}

#' Filter by allele exclusivity
#'
#' Filter `data` to keep alleles based on whether they are exclusive to one
#' group or shared between multiple groups. Two variants are available:
#'   * `filter_exclusive_alleles()` keeps alleles appearing only in one group of `group_by(data, ...)`
#'   * `filter_shared_alleles()` keeps alleles appearing in multiple groups of `group_by(data, ...)`
#'
#' @param data A data frame or tibble
#' @param ... <data-masking> Columns to determine exclusivity by
#' @returns A data frame or tibble with only exclusive
#'   (`filter_exclusive_alleles`) or shared (`filter_shared_alleles`) alleles
#' @export
filter_exclusive_alleles <- function(data, ...){

  .complete <- data %>%
    dplyr::group_by(.data$allele) %>%
    dplyr::filter(length(unique( !!!rlang::enquos(...) )) == 1) %>%
    dplyr::ungroup()

  return(.complete)
}

#' @rdname filter_exclusive_alleles
#' @export
filter_shared_alleles <- function(data, ...){

  .complete <- data %>%
    dplyr::group_by(.data$allele) %>%
    dplyr::filter(length(unique( !!!rlang::enquos(...) )) > 1) %>%
    dplyr::ungroup()

  return(.complete)
}


#' Determine allele combinations
#'
#' Determine combinations of `allele` present in each isolate, concatenate a
#' list of alleles into a `combo` column, and remove isolates with only a single
#' allele.
#'
#' @param data A data frame or tibble
#' @param filter_terms A character vector of (potentially partial) allele names
#'   to keep. Collapsed with '|' and filtered with [grepl()].
#' @export
determine_combinations <- function(data, filter_terms = NULL){
  . <- NULL # Workaround to suppress `no visible binding for global variable`
  .combos <- data %>%
    {if (!is.null(filter))  dplyr::filter(., grepl(paste(filter_terms, collapse = "|"), .data$allele, ignore.case = TRUE )) else .} %>%
    dplyr::group_by(.data$biosample) %>%
    dplyr::filter(length(unique(.data$allele)) > 1) %>%
    dplyr::summarize(combo = paste(.data$allele, collapse = ", "), .groups = "drop")

  .complete <- data %>%
    inner_join(.combos, by = "biosample")

  return(.complete)
}

########## Split Integers into Groups ##########


#' Define boundaries to split integer vector into groups of roughly equal size
#'
#' Used to define boundaries to split year into periods
#'
#' Companion function to categorize_integer_groups()
#'
#' @param values integer vector of values to group
#' @param groups integer of groups to define
#' @export
determine_nearly_equal_integer_groups_lt <- function(values, groups = 4){
  if(!is.integer(values)){
    stop("Please ensure 'values' is an integer vector.")
  }
  quantiles <- values %>%
    stats::quantile(probs = seq(from = 0, to = 1, length.out = groups + 1), type = 1)

  value_boundaries <- vector(mode = "integer", length = length(quantiles) )
  value_boundaries[1] = quantiles[[1]]

  for(i in 2:(length(quantiles) - 1)){
    last_value <- value_boundaries[[i - 1]]
    threshold = length(values) / groups
    test_value <- quantiles[[i]]

    tibble <- tibble::tibble(values) %>%
      dplyr::filter(values >= last_value)

    include <- tibble %>%
      dplyr::filter(values <= test_value ) %>%
      nrow()

    exclude <- tibble %>%
      dplyr::filter(values < test_value ) %>%
      nrow()

    if(abs(threshold - include) < abs(threshold - exclude) ){
      next_value <- test_value + 1
    }else{
      next_value <- test_value
    }
    value_boundaries[i] <- next_value
  }
  value_boundaries[length(value_boundaries)] = quantiles[[length(quantiles)]] + 1

  if(any(duplicated(value_boundaries))){
    stop("Duplicate boundaries. Ensure `data` contains sufficient rows and consider tring fewer `groups`.")
  }

  return(value_boundaries)
}

#' Assign column values to groups
#'
#' Companion function to determine_nearly_equal_integer_groups_lt()
#'
#' Either groups or split should be provided. Split takes precedence.
#'
#' @param tibble A data frame or tibble
#' @param source_col  integer column to group
#' @param label_col  (new) column for group names
#' @param label_sep seperator between upper and lower bondary in lable
#' @param groups number integer of groups to define (fed into
#'        determine_nearly_equal_integer_groups_lt() to determine boundaries)
#' @param split integer vector of boundaries, which is generally the output of
#'        determine_nearly_equal_integer_groups_lt()
#' @export

categorize_integer_groups <- function(data, source_col, label_col, label_sep = " - ", groups, split = NULL){
  source_col_sym = rlang::sym(source_col)
  label_col_sym = rlang::sym(label_col)

  if( is.null(split) ){
    bounds <- data %>%
      tidyr::drop_na(!!source_col) %>%
      dplyr::group_by(.data$biosample) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(!!source_col) %>%
      determine_nearly_equal_integer_groups_lt(groups)
  }

  lower_bounds <- bounds[1:length(bounds)-1]
  upper_bounds <- bounds[2:length(bounds)]
  bounds_labels <- paste(lower_bounds,label_sep,upper_bounds - 1, sep = "")
  bounds_conditions <- paste(source_col," >= ",lower_bounds," & ",source_col," < ",upper_bounds," ~ ","'",bounds_labels,"'",sep="")

  data %>%
    dplyr::mutate(!!label_col := dplyr::case_when(!!!rlang::parse_exprs(bounds_conditions))) %>%
    return()
}

########## Other ##########

#' Process na_last arguments
#'
#' Used internally to avoid replicating na_last logic across multiple functions
#'
#' @param data a dataframe or tibble
#' @param ... <data-masking> columns to reorder NAs
#' @param na_last argument to be processed (pass directly from parent function)
#' @export
.resolve_na_last <- function(data, ..., na_last){
  . <- NULL # Workaround to suppress `no visible binding for global variable`

  .complete <- data %>%

    {if(identical(na_last, TRUE)) dplyr::arrange(., dplyr::across(c(...), ~ is.na(.)))
      else if(identical(na_last, FALSE)) dplyr::arrange(., dplyr::across(c(...), ~ !is.na(.)))
      else if(is.na(na_last)) tidyr::drop_na(., ...)
      else dplyr::arrange(., dplyr::across(c(...), ~ is.na(.))) %>%
        pipe_inform(., "Invalid `na_last` provided. Defaulted to TRUE.") }

  return(.complete)
}

#' Apply separate_rows sequentially by column
#'
#' @description
#' Separate and all unique combinations of separated rows when
#' different columns may have different numbers of delimiters for the same row
#' (which causes a vector recycling error when applying separate_rows to
#' multiple columns). This essentially "expands" the data column by column.
#'
#' This approach is relatively slow, but serves the purpose without needing to
#' modify separate_rows directly.
#'
#' @param data a dataframe or tibble
#' @param ... <tidy-select> columns to apply separate_rows to
#' @param sep Separator delimiting collapsed values.
#' @param convert If TRUE will automatically run type.convert() on the key
#'   column. This is useful if the column types are actually numeric, integer,
#'   or logical.
#' @export
separate_rows_sequential <- function(data, ..., sep = "[^[:alnum:].]+", convert = FALSE){
  vars <- tidyselect::eval_select(expr(c(...)), data)
  .modified <- data

  for(name in names(vars)){
    .modified <- tidyr::separate_rows(.modified, any_of(name), sep = sep, convert = convert)
  }

  return(.modified)
}

#' Determine all possible unique proteins
#'
#' @description Uses `protein` accession numbers and `allele` names to create a
#'   character vector representing the complete subset of *possible* unique
#'   alleles in the datasets.
#'
#' @details This function is only for use with MicroBIGG-E datasets, as Isolates
#'   Browser datasets do not contain sufficient data to make the determination
#'   of unique unassinged alleles.
#'
#'   Every distinct assigned allele with a "WP_" prefixed accession number is
#'   deduplicated. This does NOT alter values in the original `data` and some
#'   downstream matching may require Identical Protein Groups data.
#'
#'   Every distinct accession number corresponding to an unassigned allele is
#'   included as they correspond to *possible* distinct alleles. In most
#'   datasets, there will likely be overlap with several of these accession
#'   numbers corresponding to the same distinct protein. Downstream use may
#'   require deduplication, which can be achieved using Identical Protein Groups
#'   data.
#'
#' @param data A dataframe or tibble with columns `allele` and `protein`
#' @returns A character vector of accession numbers
#' @export
possible_unique_proteins <- function(data){
  wp_accessions <- data %>%
    filter_assigned_bla() %>%
    dplyr::filter(grepl("WP_", .data$protein)) %>%
    dplyr::pull("protein") %>%
    unique()

  wp_alleles <- data %>%
    filter_assigned_bla() %>%
    dplyr::filter(grepl("WP_", .data$protein)) %>%
    dplyr:: pull("allele") %>%
    unique()

  other_accessions <- data %>%
    tidyr::drop_na("protein") %>%
    dplyr::filter(!.data$allele %in% wp_alleles) %>%
    dplyr::pull("protein") %>%
    unique()

  accessions <- c(wp_accessions, other_accessions)

  return(accessions)
}

########## File Handling ##########

#' Ensure directory exists
#'
#' Checks if the directory path provided exists and creates it (with
#' notification) if it does not.
#'
#' @param path Directory path
#' @export
ensure_directory <- function(path){
  if(!dir.exists(path)){
    dir.create(path, recursive = TRUE)
    message(paste("Directory created:", path))
  }
}
