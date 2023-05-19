# Where should NA go? TRUE at the end, FALSE at the beginning, NA dropped.
########## Flow and Piping ##########

#' Pass data along the pipeline and rlang::inform()
#'
#' @param .data a dataframe or tibble
#' @param ... parameters to pass to inform()

pipe_inform <- function(.data, ...){
  rlang::inform(...)
  return(.data)
}

#' Pass data along the pipeline and rlang::warn()
#'
#' @param .data a dataframe or tibble
#' @param ... parameters to pass to warn()

pipe_warn <- function(.data, ...){
  rlang::warn(...)
  return(.data)
}

########## Counting ##########

#' Count by Columns
#'
#' @param .data A data frame or tibble.
#' @param ... <data masking> columns to group by.
#' @param sort If TRUE, will show the largest groups at the top.
#' @param name The name of the new column in the output.
#' @param na_last Where to sort NA. If TRUE, NA is put last; if FALSE,
#' NA is put first; if NA, NA is dropped.
#' @param isolates Should isolates (BioSamples) be counted instead of alleles?
count_by_column <- function(.data, ..., isolates = FALSE, sort = TRUE,
                            name = "n", na_last = TRUE){

  .complete <- .data %>%
    {if (isolates) distinct(., biosample, .keep_all = TRUE) else . } %>%
    count(..., name = name) %>%
    {if(sort) mutate(., across(c(...), ~fct_reorder(., !!rlang::sym(name), .desc = TRUE))) %>%
        arrange(., across(c(...))) else . } %>%
    {if(identical(na_last, TRUE)) arrange(., across(c(...), ~ is.na(.)))
      else if(identical(na_last, FALSE)) arrange(., across(c(...), ~ !is.na(.)))
      else if(is.na(na_last)) drop_na(., ...)
      else arrange(., across(c(...), ~ is.na(.))) %>%
        pipe_inform(., "Invalid `na_last` provided. Defaulted to TRUE.") }

  return(.complete)
}

#' Count Alleles by Columns
#'
#' This is a convenience function for count_by_column(name = "alleles",
#' isolates = FALSE)
#'
#' @param .data A data frame or tibble.
#' @param ... <data masking> columns to group by.
#' @param sort If TRUE, will show the largest groups at the top.
#' @param na_last Where to sort NA. If TRUE, NA is put last; if FALSE,
#' NA is put first; if NA, NA is dropped.
count_alleles <- function(.data, ..., sort = TRUE, name = "alleles", na_last = TRUE){
  .complete <- count_by_column(.data = .data, ..., sort = sort,
                               isolates = FALSE, name = name,
                               na_last = na_last)
    return(.complete)
}

#' Count Isolates by Columns
#'
#' This is a convenience function for
#' count_by_column(isolates = TRUE, name = "isolates")
#'
#' @param .data A data frame or tibble.
#' @param ... <data masking> columns to group by.
#' @param sort If TRUE, will show the largest groups at the top.
#' @param na_last Where to sort NA. If TRUE, NA is put last; if FALSE,
#' NA is put first; if NA, NA is dropped.
count_isolates <- function(.data, ..., sort = TRUE, na_last = TRUE){
  .complete <- count_by_column(.data = .data, ..., sort = sort, isolates = TRUE,
                  name = "isolates", na_last = na_last)
    return(.complete)
}

########## Releveling ##########

#' Relevel factors in human readable order and sort columns
#'
#' Relevel factors in 'human' alphanumeric order(i.e. digits are sorted
#' numerically rather than as strings) and sort columns
#'
#' @param .data A dataframe or tibble
#' @param ... <data-masking> cols to relevel alphanumerically
relevel_numeric <- function(.data, ...) {

  .complete <- .data %>%
    mutate(across(c(...), ~fct_relevel(., ~ str_sort(., numeric = TRUE)))) %>%
    arrange(...)

  return(.complete)
}


#' Relevel to first/last position and sort
#'
#' Relevel character vector values to first or last position(s) in `...` columns
#'
#' Note that both `first` and `last` will remain in the order provided
#'
#' @param .data A data frame or tibble
#' @param ... <data masking> columns to reorder
#' @param first a character vector of values to move to the first position
#' @param last a character vector of values to move to the last position
#' @param na_last Where should NA go? TRUE at the end, FALSE at the beginning, NA dropped.
relevel_first_last <- function(.data, ..., first = NULL, last = NULL, na_last = TRUE){

  .complete <- .data %>%
    mutate(across(c(...), ~fct_expand(., cur_column(), first, last))) %>%
    mutate(across(c(...), ~fct_relevel(., first, after = 0))) %>%
    mutate(across(c(...), ~fct_relevel(., last, after = Inf))) %>%
    arrange(...) %>%
    .resolve_na_last(..., na_last = na_last)

  return(.complete)
}

########## Filtering ##########

#' Remove assigned bla alleles
#'
#' Filter .data to remove unassigned bla alleles
#' (i.e. those without an allele number)
#'
#' #' @param .data A data frame or tibble
remove_assigned_bla <- function(.data){
  .complete <- .data %>%
    filter(!(grepl("bla", allele, ignore.case = TRUE) & grepl("-", allele)))

  return(.complete)
}

#' Remove unassigned bla alleles
#'
#' Filter .data to remove unassigned bla alleles
#' (i.e. those with an allele number)
#'
#' @param .data A data frame or tibble
remove_unassigned_bla <- function(.data){
  .complete <- .data %>%
    filter(!(grepl("bla", allele, ignore.case = TRUE) & !grepl("-", allele)))

  return(.complete)
}

#' Filter attribute column by string
#'
#' Shorthand for filter(grepl(string, attribute, ignore.case = TRUE))
#' which looks better in pipes and follows tidyverse argument order
#'
#' @param .data A data frame or tibble
#' @param attribute <data-masking> column to filter
#' @param string string to filter by
filter_attribute <- function(.data, attribute, string){

  .complete <- .data %>%
    filter(grepl(string, {{attribute}}, ignore.case = TRUE ))

  return(.complete)
}

#' Filter alleles exclusive to only a single combination of columns
#'
#' @param .data A data frame or tibble
#' @param ... <data-masking> cols to find alleles unique to different values
determine_exclusive_alleles <- function(.data, ...){

  .complete <- .data %>%
    group_by(allele) %>%
    filter(length(unique( !!!rlang::enquos(...) )) == 1) %>%
    ungroup()

  return(.complete)
}

#' Filter alleles shared between multiple combinations of columns
#'
#' @param .data A data frame or tibble
#' @param ... <data-masking> cols to find alleles shared among different values
determine_shared_alleles <- function(.data, ...){

  .complete <- .data %>%
    group_by(allele) %>%
    filter(length(unique( !!!rlang::enquos(...) )) > 1) %>%
    ungroup()

  return(.complete)
}


#' Add `combo` column of combinations of `allele` present in each `biosample`
#' and filter out biosamples with only a single allele
#'
#' @param .data a data frame or tibble
#' @param filter_terms character vector of terms to filter alleles before
#'        determining groups
determine_combinations <- function(.data, filter_terms = NULL){
  .combos <- .data %>%
    {if (!is.null(filter))  filter(., grepl(paste(filter_terms, collapse = "|"), allele, ignore.case = TRUE )) else .} %>%
    group_by(biosample) %>%
    filter(length(unique(allele)) > 1) %>%
    summarize(combo = paste(allele, collapse = ", "), .groups = "drop")

  .complete <- .data %>%
    inner_join(.combos, by = "biosample")

  return(.complete)
}

#' Determine all possible unique proteins present in an MicroBIGG-E dataset
#'
#' Combines a single accession number for each distinct, assigned allele with
#' all distinct accession numbers for unassigned alleles into a single character
#' vector representing the maximum number of potential distinct alleles present
#' in the input dataset
#'
#' @param .data a tibble or dataframe containing MicroBIGG-E data
possible_unique_proteins_mbe <- function(.data){
  accessions_assigned <- .data %>%
    remove_unassigned_bla() %>%
    group_by(allele) %>%
    slice_head() %>%
    ungroup() %>%
    distinct(protein) %>%
    pull()

  accessions_unassigned <- .data %>%
    remove_assigned_bla() %>%
    distinct(protein) %>%
    pull()

  accessions <- c(accessions_assigned, accessions_unassigned)

  return(accessions)
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
#'
determine_nearly_equal_integer_groups_lt <- function(values, groups = 4){
  if(!is.integer(values)){
    stop("Please ensure 'values' is an integer vector.")
  }
  quantiles <- values %>%
    quantile(probs = seq(from = 0, to = 1, length.out = groups + 1), type = 1)

  value_boundaries <- vector(mode = "integer", length = length(quantiles) )
  value_boundaries[1] = quantiles[[1]]

  for(i in 2:(length(quantiles) - 1)){
    last_value <- value_boundaries[[i - 1]]
    threshold = length(values) / groups
    test_value <- quantiles[[i]]

    tibble <- tibble(values) %>%
      filter(values >= last_value)

    include <- tibble %>%
      filter(values <= test_value ) %>%
      nrow()

    exclude <- tibble %>%
      filter(values < test_value ) %>%
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
    stop("Duplicate boundaries. This is likely caused by unbalanced data or too many groups.")
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
#' @param source_vol  integer column to group
#' @param label_col  (new) column for group names
#' @param label_sep seperator between upper and lower bondary in lable
#' @param groups number integer of groups to define (fed into
#'        determine_nearly_equal_integer_groups_lt() to determine boundaries)
#' @param split integer vector of boundaries, which is generally the output of
#'        determine_nearly_equal_integer_groups_lt()

categorize_integer_groups <- function(tibble, source_col, label_col, label_sep = " - ", groups, split = NULL){
  source_col = sym(source_col)
  label_col = sym(label_col)

  if( is.null(split) ){
    bounds <- tibble %>%
      drop_na(!!source_col) %>%
      group_by(biosample) %>%
      filter(row_number() == 1) %>%
      pull(!!source_col) %>%
      determine_nearly_equal_integer_groups_lt(groups)
  }

  lower_bounds <- bounds[1:length(bounds)-1]
  upper_bounds <- bounds[2:length(bounds)]
  bounds_labels <- paste(lower_bounds,label_sep,upper_bounds - 1, sep = "")
  bounds_conditions <- paste(source_col," >= ",lower_bounds," & ",source_col," < ",upper_bounds," ~ ","'",bounds_labels,"'",sep="")

  tibble %>%
    mutate(!!label_col := case_when(!!!rlang::parse_exprs(bounds_conditions))) %>%
    return()
}

########## Other ##########

#' Process na_last arguments
#'
#' Used internally to avoid replicating na_last logic across multiple functions
#'
#' @param .data a dataframe or tibble
#' @param ... <data-masking> columns to reorder NAs
#' @param na_last argument to be processed (pass directly from parent function)
.resolve_na_last <- function(.data, ..., na_last){

  .complete <- .data %>%

    {if(identical(na_last, TRUE)) arrange(., across(c(...), ~ is.na(.)))
      else if(identical(na_last, FALSE)) arrange(., across(c(...), ~ !is.na(.)))
      else if(is.na(na_last)) drop_na(., ...)
      else arrange(., across(c(...), ~ is.na(.))) %>%
        pipe_inform(., "Invalid `na_last` provided. Defaulted to TRUE.") }

  return(.complete)
}

#' Apply separate_rows sequentially by column
#'
#' Useful for separating and all unique combinations of separated rows when
#' different columns may have different numbers of delimiters for the same row
#' (which causes a vector recycling error when applying separate_rows to
#' multiple columns). This essentially "expands" the data column by column.
#'
#' This approach is relatively slow, but serves the purpose without needing to
#' modify separate_rows directly.
#'
#' @param .data a dataframe or tibble
#' @param ... <tidy-select> columns to apply separate_rows to
#' @param sep Separator delimiting collapsed values.
#' @param convert If TRUE will automatically run type.convert() on the key
#'   column. This is useful if the column types are actually numeric, integer,
#'   or logical.
separate_rows_sequential <- function(.data, ..., sep = "[^[:alnum:].]+", convert = FALSE){
  vars <- tidyselect::eval_select(expr(c(...)), .data)
  .modified <- .data

  for(name in names(vars)){
    .modified <- separate_rows(.modified, any_of(name), sep = sep, convert = convert)
  }

  return(.modified)
}

########## File Handling ##########

#' Ensure directory exists (and make if it doesn't)
#'
#' @param path directory path
ensure_directory <- function(path){
  if(!dir.exists(path)){
    dir.create(path, recursive = TRUE)
    message(paste("Directory created:", path))
  }
}
