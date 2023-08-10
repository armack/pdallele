#' @importFrom forcats fct_expand fct_relevel fct_reorder
#' @import dplyr
#' @importFrom rlang :=
NULL

#' Rank, sort, and combine values
#'
#' @description
#' These functions are primarily designed to prepare data for use in
#' [donut_plot()] and serve three main purposes:
#'
#' 1. Determine the `n` largest (or smallest) `count`
#' 2. Combine remaining `count` into `other`
#' 3. Optionally relevel `...` using [relevel_numeric()]
#'
#' They differ in the handling of groups:
#'
#' * `top_n_overall()` determines the "top n" of the entire dataset and keeps those labels across all groups.
#' * `top_n_by_group()` determines the "top n" on a group-by-group basis and keeps the "top n" labels for each group independently.
#' * `top_n_across_groups()` determines the "top n" on a group-by-group basis and keeps labels appearing in the "top n" of any group across all groups.
#'
#' All values that aren't within the 'top n' are combined into `other`, which
#' can be dropped from the output if desired.
#'
#' No tie-breaking is performed. In the event of a tie for the `n`th position,
#' the first row by `arrange(...)` is kept.
#'
#' @section Output Length:
#'
#' Number of rows included in the output tibble willvary based on the function
#' used and the number of groups `group_by(data, ...)` creates.
#'
#' Assuming each group contains at least `n` labels (i.e. `group_by(data, ...)`
#' groups):
#'
#' * `top_n_overall()` outputs a total of `n` distinct labels with all others combined into `other`.
#' * `top_n_by_group()`  outputs a total of `n * n_groups(group_by(data, ...))` distinct labels with all others combined into`other`.
#' * `top_n_across_groups()` outputs a maximum of `n * n_groups(group_by(data, ...))` distinct labels with all others combined into `other`.
#'
#' The exact number of rows output by `top_n_across_groups()` will vary based on
#' the number of shared values of `...` between groups, but will always be
#' between `n` (e.g. complete overlap) and `n * n_groups(group_by(data, ...))`
#' (e.g. no overlap).
#'
#' @param data A data frame or tibble
#' @param ... <data-masking> Columns with 'label' values for `count`
#' @param count Column of values (typically counts) to rank
#' @param n Number of 'labels' to keep per group (all other values will be
#'   collapsed into `other`)
#' @param desc Should the n highest (`TRUE`) or lowest (`FALSE`) values of
#'   `count` be kept?
#' @param other Label for the `other` group
#' @param other_pos Where should `other` sort? Options: `first` before ranked
#'   values, `last` after ranked values, `drop` drop all non-ranked values.
#' @param relevel should the output be releveld based on ... (`name`),
#'   values (`count`), or not at all (`none`)?
#'
#' @export
top_n_overall <- function(data, ..., count, n = 10, desc = TRUE,
                          other = "Other", relevel = "name",
                          other_pos = "last"){
  . <- NULL # Workaround to suppress `no visible binding for global variable`
  dots <- rlang::enquos(...)

  if(...length() < 1){
    stop("Must provide column names in `...`")
  }else if(!all(unlist(purrr::map(dots, rlang::as_name)) %in% names(data))){
    stop("Please ensure values given in `...` are valid column names.")
  }

  if(!(relevel %in% c("name","count","none"))){
    rlang::inform("Invalid `relevel` provided. No releveling performed.")
  }

  if(!(other_pos %in% c("drop","first","last","none"))){
    rlang::inform("Invalid `other_pos` provided. No releveling performed.")
  }

  .top_values <- data %>%
    group_by(...) %>%
    summarize(total = sum({{count}}), .groups = "drop") %>%
    {if(desc) arrange(., desc("total")) else arrange(., "total")} %>%
    filter(row_number() <= n) %>%
    mutate(concat = paste0(!!!dots)) %>%
    pull("concat")

  .complete <- data %>%
    mutate(concat = paste0(!!!dots)) %>%
    mutate(group_total = if_else(row_number() == 1, sum({{count}}), NA_integer_)) %>%
    mutate(across(c(...), ~fct_expand(., cur_column(), other))) %>%
    mutate(across(c(...), ~if_else(concat %in% .top_values, ., forcats::as_factor(other)))) %>%
    group_by(..., .data$group_total, .add = TRUE) %>%
    summarize({{count}} := sum({{count}}), .groups = "drop") %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    {if (identical(relevel,"name")) relevel_numeric(., ...)
      else if(identical(relevel,"count")) mutate(., across(c(...), forcats::fct_reorder(., {{count}})))
      else .} %>%
    {if (identical(other_pos, "last")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = Inf))))
      else if(identical(other_pos, "first")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = 0L))))
      else if(identical(other_pos, "drop")) filter(., !grepl(other, !!!dots))
      else .} %>%
    arrange(...)

  return(.complete)
}


#' Keep the `n` highest (or lowest) `count` rows in each group across all
#' groups and combine the rest into `other`
#'
#' @rdname top_n_overall
#' @export
top_n_across_groups <- function(data, ..., count, n = 10,
                                    desc = TRUE, other = "Other",
                                    relevel = "name", other_pos = "last"){
  . <- NULL # Workaround to suppress `no visible binding for global variable`
  dots <- rlang::enquos(...)

  if(...length() < 1){
    stop("Must provide column names in `...`")
  }else if(!all(unlist(purrr::map(dots, rlang::as_name)) %in% names(data))){
    stop("Please ensure values given in `...` are valid column names.")
  }

  if(!(relevel %in% c("name","count","none"))){
    rlang::inform("Invalid `relevel` provided. No releveling performed.")
  }

  if(!(other_pos %in% c("drop","first","last","none"))){
    rlang::inform("Invalid `other_pos` provided. No releveling performed.")
  }

  .complete <- data %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    mutate(top = if_else(row_number() <= n, paste(!!!dots), NA_character_)) %>%
    mutate(group_total = if_else(row_number() == 1, sum({{count}}), NA_integer_)) %>%
    with_groups(NULL, mutate, keep = paste(!!!dots) %in% .data$top) %>%
    mutate(across(c(...), ~fct_expand(., cur_column(), other))) %>%
    mutate(across(c(...), ~if_else(keep, ., forcats::as_factor(other)))) %>%
    group_by(..., .data$group_total, .add = TRUE) %>%
    summarize({{count}} := sum({{count}}), .groups = "drop") %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    {if (identical(relevel,"name")) relevel_numeric(., ...)
      else if(identical(relevel,"count")) mutate(., across(c(...), ~forcats::fct_reorder(., {{count}})))
      else .} %>%
    {if (identical(other_pos, "last")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = Inf))))
      else if(identical(other_pos, "first")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = 0L))))
      else if(identical(other_pos, "drop")) filter(., !grepl(other, !!!dots))
      else .} %>%
    arrange(...)

  return(.complete)
}

#' Keep the `n` highest (or lowest) `count` rows by group and combine the rest
#' into `other`
#'
#' @rdname top_n_overall
#' @export
top_n_by_group <- function(data, ..., count, n = 10, desc = TRUE,
                           other = "Other", relevel = "name",
                           other_pos = "last"){
  . <- NULL # Workaround to suppress `no visible binding for global variable`
  dots <- rlang::enquos(...)

  if(...length() < 1){
    stop("Must provide column names in `...`")
  }else if(!all(unlist(purrr::map(dots, rlang::as_name)) %in% names(data))){
    stop("Please ensure values given in `...` are valid column names.")
  }

  if(!(relevel %in% c("name","count","none"))){
    rlang::inform("Invalid `relevel` provided. No releveling performed.")
  }

  if(!(other_pos %in% c("drop","first","last","none"))){
    rlang::inform("Invalid `other_pos` provided. No releveling performed.")
  }

  .complete <- data %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    mutate(group_total = if_else(row_number() == 1, sum({{count}}), NA_integer_)) %>%
    mutate(across(c(...), ~fct_expand(., cur_column(), other))) %>%
    mutate(across(c(...), ~if_else(row_number() <= n, ., forcats::as_factor(other)))) %>%
    group_by(..., .data$group_total, .add = TRUE) %>%
    summarize({{count}} := sum({{count}}), .groups = "drop") %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    {if (identical(relevel,"name")) relevel_numeric(., ...)
      else if(identical(relevel,"count")) mutate(., across(c(...), ~forcats::fct_reorder(., {{count}})))
      else .} %>%
    {if (identical(other_pos, "last")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = Inf))))
      else if(identical(other_pos, "first")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = 0L))))
      else if(identical(other_pos, "drop")) filter(., !grepl(other, !!!dots))
      else .} %>%
    arrange(...)

  return(.complete)
}

