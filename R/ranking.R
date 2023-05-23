#' @importFrom forcats fct_expand fct_relevel fct_reorder
#' @import dplyr
#' @importFrom rlang :=
NULL
#' Keep the `n` highest (or lowest) `count` rows in each group across all
#' groups and combine the rest into `other`
#'
#' This function keeps the "top n" in each group across all groups. See the
#' companion function top_n_by_group() to keep the "top n" in each group.
#'
#' Note that no tie-breaking is performed. In the event of a tie for the nth
#' row, the first row by arrange() is kept.
#'
#' @param data a dataframe or tibble
#' @param ... <data-masking> columns of labels to be ranked
#' @param count column of values (typically counts) to rank
#' @param n number of rows to keep per group (all other values will be collapsed
#' into an `other` row, meaning n+1 total rows will appear in out unless
#' `other_pos` is set to "drop"
#' @param desc should the n highest values of `count` be kept?
#' @param other character label for the "other" group
#' @param other_pos where should "other" sort?
#' @param relevel should factors in `...` be releveld in alphanumeric order?
#' @export
new_top_n_across_groups <- function(data, ..., count, n = 10,
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
    mutate(across(c(...), ~if_else(keep, ., as_factor(other)))) %>%
    group_by(..., "group_total", .add = TRUE) %>%
    summarize({{count}} := sum({{count}}), .groups = "drop") %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    {if (identical(relevel,"name")) relevel_numeric(., ...)
      else if(identical(relevel,"count")) mutate(., across(c(...), ~fct_reorder(., {{count}})))
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
#' This function keeps the "top n" with each group. See the companion function
#' top_n_across_groups() to keep the "top n" in each group across all groups.
#'
#' Note that no tie-breaking is performed. In the event of a tie for the nth
#' row, the first row by arrange() is kept.
#'
#' @param data a dataframe or tibble
#' @param ... <data-masking> columns of labels to be ranked
#' @param count column of values (typically counts) to rank
#' @param n number of rows to keep per group (all other values will be collapsed
#' into an `other` row, meaning n+1 total rows will appear in out unless
#' `other_pos` is set to "drop"
#' @param desc should the n highest values of `count` be kept?
#' @param other character label for the "other" group
#' @param other_pos where should "other" sort?
#' @param relevel should factors in `...` be releveled in alphanumeric order?
#'
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
    mutate(across(c(...), ~if_else(row_number() <= n, ., as_factor(other)))) %>%
    group_by(..., "group_total", .add = TRUE) %>%
    summarize({{count}} := sum({{count}}), .groups = "drop") %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    {if (identical(relevel,"name")) relevel_numeric(., ...)
      else if(identical(relevel,"count")) mutate(., across(c(...), ~fct_reorder(., {{count}})))
      else .} %>%
    {if (identical(other_pos, "last")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = Inf))))
      else if(identical(other_pos, "first")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = 0L))))
      else if(identical(other_pos, "drop")) filter(., !grepl(other, !!!dots))
      else .} %>%
    arrange(...)

  return(.complete)
}

#' Keep the `n` highest (or lowest) `count` rows overall for each group and
#' combine the rest into `other`
#'
#' This function keeps the "top n" overall for each group. See the companion function
#' top_n_across_groups() to keep the "top n" in each group across all groups.
#'
#' Note that no tie-breaking is performed. In the event of a tie for the nth
#' row, the first row by arrange() is kept.
#'
#' @param data a dataframe or tibble
#' @param ... <data-masking> columns of labels to be ranked
#' @param count column of values (typically counts) to rank
#' @param n number of rows to keep per group (all other values will be collapsed
#' into an `other` row, meaning n+1 total rows will appear in out unless
#' `other_pos` is set to "drop"
#' @param desc should the n highest values of `count` be kept?
#' @param other character label for the "other" group
#' @param other_pos where should "other" sort?
#' @param relevel should factors in `...` be releveled in alphanumeric order?
#'
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
    mutate(across(c(...), ~if_else(concat %in% .top_values, ., as_factor(other)))) %>%
    group_by(..., "group_total", .add = TRUE) %>%
    summarize({{count}} := sum({{count}}), .groups = "drop") %>%
    {if(desc) arrange(., desc({{count}})) else arrange(., {{count}})} %>%
    {if (identical(relevel,"name")) relevel_numeric(., ...)
      else if(identical(relevel,"count")) mutate(., across(c(...), fct_reorder(., {{count}})))
      else .} %>%
    {if (identical(other_pos, "last")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = Inf))))
      else if(identical(other_pos, "first")) mutate(., across(c(...), ~ suppressWarnings(fct_relevel(., other, after = 0L))))
      else if(identical(other_pos, "drop")) filter(., !grepl(other, !!!dots))
      else .} %>%
    arrange(...)

  return(.complete)
}
