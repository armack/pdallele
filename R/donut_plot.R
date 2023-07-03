#' @importFrom dplyr %>%
NULL

########## Donut Plot Functions ##########

#' Generate donut plots
#'
#' @description Generates a donut plot using `fill` as categories and `count` as
#'   counts.
#'
#'   * `donut_plot()` generates a single donut plot
#'   * `donut_plot_panel()` generates a panel of donut plots faceted by `facet`
#'
#'   Data can be prepared for plotting using a ranking function
#'   [top_n_overall()], [top_n_across_groups()], or [top_n_by_group()].
#'
#' @param data A data table or tibble
#' @param fill Column containing levels/labels (e.g. alleles names) (rich text)
#' @param counts Column containing counts for the levels
#' @param count_label Label for counts (rich text)
#' @param legend_title Label for the legend (rich text)
#' @param labels Should callouts with count & percentage be shown for each arc?
#' @param total Should total of `count` be displayed inside donuts?
#' @param palette A [ggthemes::tableau_color_pal()] color palette name. By
#'   default, uses '`Tableau 10`' for ≤ 10 categories, '`Tableau 20`' for 11 to
#'   20 categories, and the standard `ggplot2` color palette for > 20 categories
#'
#' @returns A `ggplot2` plot object


donut_plot <- function(data, fill, counts,
                       count_label = "Alleles", legend_title = "**Allele**",
                       labels = TRUE, total = TRUE, palette = NULL){
  . <- NULL # Workaround to suppress `no visible binding for global variable`

  .donut_data <- data %>%
    dplyr::mutate(y_proportion = {{counts}} / sum({{counts}})) %>%
    dplyr::mutate(y_position = 1 - cumsum(.data$y_proportion) + 0.5 * .data$y_proportion) %>%
    dplyr::mutate(y_label =
                    paste0({{counts}}, " ", count_label, "\n (",
                           format(round(100 * {{counts}} / sum({{counts}}), digits = 1), nsmall = 1),
                           "%)")) %>%
    dplyr::mutate(group_label = dplyr::if_else(is.na(.data$group_total), NA_character_,
                                               paste0(.data$group_total,"<br />", count_label)))

  plot <- ggplot2::ggplot(.donut_data, ggplot2::aes(y = .data$y_proportion)) +
    theme_donut() +
    ggplot2::geom_col(ggplot2::aes(fill = {{fill}}, x = 2)) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(0, 4) +
    ggplot2::labs(fill = legend_title)

  # Add labels for each arc
  if (labels){
    plot <- plot + ggrepel::geom_label_repel(
      ggplot2::aes(x = 2.5, y = .data$y_position, label = .data$y_label), min.segment.length = 0,
      nudge_x = 1.25, force = 2, force_pull = 0, max.overlaps = 200
    )
  }

  # Add center total label
  if(total){
    plot <- plot +
      ggtext::geom_richtext(ggplot2::aes(label = .data$group_label, x = 0, y = 0),
                            fill = NA, label.colour = NA, na.rm = TRUE)
  }

  # Number of arcs for choosing colors and colors
  arcs <- length(unique(pull(.donut_data, {{fill}})))

  # Default to Tableau color palettes if not too long. Allow selection too.
  if (!rlang::is_null(palette))
    # No error checking - allow scale_fill_tableau to handle
    plot <- plot +
    ggthemes::scale_fill_tableau(palette = palette) +
    ggthemes::scale_color_tableau(palette = palette)
  else{
    if (arcs <= 10)
      plot <- plot +
        ggthemes::scale_fill_tableau("Tableau 10") +
        ggthemes::scale_color_tableau("Tableau 10")
    else if (arcs > 10 & arcs <= 20)
      plot <- plot +
        ggthemes::scale_fill_tableau("Tableau 20") +
        ggthemes::scale_color_tableau("Tableau 20")
    # don't set a palette for > 20 arcs, just use default ggplot rainbow palette
  }

  return(plot)
}

#' Generate a faceted panel of donut plots
#'
#' @rdname donut_plot
#'
#' @param facet Column to facet plots by
#' @param nrow Number of facet rows
#' @param ncol Number of facet columns

donut_plot_panel <- function(data, fill, counts, facet,
                             count_label = "Alleles", legend_title = "**Allele**",
                             total = TRUE, palette = NULL, nrow = NULL, ncol = NULL){
  .donut_data <- data %>%
    relevel_numeric({{facet}}) %>%
    dplyr::group_by( {{facet}} ) %>%
    dplyr::mutate(y_proportion = {{counts}} / sum({{counts}})) %>%
    dplyr::mutate(group_label = dplyr::if_else(is.na(.data$group_total), NA_character_,
                                               paste0(.data$group_total,"<br />", count_label)))

  plot <- ggplot2::ggplot(.donut_data, ggplot2::aes(y = .data$y_proportion)) +
    theme_donut() +
    ggplot2::geom_col(ggplot2::aes(fill = {{fill}}, x = 2)) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(0, 2.5) +
    ggplot2::labs(fill = legend_title) +
    ggplot2::facet_wrap(facets = vars( {{facet}} ), nrow = nrow, ncol = ncol)

  if(total){
    plot <- plot +
      ggtext::geom_richtext(ggplot2::aes(label = .data$group_label, x = 0, y = 0),
                            fill = NA, label.colour = NA, na.rm = TRUE)
  }

  # Number of arcs for choosing colors and colors
  arcs <- length(unique(pull(.donut_data, {{fill}})))

  # Default to Tableau color palettes if not too long. Allow selection too.
  if (!rlang::is_null(palette))
    # No error checking - allow scale_fill_tableau to handle
    plot <- plot +
    ggthemes::scale_fill_tableau(palette = palette) +
    ggthemes::scale_color_tableau(palette = palette)
  else{
    if (arcs <= 10)
      plot <- plot +
        ggthemes::scale_fill_tableau("Tableau 10") +
        ggthemes::scale_color_tableau("Tableau 10")
    else if (arcs > 10 & arcs <= 20)
      plot <- plot +
        ggthemes::scale_fill_tableau("Tableau 20") +
        ggthemes::scale_color_tableau("Tableau 20")
    # don't set a palette for > 20 arcs, just use default ggplot rainbow palette
  }

  return(plot)
}
