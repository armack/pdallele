########## Bar Plot Functions ##########

########## Donut Plot Functions ##########

#' Generates a panel of donut plots of `fill` and `counts`
#'
#' See companion function donut_plot_panel() for a variation better suited for
#' comparison between groups of data
#'
#' @param .data a data table or tibble
#' @param fill column containing levels/labels (e.g. alleles names) (rich text)
#' @param counts column containing counts for the levels
#' @param count_label label for counts (e.g. Alleles) (rich text)
#' @param legend_title label for the legend (rich text)
#' @param labels should callouts with count & percentage be shown for each arc?
#' @param total should total of `count` be displayed inside the donut?
#' @param palette specify a ggthemes::tableau_color_pal() color palette name


donut_plot <- function(.data, fill = allele_markdown, counts = alleles,
                       count_label = "Alleles", legend_title = "**Allele**",
                       labels = TRUE, total = TRUE, palette = NULL){
  .donut_data <- .data %>%
    mutate(y_proportion = {{counts}} / sum({{counts}})) %>%
    mutate(y_position = 1 - cumsum(y_proportion) + 0.5 * y_proportion) %>%
    mutate(y_label =
             paste0({{counts}}, " ", count_label, "\n (",
                    format(round(100 * {{counts}} / sum({{counts}}), digits = 1), nsmall = 1),
                    "%)")) %>%
    mutate(group_label = if_else(is.na(group_total), NA_character_,
                                 paste0(group_total,"<br />", count_label)))

  plot <- ggplot(.donut_data, aes(y = y_proportion)) +
    theme_donut() +
    geom_col(aes(fill = {{fill}}, x = 2)) +
    coord_polar(theta = "y") +
    xlim(0, 4) +
    labs(fill = legend_title)

  # Add labels for each arc
  if (labels){
    plot <- plot + ggrepel::geom_label_repel(
      aes(x = 2.5, y = y_position, label = y_label), min.segment.length = 0,
      nudge_x = 1.25, force = 2, force_pull = 0, max.overlaps = 200
    )
  }

  # Add center total label
  if(total){
    plot <- plot +
      ggtext::geom_richtext(aes(label = group_label, x = 0, y = 0),
                            fill = NA, label.colour = NA, na.rm = TRUE)
  }

  # Number of arcs for choosing colors and colors
  arcs <- length(unique(pull(.donut_data, {{fill}})))

  # Default to Tableau color palettes if not too long. Allow selection too.
  if (!is_null(palette))
    # No error checking - allow scale_fill_tableau to handle
    plot <- plot +
    ggthemes::scale_fill_tableau(palette = palette, type = palette_type) +
    ggthemes::scale_color_tableau(palette = palette, type = palette_type)
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

#' Generates a panel of donut plots of `fill` and `counts` faceted by `facet`
#'
#' See companion function donut_plot() for a variation better suited for
#' a single condition without faceting
#'
#' @param .data a data table or tibble
#' @param fill column containing levels/labels (e.g. alleles names) (rich text)
#' @param counts column containing counts for the levels
#' @param facet column to facet plots by (e.g. region)
#' @param count_label label for counts (e.g. Alleles) (rich text)
#' @param legend_title label for the legend (rich text)
#' @param total should total of `count` be displayed inside the donut?
#' @param palette specify a ggthemes::tableau_color_pal() color palette name
#' @param nrow number of facet rows
#' @param ncol number of facet columns

donut_plot_panel <- function(.data, fill = allele_markdown, counts = alleles, facet,
                             count_label = "Alleles", legend_title = "**Allele**",
                             total = TRUE, palette = NULL, nrow = NULL, ncol = NULL){
  .donut_data <- .data %>%
    group_by( {{facet}} ) %>%
    mutate(y_proportion = alleles / sum(alleles)) %>%
    mutate(group_label = if_else(is.na(group_total), NA_character_,
                                 paste0(group_total,"<br />", count_label)))

  plot <- ggplot(.donut_data, aes(y = y_proportion)) +
    theme_donut() +
    geom_col(aes(fill = {{fill}}, x = 2)) +
    coord_polar(theta = "y") +
    xlim(0, 2.5) +
    labs(fill = legend_title) +
    facet_wrap(facets = vars( {{facet}} ), nrow = nrow, ncol = ncol)

  if(total){
    plot <- plot +
      ggtext::geom_richtext(aes(label = group_label, x = 0, y = 0),
                            fill = NA, label.colour = NA, na.rm = TRUE)
  }

  # Number of arcs for choosing colors and colors
  arcs <- length(unique(pull(.donut_data, {{fill}})))

  # Default to Tableau color palettes if not too long. Allow selection too.
  if (!is_null(palette))
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

########## Heatmap Functions ##########

#' Prepare data for rank and proprtion heatmap
#'
#' Count all instances of `count` by `group`, filter `n` most frequent `count`
#' either overall or `by_group`. Optionally prepare to facet by `facet` and
#' remove NA values.
#'
#' Proportion refers to the proportion of isolates in a given `group` associated
#' with each `count`. Ranks are the `dense_rank(desc())` of the number of
#' isolates associated with each `count`.
#'
#' Note: The default setting `by_group = TRUE` will keep all `count` that appear
#' in the `n` most frequent in any `group`, typically leading to more than `n`
#' total values being kept. Conversely, `by_group = FALSE` will keep only the
#' `n` most frequent from the entire tibble, ensuring only `n` values in the
#' output.
#'
#' @param .data a dataframe or tibble
#' @param count <data-masking> column to count (long axis)
#' @param group <data-masking> column to compare across (short axis)
#' @param n number of most frequent `count` to include (i.e. top n =)
#' @param facet <data-masking>column to facet by
#' @param by_group should `n` be aggregated by group or for the overall dataset?
#' @param na.rm should NA values be removed from `count`, `group`, and `facet`?

heatmap_data <- function(.data, count = allele_markdown, group, n = 5, facet,
                         by_group = TRUE, na.rm = TRUE) {

  facet_var <- rlang::enquo(facet)
  has_facets <- !rlang::quo_is_missing(facet_var)

  .complete <- .data %>%
    ungroup() %>%
    { if(has_facets) group_by(., {{facet}}) else . } %>%
    {if(na.rm) drop_na(., c({{count}}, {{group}},
                            if(has_facets) any_of(rlang::as_name(facet_var)) else NULL)) else . } %>%
    mutate(total_isolates = length(unique(biosample))) %>%
    group_by({{group}}, .add = TRUE) %>%
    mutate(group_isolates = length(unique(biosample))) %>%
    select(c({{count}}, {{group}}, group_isolates, total_isolates,
             if(has_facets) any_of(rlang::as_name(facet_var)) else NULL)) %>%
    group_by({{count}}, .add = TRUE) %>%
    mutate(group_count = n(), group_prop = group_count/group_isolates) %>%
    ungroup({{group}}) %>%
    mutate(total_count = n(), total_prop = total_count/total_isolates) %>%
    group_by({{group}}, .add = TRUE) %>%
    slice_head() %>%
    ungroup({{count}}) %>%
    mutate(group_rank = dense_rank((desc(group_count)))) %>%
    ungroup({{group}}) %>%
    mutate(total_rank = dense_rank((desc(total_count)))) %>%
    {if(identical(by_group, TRUE)) filter(., {{count}} %in% pull(filter(., group_rank <= n), {{count}}))
      else filter(., {{count}} %in% pull(filter(., total_rank <= n), {{count}})) } %>%
    group_by({{count}}, .add = TRUE) %>%
    mutate(group_sort = length(unique({{group}})) / mean(group_rank)) %>%
    mutate(total_sort =  mean(total_rank)) %>%
    ungroup() %>%
    {if(identical(by_group, TRUE)) mutate(., {{count}} := fct_reorder({{count}}, group_sort, .desc = F))
      else mutate(., {{count}} := fct_reorder({{count}}, total_sort, .desc = T)) }

  if(has_facets){
    levels <- .complete %>%
      select({{count}}, {{facet}}) %>%
      relevel_numeric({{facet}}) %>%
      arrange(desc({{facet}}), {{count}}) %>%
      distinct({{count}}) %>%
      pull({{count}}) %>%
      as.character()

    .complete <- .complete %>%
      mutate({{count}} := fct_relevel({{count}}, levels))
  }

  return(.complete)
}

#' Prepare ggplot object of heatmap
#'
#' Notes on Plot Orientation:
#'
#' `wide = FALSE` places `count` on the Y-axis and `group` on the X-axis which
#' will result in "longer" plots when `count` has more levels than `group`.
#' Plots will facet vertically.
#'
#' `wide = TRUE` swaps the axes and will result in "wider" plots when `count`
#' has more levels than `group`. Plots will facet horizontally.
#'
#' Facet directions are fixed with plot orientation to avoid sections or to
#' avoid unexpected, unintended, and undesirable "tiling" of plots
#'
#' @param .data a dataframe or tibble
#' @param count <data-masking> column to count (long axis)
#' @param group <data-masking> column to compare across (short axis)
#' @param wide should the output be wide or tall? See description.
#' @param rank title for the rank (size) legend
#' @param prop title for the proportion (color) legend
#' @param prop_limits Limits of the proportion color scale. Passed to
#'   `ggplot2::continuous_scale(limits = )`
#' @param prop_breaks Proportion values to be displayed in the legend. Passed to
#'   `ggplot2::continuous_scale(breaks = )`
#' @param prop_labels Labes corresponding to `prop_breaks`. Passed to
#'   `ggplot2::continuous_scale(labels = )`
#' @param prop_barheight Height of the bar for the proportion legend Passed to
#'   `ggplot2::guide_colorbar(barheight = )`
#' @param rank_breaks Rank values to be displayed in the legend. Passed to
#'   `ggplot2::scale_size(breaks = )`.
#' @param rank_limits Limits of the rank size scale. Passed to
#'   `ggplot2::scale_size(limits = )`
#' @param rank_range Minimum and maximum size of the circle. Passed to
#'   `ggplot2::scale_size(range = )`
#' @param text_size Text size in plot. Passed to `ggplot2::element_text()`.
#' @param option Color palette, passed to `viridis::scale_color_viridis()`.

heatmap_plot <- function(.data, count = allele_markdown, group, facet,
                         option = "viridis", wide = FALSE,
                         rank = "**Rank**", prop = "**Proportion**",
                         prop_breaks = c(0.0001,0.0005,0.001,0.0025,
                                         0.005,0.01,0.025,0.05,0.1,
                                         0.2,0.3,0.5,0.75,.90),
                         prob_labels = c("0.01%","0.05%","0.1%","0.25%",
                                         "0.5%","1%","2.5%","5%","10%",
                                         "20%","30%","50%","75%","90%"),
                         prop_barheight = 15, prop_limits = NULL,
                         rank_breaks = c(1,5,10,25,50,75,100,150,200),
                         rank_limits = NULL, rank_range = c(12,1),
                         text_size = 16){

  has_facets <- !rlang::quo_is_missing(rlang::enquo(facet))

  if(wide)
    .data <- mutate(.data, {{count}} := fct_rev({{count}}))

  if(!wide)
    .start <- ggplot(.data, aes(x = {{group}}, y = {{count}}))
  else
    .start <- ggplot(.data, aes(x = {{count}}, y = {{group}}))

  .plot <- .start  +
    theme_heatmap() %+replace%
    theme(text = element_text(size = text_size)) +
    geom_point(aes(color = group_prop, size = group_rank), show.legend = TRUE) +
    viridis::scale_color_viridis(option = option, trans = "log", breaks = prop_breaks,
                                 labels = prob_labels, limits = prop_limits) +
    scale_size(range = rank_range, limits = rank_limits, breaks = rank_breaks) +
    labs(size = rank, color = prop) +
    guides(color = guide_colorbar(barheight = prop_barheight, order = 2),
           size = guide_legend(order = 1))

  if(has_facets){
    if(!wide)
      .plot <- .plot +
        facet_grid(rows = vars({{facet}}), scales = "free_y", space = "free_y")
    else
      .plot <- .plot +
        facet_grid(cols = vars({{facet}}), scales = "free_x", space = "free_x")
  }

  return(.plot)
}
