#' @importFrom ggplot2 %+replace%
NULL

#' Base theme for `pdallele` plots
#' @inheritParams ggplot2::theme_gray
theme_simple <- function(base_size = 11, base_family = "",
                     base_line_size = base_size / 22,
                     base_rect_size = base_size / 22) {
  # Starts with theme_void and then modify some parts
  ggplot2::theme_minimal(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    ggplot2::theme(
      # Add ggtext::element_markdown() to in-use text
      legend.text      = ggtext::element_markdown(),
      legend.title     = ggtext::element_markdown(),
      axis.text.x      = ggtext::element_markdown(),
      axis.text.y      = ggtext::element_markdown(),

      # Bold some useful text
      strip.text       = ggtext::element_markdown(face = "bold"),

      complete = TRUE
    )
}

#' Additional theme elements for `pdallele` donut plots
#' @inheritParams theme_simple
theme_donut <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  # Starts with theme_simple and then modify some parts
  theme_simple(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),

      complete = TRUE
    )
}

#' Additional theme elements for `pdallele` heatmaps
#' @inheritParams theme_simple
theme_heatmap <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  # Starts with theme_simple and then modify some parts
  theme_simple(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5),
      axis.text.y = ggtext::element_markdown(),
      axis.ticks = ggplot2::element_blank(),

      panel.background = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.5, "in"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),

      complete = TRUE
    )
}
