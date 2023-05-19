theme_simple <- function(base_size = 11, base_family = "",
                     base_line_size = base_size / 22,
                     base_rect_size = base_size / 22) {
  # Starts with theme_void and then modify some parts
  theme_minimal(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    theme(
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
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      
      complete = TRUE
    )
}

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
    theme(
      axis.title = element_blank(),
      axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5),
      axis.text.y = ggtext::element_markdown(),
      axis.ticks = element_blank(),
      
      panel.background = element_blank(),
      panel.spacing = unit(0.5, "in"),
      strip.background = element_blank(),
      strip.text = element_blank(),
  
      complete = TRUE
    )
}