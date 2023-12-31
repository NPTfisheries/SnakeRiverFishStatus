# Purpose: map theme

theme_map <- function(...) {
  theme_void() +
    theme(
      text = element_text(family = "serif", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      #panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
      #legend.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      ...
    )
}
