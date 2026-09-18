# Shared formatting for true log10 axes in current primary reporting.
# Minor lines mark linear multiples within each decade: 2, 3, ..., 9
# times each power of ten (e.g. 20, 30, ..., 90 between 10 and 100).
# Supply original-scale limits, as ggplot2 does for minor_breaks functions.
cbh_log10_minor_breaks <- function(limits) {
  limits <- limits[is.finite(limits) & limits > 0]
  if (length(limits) < 2L) return(numeric())
  decades <- seq.int(floor(log10(min(limits))), floor(log10(max(limits))))
  breaks <- sort(as.vector(outer(2:9, 10^decades, `*`)))
  breaks[breaks >= min(limits) & breaks <= max(limits)]
}

cbh_log10_grid_theme <- function() {
  ggplot2::theme(panel.grid.minor = ggplot2::element_line(
    colour = "grey88", linewidth = 0.25))
}
