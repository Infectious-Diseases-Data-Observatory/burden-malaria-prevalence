# Shared formatting for true log10 axes in current primary reporting.
# Minor lines are 0.1 log10 units apart (nine inside each decade).
# Supply original-scale limits, as ggplot2 does for minor_breaks functions.
cbh_log10_minor_breaks <- function(limits) {
  limits <- limits[is.finite(limits) & limits > 0]
  if (length(limits) < 2L) return(numeric())
  first <- ceiling(10 * log10(min(limits)) - 1e-10)
  last <- floor(10 * log10(max(limits)) + 1e-10)
  if (first > last) return(numeric())
  10^(seq.int(first, last) / 10)
}

cbh_log10_grid_theme <- function() {
  ggplot2::theme(panel.grid.minor = ggplot2::element_line(
    colour = "grey88", linewidth = 0.25))
}
