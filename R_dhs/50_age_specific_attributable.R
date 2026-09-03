# =============================================================================
# 50_age_specific_attributable.R — all-cause mortality attributable to malaria
# by age band as a function of PfPR2-10, from the six age-band person-time
# models (script 42, smooth dose-response, 4-month split).
#
# Two quantities, each with a 95% interval from 1,000 draws of the fitted
# coefficients:
#   attributable fraction  AF_a(p) = 1 - rate_a(1%) / rate_a(p): the share of
#       all-cause deaths in band a that the model attributes to prevalence p
#       rather than the 1% reference, holding covariates, calendar year and the
#       random effects at their reference values.
#   attributable rate      rate_a(p) - rate_a(1%), in deaths per 1,000
#       child-years, where rate_a(p) is the model's predicted all-cause death
#       rate in band a at prevalence p for an average cell (covariates at their
#       means, year centred, first window, random effects at zero). Segments
#       within a band are weighted by their observed person-time.
#
# Reads the model bundle script 42 writes; run script 42 first.
#
# Outputs (results/dhs_rebuild)
#   person_time_age_specific_attributable.csv   band x PfPR grid: AF and attributable rate
#   person_time_age_specific_anchors.csv        the same at 10, 30, 50% prevalence
#   figure36_age_specific_attributable.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "MASS", "ggplot2", "patchwork"))

N_DRAWS <- 1000L
GRID <- seq(0, 80, by = 1)
ANCHORS <- c(10, 30, 50)
set.seed(20260903)

bundle <- readRDS(file.path(DERIVED_DIR, "person_time_model_bundle.rds"))
fits <- bundle$bands6b_smooth_fits
model_data <- read.csv(file.path(RESULTS_DIR, "person_time_model_data.csv"), stringsAsFactors = FALSE)
country_levels <- levels(factor(model_data$iso3))
survey_levels <- levels(factor(model_data$svkey))

band_curves <- function(band) {
  fit <- fits[[band]]
  band_data <- model_data[model_data$age6b == band, ]
  segment_levels <- levels(droplevels(factor(band_data$segment)))
  # person-time share of each segment within the band, to average segment rates
  segment_share <- tapply(band_data$person_months_n, factor(band_data$segment, levels = segment_levels), sum)
  segment_share <- segment_share / sum(segment_share)
  frame <- function(p, segment) {
    out <- data.frame(pfpr10 = p / 10,
                      segment_f = factor(segment, levels = segment_levels),
                      window_f = factor("1", levels = as.character(1:5)), year_c = 0, log_pm = 0,
                      country = factor(country_levels[1], levels = country_levels),
                      survey = factor(survey_levels[1], levels = survey_levels))
    G <- matrix(0, nrow(out), length(grep("^G", names(coef(fit)))))
    colnames(G) <- sub("^G", "", grep("^G", names(coef(fit)), value = TRUE))
    out$G <- G
    out
  }
  draws <- MASS::mvrnorm(N_DRAWS, coef(fit), vcov(fit))
  zero_re <- function(X) { X[, grep("^s\\(country\\)|^s\\(survey\\)", colnames(X))] <- 0; X }
  # deaths per person-month for each segment on the grid and at the reference
  rate_grid <- 0; rate_ref <- 0
  for (s in seq_along(segment_levels)) {
    Xg <- zero_re(predict(fit, frame(GRID, segment_levels[s]), type = "lpmatrix", discrete = FALSE))
    Xr <- zero_re(predict(fit, frame(AF_REFERENCE, segment_levels[s]), type = "lpmatrix", discrete = FALSE))
    rate_grid <- rate_grid + segment_share[[s]] * exp(draws %*% t(Xg))
    rate_ref <- rate_ref + segment_share[[s]] * exp(draws %*% t(Xr))
  }
  rate_ref <- as.numeric(rate_ref)
  af <- 1 - rate_ref / rate_grid
  attributable <- (rate_grid - rate_ref) * 12000   # per 1,000 child-years
  rate <- rate_grid * 12000
  q <- function(m, f) apply(m, 2, f)
  data.frame(age_band = band, pfpr = GRID,
             af = colMeans(af), af_lo = q(af, function(v) stats::quantile(v, 0.025)),
             af_hi = q(af, function(v) stats::quantile(v, 0.975)),
             attributable_per_1000cy = colMeans(attributable),
             attributable_lo = q(attributable, function(v) stats::quantile(v, 0.025)),
             attributable_hi = q(attributable, function(v) stats::quantile(v, 0.975)),
             allcause_rate_per_1000cy = colMeans(rate),
             rate_lo = q(rate, function(v) stats::quantile(v, 0.025)),
             rate_hi = q(rate, function(v) stats::quantile(v, 0.975)),
             reference_rate_per_1000cy = mean(rate_ref) * 12000,
             stringsAsFactors = FALSE)
}
curves <- do.call(rbind, lapply(AGE6B, band_curves))
rownames(curves) <- NULL
curves$age_band <- factor(curves$age_band, levels = AGE6B)
# the data support: the prevalence range actually observed in each band
support <- aggregate(pfpr ~ age6b, data = model_data, FUN = function(v) stats::quantile(v, c(0.01, 0.99)))
support <- data.frame(age_band = factor(support$age6b, levels = AGE6B), p01 = support$pfpr[, 1], p99 = support$pfpr[, 2])
write.csv(curves, file.path(RESULTS_DIR, "person_time_age_specific_attributable.csv"), row.names = FALSE)

anchors <- curves[curves$pfpr %in% ANCHORS, ]
write.csv(anchors, file.path(RESULTS_DIR, "person_time_age_specific_anchors.csv"), row.names = FALSE)
message("Attributable fraction of all-cause mortality (%) and attributable deaths per 1,000 child-years, by age band:")
print(transform(anchors,
                af = sprintf("%.1f (%.1f-%.1f)", 100 * af, 100 * af_lo, 100 * af_hi),
                attributable = sprintf("%.1f (%.1f-%.1f)", attributable_per_1000cy, attributable_lo, attributable_hi),
                allcause = sprintf("%.0f", allcause_rate_per_1000cy),
                reference = sprintf("%.0f", reference_rate_per_1000cy))[
                  , c("age_band", "pfpr", "af", "attributable", "allcause", "reference")],
      row.names = FALSE, right = FALSE)

## ---- figure ---------------------------------------------------------------------------------------------
palette <- c("#7570B3", "#1B9E77", "#66A61E", "#E6AB02", "#D95F02", "#A6761D")
names(palette) <- AGE6B
plot_af <- ggplot2::ggplot(curves, ggplot2::aes(pfpr, 100 * af, colour = age_band, fill = age_band)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * af_lo, ymax = 100 * af_hi), alpha = 0.12, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::scale_colour_manual(values = palette, name = "Age band") +
  ggplot2::scale_fill_manual(values = palette, name = "Age band") +
  ggplot2::labs(x = expression(PfPR[2-10]~"(%)"), y = "Attributable fraction of all-cause mortality (%)",
                title = "Share of all-cause deaths attributable to malaria, by age",
                subtitle = "Relative to 1% prevalence; smooth dose-response, one model per age band; 95% intervals") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "right")
plot_rate <- ggplot2::ggplot(curves, ggplot2::aes(pfpr, attributable_per_1000cy, colour = age_band, fill = age_band)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = attributable_lo, ymax = attributable_hi), alpha = 0.12, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::scale_colour_manual(values = palette, name = "Age band") +
  ggplot2::scale_fill_manual(values = palette, name = "Age band") +
  ggplot2::coord_cartesian(ylim = c(-2, 16)) +
  ggplot2::labs(x = expression(PfPR[2-10]~"(%)"), y = "Attributable deaths per 1,000 child-years",
                title = "Malaria-attributable all-cause mortality rate, by age",
                subtitle = "Rate at prevalence p minus the rate at 1%, for an average cell.\nThe neonatal band's interval (about -25 to +55) is cut off: null association, high baseline rate.") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "right")
plot_faceted <- ggplot2::ggplot(curves, ggplot2::aes(pfpr, 100 * af)) +
  ggplot2::geom_rect(data = support, ggplot2::aes(xmin = p01, xmax = p99, ymin = -Inf, ymax = Inf),
                     inherit.aes = FALSE, fill = "grey92") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * af_lo, ymax = 100 * af_hi, fill = age_band), alpha = 0.25, colour = NA) +
  ggplot2::geom_line(ggplot2::aes(colour = age_band), linewidth = 0.9) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::facet_wrap(~ age_band, nrow = 1) +
  ggplot2::scale_colour_manual(values = palette, guide = "none") +
  ggplot2::scale_fill_manual(values = palette, guide = "none") +
  ggplot2::labs(x = expression(PfPR[2-10]~"(%)"), y = "Attributable fraction (%)",
                subtitle = "Grey band: the central 98% of the prevalence values observed in that age band") +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(file.path(RESULTS_DIR, "figure36_age_specific_attributable.png"),
                (plot_af | plot_rate) / plot_faceted + patchwork::plot_layout(heights = c(1.2, 0.8)),
                width = 14, height = 10, dpi = 200, bg = "white")
message("Wrote person_time_age_specific_attributable.csv, person_time_age_specific_anchors.csv and figure36_age_specific_attributable.png")
