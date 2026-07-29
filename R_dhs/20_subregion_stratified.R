# =============================================================================
# 20_subregion_stratified.R — does the prevalence-mortality relationship differ
# by sub-region, and is the time interaction regional?
#
# Scripts 17 and 19 established that a PfPR x calendar-year interaction is
# significant once the temporal resolution is sharpened, and that it is not
# explained by which countries are surveyed. This script asks WHERE it comes
# from, by refitting the model separately in West Africa, Central Africa and
# East & Southern Africa (UN M49 groupings, Eastern and Southern combined).
#
# Specification as in script 19: Period-24 chmort post-neonatal mortality and
# PfPR2-10 averaged over the previous 2 years, with the interaction also tested
# under the published Period-60 / survey-year specification for comparability.
# Unlike the period-stratified fits, every sub-region spans the full 2000-2024
# window, so the curves can be evaluated at the common reference point
# (year_c = 0, covariates at the pooled mean) and remain directly comparable.
#
# HEADLINE. The pooled interaction is an average of opposing regional trends:
# strongly positive in West Africa (the slope steepening over time) and, under the
# published specification, significantly NEGATIVE in East & Southern Africa (the
# slope flattening). Central Africa has too few surveys to say. A single pooled
# time interaction is therefore a misspecification in both directions, and so is
# assuming a single time-stable slope.
#
# Outputs (results/dhs_rebuild/):
#   subregion_stratified_summary.csv       dose-response by sub-region
#   subregion_stratified_interaction.csv   interaction tests by sub-region
#   subregion_stratified_splines.png       overlaid dose-response and AF
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2", "patchwork", "scales"))

if (!file.exists(MODEL_BUNDLE_RDS)) stop("Run script 04 first.")
bundle <- readRDS(MODEL_BUNDLE_RDS)
catalog <- bundle$catalog
catalog$included_in_main <- as.logical(catalog$included_in_main)
analysis <- read_analysis_data()
analysis$k <- paste(analysis$svkey, analysis$regkey, sep = "|")

PFPR_LONG_CSV <- file.path(DERIVED_DIR, "map_pfpr_lagged_long.csv")
MORT_LONG_CSV <- file.path(DERIVED_DIR, "mortality_by_period_long.csv")
if (!file.exists(PFPR_LONG_CSV) || !file.exists(MORT_LONG_CSV)) {
  message("Run R_dhs/17_sensitivity_timing.R first; skipping.")
  quit(save = "no", status = 0)
}
pfpr_long <- read.csv(PFPR_LONG_CSV, stringsAsFactors = FALSE)
pfpr_long$k <- paste(pfpr_long$svkey, pfpr_long$regkey, sep = "|")
mort_long <- read.csv(MORT_LONG_CSV, stringsAsFactors = FALSE)
mort_long$k <- paste(mort_long$svkey, mort_long$regkey, sep = "|")
window2 <- local({
  sub <- pfpr_long[pfpr_long$lag <= 1, ]
  tapply(sub$pfpr2_10, sub$k, mean)
})

# UN M49 groupings, Eastern and Southern Africa combined; stated explicitly
# because countrycode's sub-region destinations are absent in some versions
SUBREGION_MAP <- c(
  BEN = "West Africa", BFA = "West Africa", CIV = "West Africa", CPV = "West Africa",
  GHA = "West Africa", GIN = "West Africa", GMB = "West Africa", GNB = "West Africa",
  LBR = "West Africa", MLI = "West Africa", MRT = "West Africa", NER = "West Africa",
  NGA = "West Africa", SEN = "West Africa", SLE = "West Africa", TGO = "West Africa",
  AGO = "Central Africa", CAF = "Central Africa", CMR = "Central Africa",
  COD = "Central Africa", COG = "Central Africa", GAB = "Central Africa",
  GNQ = "Central Africa", STP = "Central Africa", TCD = "Central Africa",
  BDI = "East & Southern", BWA = "East & Southern", COM = "East & Southern",
  DJI = "East & Southern", ERI = "East & Southern", ETH = "East & Southern",
  KEN = "East & Southern", LSO = "East & Southern", MDG = "East & Southern",
  MOZ = "East & Southern", MWI = "East & Southern", NAM = "East & Southern",
  RWA = "East & Southern", SDN = "East & Southern", SOM = "East & Southern",
  SSD = "East & Southern", SWZ = "East & Southern", TZA = "East & Southern",
  UGA = "East & Southern", ZAF = "East & Southern", ZMB = "East & Southern",
  ZWE = "East & Southern")
SUBREGIONS <- c("West Africa", "Central Africa", "East & Southern")
SPECS <- c("linear_no_interaction", "spline_no_interaction", "linear_time_interaction",
           "spline_time_interaction", "full_te_surface")

# assemble the panel under either the sharp or the published specification
build_panel <- function(period, use_window) {
  dd <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
  if (!is.null(period)) {
    sub <- mort_long[mort_long$period == period,
                     c("k", "postneonatal_mortality", "exposure")]
    idx <- match(dd$k, sub$k)
    dd$postneonatal_mortality <- sub$postneonatal_mortality[idx]
    dd$exposure <- sub$exposure[idx]
  }
  if (use_window) {
    prev <- as.numeric(window2[dd$k])
    dd$pfpr2_10 <- prev
    dd$pfpr10 <- prev / 10
  }
  dd$subregion <- unname(SUBREGION_MAP[dd$iso3])
  dd[is.finite(dd$postneonatal_mortality) & dd$postneonatal_mortality > 0 &
       is.finite(dd$exposure) & dd$exposure > 0 &
       is.finite(dd$pfpr10) & dd$pfpr10 > 0 & !is.na(dd$subregion), , drop = FALSE]
}
fit_spec <- function(dd, spec, method = "REML") {
  tryCatch(fit_ridge_gam(dd, "postneonatal_mortality", catalog, spec,
                         method = method, preprocessing = bundle$preprocessing),
           error = function(e) NULL)
}

## ---- A. dose-response fitted separately in each sub-region -----------------
sharp <- build_panel(24, TRUE)
fits <- list(); summary_rows <- list()
for (region in SUBREGIONS) {
  z <- sharp[sharp$subregion == region, , drop = FALSE]
  spline <- fit_spec(z, "spline_no_interaction")
  if (is.null(spline)) { message("  ", region, ": spline failed to fit"); next }
  linear <- fit_spec(z, "linear_no_interaction")
  lin <- if (is.null(linear)) NULL else model_summary_row(linear)
  st <- summary(spline$model)$s.table
  sr <- grep("^s\\(pfpr10\\)", rownames(st))
  anchors <- af_from_model(spline$model, c(10, 30, 50), year_c = 0)
  fits[[region]] <- spline
  summary_rows[[region]] <- data.frame(
    subregion = region, n = nrow(z), countries = length(unique(z$iso3)),
    pfpr_median = median(z$pfpr2_10), pfpr_max = max(z$pfpr2_10),
    smooth_edf = st[sr, "edf"], smooth_p = st[sr, "p-value"],
    pct_change_per_10 = if (is.null(lin)) NA_real_ else lin$pct_change_per_10,
    lo = if (is.null(lin)) NA_real_ else lin$pct_change_lo,
    hi = if (is.null(lin)) NA_real_ else lin$pct_change_hi,
    af10 = 100 * anchors$af[1], af30 = 100 * anchors$af[2], af50 = 100 * anchors$af[3])
}
summary_res <- do.call(rbind, summary_rows)
write.csv(summary_res, file.path(RESULTS_DIR, "subregion_stratified_summary.csv"), row.names = FALSE)
cat("=== A. Dose-response by sub-region (Period 24 + 2-year window) ===\n")
print(within(summary_res, {
  pfpr_median <- round(pfpr_median, 1); pfpr_max <- round(pfpr_max, 1)
  smooth_edf <- round(smooth_edf, 2); smooth_p <- signif(smooth_p, 2)
  pct_change_per_10 <- round(pct_change_per_10, 1); lo <- round(lo, 1); hi <- round(hi, 1)
  af10 <- round(af10); af30 <- round(af30); af50 <- round(af50)
}), row.names = FALSE)

## ---- B. is the interaction significant within each sub-region? -------------
interaction_rows <- list()
for (variant in list(list(lab = "Period 24 + 2-year window", period = 24, window = TRUE),
                     list(lab = "Period 60 + survey-year", period = NULL, window = FALSE))) {
  panel <- build_panel(variant$period, variant$window)
  for (region in SUBREGIONS) {
    z <- panel[panel$subregion == region, , drop = FALSE]
    aic <- vapply(SPECS, function(s) {
      f <- fit_spec(z, s, "ML"); if (is.null(f)) NA_real_ else AIC(f$model) }, 0)
    fl <- fit_spec(z, "linear_time_interaction", "ML")
    fs <- fit_spec(z, "spline_time_interaction", "ML")
    linear_row <- if (is.null(fl)) rep(NA_real_, 4) else
      summary(fl$model)$p.table["pfpr10:year_c", ]
    ti_edf <- ti_p <- NA_real_
    if (!is.null(fs)) {
      st <- summary(fs$model)$s.table
      rr <- grep("^ti\\(pfpr10", rownames(st))
      if (length(rr)) { ti_edf <- st[rr, "edf"]; ti_p <- st[rr, "p-value"] }
    }
    interaction_rows[[paste(variant$lab, region)]] <- data.frame(
      specification = variant$lab, subregion = region, n = nrow(z),
      aic_best = names(which.min(aic)),
      pfpr_year_coef = linear_row[1], pfpr_year_p = linear_row[4],
      ti_edf = ti_edf, ti_p = ti_p)
  }
}
interaction_res <- do.call(rbind, interaction_rows)
write.csv(interaction_res, file.path(RESULTS_DIR, "subregion_stratified_interaction.csv"), row.names = FALSE)
cat("\n=== B. PfPR x calendar-year interaction within each sub-region ===\n")
print(within(interaction_res, {
  pfpr_year_coef <- round(pfpr_year_coef, 5); pfpr_year_p <- signif(pfpr_year_p, 3)
  ti_edf <- round(ti_edf, 2); ti_p <- signif(ti_p, 3)
}), row.names = FALSE)
cat("\nA positive coefficient means the prevalence-mortality slope steepened over\n",
    "time; a negative one that it flattened. The signs differ by sub-region, so a\n",
    "single pooled interaction averages opposing regional trends.\n", sep = "")

## ---- overlaid dose-response and attributable fraction ----------------------
# each sub-region spans the whole study period, so a common evaluation point is
# valid; curves are still clipped to each sub-region's observed prevalence range
COLOURS <- setNames(c("#08519c", "#e08214", "#238b45"), SUBREGIONS)
rate_curves <- list(); af_curves <- list()
for (region in names(fits)) {
  z <- sharp[sharp$subregion == region, ]
  grid <- exp(seq(log(max(1, min(z$pfpr2_10))), log(max(z$pfpr2_10)), length.out = 200))
  prediction <- link_prediction(fits[[region]]$model,
                                newdata_at_mean(fits[[region]]$model, grid / 10, year_c = 0))
  rate_curves[[region]] <- data.frame(subregion = region, prevalence = grid,
    rate = 1000 * exp(prediction$fit),
    lo = 1000 * exp(prediction$fit - 1.96 * prediction$se),
    hi = 1000 * exp(prediction$fit + 1.96 * prediction$se))
  af_curves[[region]] <- data.frame(subregion = region, prevalence = grid,
    af = 100 * af_from_model(fits[[region]]$model, grid, year_c = 0)$af)
}
rate_df <- do.call(rbind, rate_curves); af_df <- do.call(rbind, af_curves)
rate_df$subregion <- factor(rate_df$subregion, levels = SUBREGIONS)
af_df$subregion <- factor(af_df$subregion, levels = SUBREGIONS)
sharp$subregion <- factor(sharp$subregion, levels = SUBREGIONS)
x_log <- ggplot2::scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50))
base_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
panel_rate <- ggplot2::ggplot(rate_df, ggplot2::aes(prevalence, rate, colour = subregion, fill = subregion)) +
  ggplot2::geom_point(data = sharp, ggplot2::aes(pfpr2_10, postneonatal_mortality, colour = subregion),
                      inherit.aes = FALSE, alpha = 0.16, size = 0.8) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
  ggplot2::geom_line(linewidth = 1.1) + x_log + ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = COLOURS, name = NULL) +
  ggplot2::scale_fill_manual(values = COLOURS, guide = "none") +
  ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%, 2-year mean), log scale"),
                y = "Post-neonatal mortality\n(per 1000, log scale)",
                title = "A  Dose-response fitted separately by sub-region") +
  base_theme +
  ggplot2::theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
                 legend.background = ggplot2::element_rect(
                   fill = scales::alpha("white", 0.75), colour = NA))
panel_af <- ggplot2::ggplot(af_df, ggplot2::aes(prevalence, af, colour = subregion)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
  ggplot2::geom_line(linewidth = 1.1) + x_log +
  ggplot2::scale_colour_manual(values = COLOURS, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(0, 60)) +
  ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%, 2-year mean), log scale"),
                y = "Malaria-attributable share of\npost-neonatal deaths (%, versus 1%)",
                title = "B  Attributable fraction by sub-region") +
  base_theme
ggplot2::ggsave(file.path(RESULTS_DIR, "subregion_stratified_splines.png"),
                panel_rate | panel_af, width = 11.5, height = 4.8, dpi = 320, bg = "white")
cat("\nsaved: subregion_stratified_{summary,interaction}.csv + subregion_stratified_splines.png\n")
