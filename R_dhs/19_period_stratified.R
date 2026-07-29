# =============================================================================
# 19_period_stratified.R — is the prevalence-mortality relationship stable over
# time? A NON-PARAMETRIC test by stratification.
#
# The time-varying question is usually asked by adding a PfPR x calendar-year
# interaction, but the answer then depends on how that interaction is
# parameterised (script 17, Arm 4: the 2000-2024 burden decline is 52-63% without
# an interaction and 1-38% with one). Here the study period is instead SPLIT into
# four equal intervals and the selected model is refitted independently within
# each, so any change in the dose-response emerges from the data rather than from
# an assumed functional form.
#
# Specification, following the timing analysis:
#   outcome    all-cause post-neonatal mortality from chmort with Period = 24
#              months, i.e. the 2-year lookback that keeps the mortality window
#              close to the survey date;
#   exposure   regional PfPR2-10 averaged over the previous 2 years (lags 0-1).
#              NOTE this sits 0.5 years later than the Period-24 window centre
#              (exact matching would use lags 0-2, script 17); the offset is
#              identical in every period, so it cancels in the BETWEEN-period
#              comparison that is the point of this script.
#   model      spline in PfPR2-10, no time interaction, ridge covariate block,
#              country random intercept and random PfPR slope, log-exposure
#              offset -- identical across periods, with the calendar-year smooth
#              reduced to k = 4 because each period spans only 6-7 years.
# The ridge standardisation is taken from the full-sample model bundle so that
# covariate coefficients are on a common scale and comparable across periods.
#
# TWO THREATS, both addressed rather than assumed away:
#   * COMPOSITION. Different countries are surveyed in different eras (only 8 of
#     the countries appear in all four periods). A difference between periods may
#     therefore reflect which countries contribute rather than calendar time, so
#     the analysis is repeated on countries present in at least three periods.
#   * EXTRAPOLATION. High-prevalence region-years become scarce later (the
#     maximum falls from about 81% to about 58%), so each fitted curve is drawn
#     only across the prevalence range actually observed in its own period.
#
# Outputs (results/dhs_rebuild/):
#   period_stratified_summary.csv            per-period fit summary
#   period_stratified_covariates.csv         ridge covariate effects by period
#   period_stratified_splines.png            overlaid dose-response and AF
#   period_stratified_scatter.png            four-panel prevalence vs mortality
#   period_stratified_covariate_forest.png   covariate effects by period
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
  message("Run R_dhs/17_sensitivity_timing.R first to build the lagged prevalence ",
          "and short-horizon mortality intermediates; skipping.")
  quit(save = "no", status = 0)
}
pfpr_long <- read.csv(PFPR_LONG_CSV, stringsAsFactors = FALSE)
pfpr_long$k <- paste(pfpr_long$svkey, pfpr_long$regkey, sep = "|")
mort_long <- read.csv(MORT_LONG_CSV, stringsAsFactors = FALSE)
mort_long$k <- paste(mort_long$svkey, mort_long$regkey, sep = "|")

## ---- assemble: Period-24 mortality + 2-year mean prevalence -----------------
window2 <- local({
  sub <- pfpr_long[pfpr_long$lag <= 1, ]
  tapply(sub$pfpr2_10, sub$k, mean)
})
short <- mort_long[mort_long$period == 24,
                   c("k", "postneonatal_mortality", "exposure")]
d <- analysis[as.logical(analysis$main_sample), , drop = FALSE]
idx <- match(d$k, short$k)
d$postneonatal_mortality <- short$postneonatal_mortality[idx]
d$exposure <- short$exposure[idx]
prev <- as.numeric(window2[d$k])
d$pfpr2_10 <- prev
d$pfpr10 <- prev / 10
d <- d[is.finite(d$postneonatal_mortality) & d$postneonatal_mortality > 0 &
         is.finite(d$exposure) & d$exposure > 0 &
         is.finite(d$pfpr10) & d$pfpr10 > 0, , drop = FALSE]

PERIOD_BREAKS <- c(1999.5, 2006.5, 2012.5, 2018.5, 2024.5)
PERIOD_LABELS <- c("2000-2006", "2007-2012", "2013-2018", "2019-2024")
d$period <- cut(d$year, breaks = PERIOD_BREAKS, labels = PERIOD_LABELS)
cat(sprintf("period-stratified sample: %d region-years, %d-%d\n",
            nrow(d), min(d$year), max(d$year)))

# countries contributing to at least three periods, for the composition check
country_periods <- tapply(as.character(d$period), d$iso3, function(x) length(unique(x)))
balanced_iso <- names(country_periods)[country_periods >= 3]
cat(sprintf("countries in all four periods: %d; in at least three: %d (%d of %d region-years)\n",
            sum(country_periods == 4), length(balanced_iso),
            sum(d$iso3 %in% balanced_iso), nrow(d)))

## ---- fit the selected model independently within each period ---------------
# year_k is reduced because a 6-7 year period cannot support the full calendar
# smooth; everything else matches the primary specification
fit_period <- function(rows, spec = "spline_no_interaction") {
  fit_ridge_gam(d[rows, , drop = FALSE], "postneonatal_mortality", catalog, spec,
                method = "REML", preprocessing = bundle$preprocessing, year_k = 4)
}

fits <- list(); summary_rows <- list()
for (p in PERIOD_LABELS) {
  rows <- d$period == p
  spline <- try(fit_period(rows), silent = TRUE)
  if (inherits(spline, "try-error")) {
    message("  period ", p, " failed to fit: ", conditionMessage(attr(spline, "condition")))
    next
  }
  linear <- try(fit_period(rows, "linear_no_interaction"), silent = TRUE)
  lin <- if (inherits(linear, "try-error")) NULL else model_summary_row(linear)
  st <- summary(spline$model)$s.table
  smooth_row <- grep("^s\\(pfpr10\\)", rownames(st))
  anchors <- af_from_model(spline$model, c(10, 30), year_c = 0)
  fits[[p]] <- spline
  summary_rows[[p]] <- data.frame(
    period = p, n = nrow(spline$data),
    countries = nlevels(droplevels(spline$data$country)),
    pfpr_min = min(d$pfpr2_10[rows]), pfpr_max = max(d$pfpr2_10[rows]),
    median_mortality = median(d$postneonatal_mortality[rows]),
    smooth_edf = if (length(smooth_row)) st[smooth_row[1], "edf"] else NA_real_,
    smooth_p = if (length(smooth_row)) st[smooth_row[1], "p-value"] else NA_real_,
    pct_change_per_10 = if (is.null(lin)) NA_real_ else lin$pct_change_per_10,
    lo = if (is.null(lin)) NA_real_ else lin$pct_change_lo,
    hi = if (is.null(lin)) NA_real_ else lin$pct_change_hi,
    af10 = 100 * anchors$af[1], af30 = 100 * anchors$af[2])
}
if (!length(fits)) stop("No period fitted.")
summary_res <- do.call(rbind, summary_rows)
write.csv(summary_res, file.path(RESULTS_DIR, "period_stratified_summary.csv"), row.names = FALSE)
cat("\n=== Model refitted independently within each period ===\n")
print(within(summary_res, {
  pfpr_min <- round(pfpr_min, 1); pfpr_max <- round(pfpr_max, 1)
  median_mortality <- round(median_mortality, 1); smooth_edf <- round(smooth_edf, 2)
  smooth_p <- signif(smooth_p, 2); pct_change_per_10 <- round(pct_change_per_10, 1)
  lo <- round(lo, 1); hi <- round(hi, 1); af10 <- round(af10); af30 <- round(af30)
}), row.names = FALSE)

# composition check: same fits restricted to countries present in >=3 periods
balanced_rows <- list()
for (p in PERIOD_LABELS) {
  rows <- d$period == p & d$iso3 %in% balanced_iso
  if (sum(rows) < 60) next
  linear <- try(fit_period(rows, "linear_no_interaction"), silent = TRUE)
  if (inherits(linear, "try-error")) next
  lin <- model_summary_row(linear)
  balanced_rows[[p]] <- data.frame(period = p, n = lin$n,
    pct_change_per_10 = lin$pct_change_per_10, lo = lin$pct_change_lo, hi = lin$pct_change_hi)
}
if (length(balanced_rows)) {
  balanced_res <- do.call(rbind, balanced_rows)
  write.csv(balanced_res, file.path(RESULTS_DIR, "period_stratified_balanced.csv"), row.names = FALSE)
  cat("\n=== Composition check: countries present in at least three periods ===\n")
  print(within(balanced_res, { pct_change_per_10 <- round(pct_change_per_10, 1)
    lo <- round(lo, 1); hi <- round(hi, 1) }), row.names = FALSE)
}

## ---- overlaid dose-response and attributable fraction ----------------------
# each curve is drawn only over the prevalence range observed in its own period,
# so the comparison never relies on extrapolation
PERIOD_COLOURS <- setNames(c("#08519c", "#3690c0", "#d95f0e", "#a50f15"), PERIOD_LABELS)
# A within-period model must be evaluated INSIDE its own period. Using the pooled
# centre (year_c = 0, i.e. about 2012) and the pooled covariate mean (G = 0) would
# extrapolate the period's calendar smooth and ridge block far outside their data
# and produce absurd levels, so each curve is evaluated at that period's own mean
# year and mean covariate profile. The attributable fraction is a contrast at
# fixed year and covariates, so it is unaffected by this choice.
newdata_in_period <- function(fit, pfpr10, year_c) {
  model <- fit$model
  n <- length(pfpr10)
  G <- fit$data$G
  out <- data.frame(pfpr10 = pfpr10, year_c = rep(year_c, length.out = n), exposure = 1,
                    country = factor(levels(model$model$country)[1],
                                     levels = levels(model$model$country)))
  out$G <- matrix(rep(colMeans(G), each = n), nrow = n, dimnames = list(NULL, colnames(G)))
  out
}
rate_curves <- list(); af_curves <- list()
for (p in names(fits)) {
  rows <- d$period == p
  centre_year_c <- mean(fits[[p]]$data$year_c)
  grid <- exp(seq(log(max(1, min(d$pfpr2_10[rows]))), log(max(d$pfpr2_10[rows])), length.out = 200))
  prediction <- link_prediction(fits[[p]]$model,
                                newdata_in_period(fits[[p]], grid / 10, centre_year_c))
  rate_curves[[p]] <- data.frame(period = p, prevalence = grid,
    rate = 1000 * exp(prediction$fit),
    lo = 1000 * exp(prediction$fit - 1.96 * prediction$se),
    hi = 1000 * exp(prediction$fit + 1.96 * prediction$se))
  af_curves[[p]] <- data.frame(period = p, prevalence = grid,
    af = 100 * af_from_model(fits[[p]]$model, grid, year_c = 0)$af)
}
rate_df <- do.call(rbind, rate_curves); af_df <- do.call(rbind, af_curves)
rate_df$period <- factor(rate_df$period, levels = PERIOD_LABELS)
af_df$period <- factor(af_df$period, levels = PERIOD_LABELS)
x_log <- ggplot2::scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50))
base_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

panel_rate <- ggplot2::ggplot(rate_df, ggplot2::aes(prevalence, rate, colour = period, fill = period)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(data = d, ggplot2::aes(pfpr2_10, postneonatal_mortality, colour = period),
                      alpha = 0.16, size = 0.8, inherit.aes = FALSE) +
  x_log + ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = PERIOD_COLOURS, name = NULL) +
  ggplot2::scale_fill_manual(values = PERIOD_COLOURS, guide = "none") +
  ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%, 2-year mean), log scale"),
                y = "Post-neonatal mortality\n(per 1000, log scale)",
                title = "A  Dose-response fitted separately within each period") +
  base_theme +
  ggplot2::theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
                 legend.background = ggplot2::element_rect(
                   fill = scales::alpha("white", 0.75), colour = NA))
panel_af <- ggplot2::ggplot(af_df, ggplot2::aes(prevalence, af, colour = period)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
  ggplot2::geom_line(linewidth = 1.1) + x_log +
  ggplot2::scale_colour_manual(values = PERIOD_COLOURS, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(0, 60)) +
  ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%, 2-year mean), log scale"),
                y = "Malaria-attributable share of\npost-neonatal deaths (%, versus 1%)",
                title = "B  Attributable fraction by period") +
  base_theme
ggplot2::ggsave(file.path(RESULTS_DIR, "period_stratified_splines.png"),
                panel_rate | panel_af, width = 11.5, height = 4.8, dpi = 320, bg = "white")

## ---- four-panel scatter: prevalence versus mortality within each period -----
# the same data as panel A above, but separated so the actual spread of points
# behind each fitted curve is visible; the grey dashed line is the pooled fit
# across all periods, so a panel departing from it is where the relationship moved
pooled <- try(fit_ridge_gam(d, "postneonatal_mortality", catalog, "spline_no_interaction",
                            method = "REML", preprocessing = bundle$preprocessing), silent = TRUE)
pooled_curve <- NULL
if (!inherits(pooled, "try-error")) {
  grid <- exp(seq(log(max(1, min(d$pfpr2_10))), log(max(d$pfpr2_10)), length.out = 200))
  pooled_fit <- link_prediction(pooled$model, newdata_at_mean(pooled$model, grid / 10, year_c = 0))
  pooled_curve <- data.frame(prevalence = grid, rate = 1000 * exp(pooled_fit$fit))
}
scatter <- ggplot2::ggplot(d, ggplot2::aes(pfpr2_10, postneonatal_mortality)) +
  ggplot2::geom_point(ggplot2::aes(size = exposure, colour = period), alpha = 0.35) +
  { if (!is.null(pooled_curve))
      ggplot2::geom_line(data = pooled_curve, ggplot2::aes(prevalence, rate),
                         inherit.aes = FALSE, colour = "grey35",
                         linetype = "dashed", linewidth = 0.7) } +
  ggplot2::geom_ribbon(data = rate_df, ggplot2::aes(prevalence, ymin = lo, ymax = hi, fill = period),
                       inherit.aes = FALSE, alpha = 0.18) +
  ggplot2::geom_line(data = rate_df, ggplot2::aes(prevalence, rate, colour = period),
                     inherit.aes = FALSE, linewidth = 1.1) +
  ggplot2::facet_wrap(~ period, nrow = 2) +
  ggplot2::scale_colour_manual(values = PERIOD_COLOURS, guide = "none") +
  ggplot2::scale_fill_manual(values = PERIOD_COLOURS, guide = "none") +
  ggplot2::scale_size_area(max_size = 4.5, name = "Births", labels = scales::comma) +
  x_log + ggplot2::scale_y_log10() +
  ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%, 2-year mean), log scale"),
                y = "All-cause post-neonatal mortality (per 1000 live births, log scale)",
                title = "Prevalence and post-neonatal mortality by period",
                subtitle = paste("Coloured line, fit within the period; grey dashed line,",
                                 "pooled fit across all periods")) +
  base_theme +
  ggplot2::theme(legend.position = "bottom",
                 strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
                 strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(RESULTS_DIR, "period_stratified_scatter.png"),
                scatter, width = 10, height = 8, dpi = 320, bg = "white")

## ---- covariate effects by period -------------------------------------------
# the ridge block enters as a single matrix term "G", so its per-variable
# coefficients are the G rows of the parametric table
covariate_rows <- list()
for (p in names(fits)) {
  ptab <- summary(fits[[p]]$model)$p.table
  gr <- grep("^G", rownames(ptab), value = TRUE)
  if (!length(gr)) next
  est <- ptab[gr, "Estimate"]; se <- ptab[gr, "Std. Error"]
  covariate_rows[[p]] <- data.frame(
    period = p, variable = sub("^G", "", gr),
    pct_change = 100 * (exp(est) - 1),
    lo = 100 * (exp(est - 1.96 * se) - 1),
    hi = 100 * (exp(est + 1.96 * se) - 1),
    stringsAsFactors = FALSE)
}
cov_df <- do.call(rbind, covariate_rows)
# same display names as script 05, so the panels are directly comparable
covariate_labels <- c(
  pct_urban = "Urban residence", dtp3_reg = "DTP3 coverage (DHS region)",
  measles = "Measles vaccination", facility = "Facility delivery",
  educ_yrs = "Maternal education", excl_bf = "Exclusive breastfeeding",
  birth_int = "Birth interval <24 months", mage1 = "Maternal age at first birth",
  imp_water = "Improved water", imp_sanit = "Improved sanitation",
  elec_dhs = "Household electricity", hib3_wuenic = "Hib3 coverage (WUENIC)",
  pcv3_wuenic = "PCV completion (WUENIC)", rotac_wuenic = "Rotavirus completion (WUENIC)",
  log_hiv_prev = "Child HIV prevalence (log, UNAIDS)", dtp3 = "DTP3 coverage (national)",
  log_gdp = "GDP per capita", hexp_gdp = "Health expenditure, % GDP",
  log_hexp_pc = "Health expenditure per capita", elec = "Electricity access (national)")
cov_df$label <- unname(covariate_labels[cov_df$variable])
cov_df$label[is.na(cov_df$label)] <- cov_df$variable[is.na(cov_df$label)]
cov_df$period <- factor(cov_df$period, levels = PERIOD_LABELS)
write.csv(cov_df, file.path(RESULTS_DIR, "period_stratified_covariates.csv"), row.names = FALSE)

# order covariates by their average effect so the panel reads top to bottom
ordering <- tapply(cov_df$pct_change, cov_df$label, mean, na.rm = TRUE)
cov_df$label <- factor(cov_df$label, levels = names(sort(ordering)))
forest <- ggplot2::ggplot(cov_df, ggplot2::aes(pct_change, label, colour = period)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), orientation = "y",
                         width = 0, linewidth = 0.5,
                         position = ggplot2::position_dodge(width = 0.7)) +
  ggplot2::geom_point(size = 1.9, position = ggplot2::position_dodge(width = 0.7)) +
  ggplot2::scale_colour_manual(values = PERIOD_COLOURS, name = NULL) +
  ggplot2::labs(x = "Change in post-neonatal mortality per 1 SD of the covariate (%, 95% CI)",
                y = NULL,
                title = "Covariate effects estimated separately within each period",
                subtitle = "Ridge-penalised conditional associations; shared standardisation across periods") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(colour = "grey93"),
                 legend.position = "top",
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(RESULTS_DIR, "period_stratified_covariate_forest.png"),
                forest, width = 9.5, height = 7.5, dpi = 320, bg = "white")
cat("\nsaved: period_stratified_{summary,covariates,balanced}.csv +",
    "period_stratified_{splines,scatter,covariate_forest}.png\n")
