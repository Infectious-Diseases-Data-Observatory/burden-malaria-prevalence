# =============================================================================
# 17_sensitivity_timing.R — prevalence/mortality TIMING sensitivity analysis.
#
# THE PROBLEM. DHS.rates::chmort computes a synthetic-cohort PERIOD life table:
# with the default Period = 60 months the age-specific death probabilities come
# from deaths and child-time occurring in the 5 years BEFORE the interview, so
# the mortality estimate refers to a window centred ~2.5 years pre-survey. But
# each region-year is stamped with the SURVEY year and paired with the MAP
# PfPR2-10 surface for the SURVEY year. Because transmission declined over
# 2000-2024, survey-year prevalence is systematically LOWER than the prevalence
# that actually prevailed during the mortality window.
#
# TWO CORRECTIONS, applied separately and jointly:
#   (1) exposure side  — average/lag the annual MAP PfPR2-10 over the years the
#                        mortality window spans (re-extracted over the same DHS
#                        admin-1 polygons, same fixed GPW weighting);
#   (2) outcome side    — recompute chmort with Period = 24 or 36 months, moving
#                        the mortality window towards the survey date.
#
# CENTRE-MATCHING. A P-month mortality window is centred P/24 years before the
# interview; averaging annual prevalence over lags 0..L is centred L/2 years
# before. Matching the two therefore requires L = P/12:
#   Period 60 -> lags 0-5 (0-4 used; MAP starts in 2000), Period 36 -> lags 0-3,
#   Period 24 -> lags 0-2. Pairing a short mortality window with a long
#   prevalence window (or vice versa) re-introduces a mismatch in the opposite
#   direction and is reported here only to show the bracket.
#
# Both outcomes are analysed: post-neonatal (primary) and neonatal (negative
# control). Heavy inputs (cached annual MAP rasters, DHS boundaries, raw Births
# Recodes) are required only to (re)build the two long intermediates; both are
# cached per survey and skipped when the long CSVs already exist.
#
# Outputs (results/dhs_rebuild/):
#   sensitivity_prevalence_timing_summary.csv        exposure-side refits
#   sensitivity_mortality_horizon_summary.csv        outcome-side refits
#   sensitivity_timing_matched_combinations.csv      centre-matched combinations
#   sensitivity_prevalence_timing.png                exposure-side forest
#   sensitivity_prevalence_timing_splines.png        exposure-side dose-response
#   sensitivity_prevalence_timing_neonatal.png       negative control
#   sensitivity_mortality_horizon.png                outcome-side dose-response
#   sensitivity_timing_matched_combinations.png      matched-combination forest
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
PFPR_LAG_CACHE <- file.path(DATA_DIR, "map_lag_cache")
MORT_CACHE     <- file.path(DATA_DIR, "mort_period_cache")
MAX_LAG <- 4L                       # MAP has no pre-2000 surfaces
PERIODS <- c(24L, 36L, 60L)         # chmort reference periods, months

analysis_surveys <- function() {
  registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
  keep <- unique(analysis$svkey[as.logical(analysis$main_sample)])
  registry[registry$svkey %in% keep, , drop = FALSE]
}

## ===========================================================================
## Part 1 — lagged / windowed regional PfPR2-10 (exposure side)
## ===========================================================================
build_lagged_pfpr <- function() {
  required_packages(c("terra", "sf"))
  dir.create(PFPR_LAG_CACHE, showWarnings = FALSE)
  registry <- analysis_surveys()
  years <- sort(unique(unlist(lapply(registry$year, function(y) (y - MAX_LAG):y))))
  years <- years[years >= DHS_START_YEAR & years <= 2024L]
  rasters <- setNames(lapply(years, function(y)
    terra::rast(file.path(MAP_RASTER_DIR, sprintf("pfpr2_10_%d.tif", y)))), as.character(years))
  reference <- rasters[[1]]
  # one shared population-weight grid: the GPW surface is time-invariant, exactly
  # as in script 02, so weights are identical across reference years
  weights_full <- terra::resample(
    terra::crop(terra::rast(GPW_TIF), reference), reference, method = "bilinear")

  extract_survey <- function(svrow) {
    cache <- file.path(PFPR_LAG_CACHE, paste0(svrow$svkey, ".rds"))
    if (file.exists(cache)) return(readRDS(cache))
    boundary_file <- file.path(BOUNDARY_DIR, paste0(svrow$SurveyId, ".rds"))
    if (!file.exists(boundary_file)) stop("boundary unavailable")
    boundary <- sf::st_make_valid(readRDS(boundary_file))
    if (!"DHSREGEN" %in% names(boundary)) stop("boundary lacks DHSREGEN")
    polygons <- terra::makeValid(terra::vect(boundary))
    refs <- (svrow$year - MAX_LAG):svrow$year
    refs <- refs[as.character(refs) %in% names(rasters)]
    if (!length(refs)) stop("no MAP surface in window")
    # stack the window years and extract once: the exact fractional-coverage
    # step is the expensive part and is identical across years
    stack <- terra::rast(rasters[as.character(refs)])
    w <- terra::mask(weights_full, stack[[terra::nlyr(stack)]])
    denominator <- terra::extract(
      w, polygons, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)[[1]]
    numerator <- terra::extract(
      stack * w, polygons, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)
    regkeys <- rkey(boundary$DHSREGEN)
    out <- do.call(rbind, lapply(seq_along(refs), function(j) data.frame(
      svkey = svrow$svkey, iso3 = svrow$iso3, survey_year = svrow$year,
      ref_year = refs[j], lag = svrow$year - refs[j], regkey = regkeys,
      pfpr2_10 = 100 * numerator[[j]] / denominator, stringsAsFactors = FALSE)))
    out <- out[nzchar(out$regkey) & is.finite(out$pfpr2_10), , drop = FALSE]
    saveRDS(out, cache)
    out
  }

  message("Extracting lagged PfPR2-10 for ", nrow(registry), " surveys ...")
  rows <- list()
  for (i in seq_len(nrow(registry))) {
    svrow <- registry[i, , drop = FALSE]
    result <- tryCatch(extract_survey(svrow), error = function(e) {
      message("  skip ", svrow$svkey, ": ", conditionMessage(e)); NULL })
    if (!is.null(result) && nrow(result)) rows[[svrow$svkey]] <- result
  }
  long <- do.call(rbind, rows)
  rownames(long) <- NULL
  write.csv(long, PFPR_LONG_CSV, row.names = FALSE)
  long
}

## ===========================================================================
## Part 2 — region mortality at shorter chmort reference periods (outcome side)
## ===========================================================================
build_period_mortality <- function() {
  required_packages("DHS.rates")
  dir.create(MORT_CACHE, showWarnings = FALSE)
  registry <- analysis_surveys()
  # reformatted recodes may sit in the project cache or the rdhs cache
  candidates <- c(
    list.files(file.path(DATA_DIR, "dhs"), pattern = "rds$", full.names = TRUE),
    list.files(path.expand("~/.rdhs_cache"), pattern = "rds$",
               full.names = TRUE, recursive = TRUE))
  stems <- toupper(sub("\\.rds$", "", basename(candidates)))
  recode_path <- function(no_extension) {
    hit <- which(stems == toupper(no_extension))
    if (length(hit)) candidates[hit[1]] else NA_character_
  }
  # same extraction as the pipeline's mortality helper: the U5MR and NNMR rows
  # of chmort, with WN as the birth-exposure denominator
  mortality_at <- function(br, period) {
    rates <- suppressMessages(DHS.rates::chmort(br, Class = "v024", Period = period))
    u5 <- rates[grepl("^U5MR", rownames(rates)), c("Class", "R", "WN")]
    names(u5)[2:3] <- c("u5mr", "exposure")
    nn <- rates[grepl("^NNMR", rownames(rates)), c("Class", "R")]
    names(nn)[2] <- "nnmr"
    merged <- merge(u5, nn, by = "Class")
    data.frame(regkey = rkey(merged$Class), period = period,
               u5mr = merged$u5mr, nnmr = merged$nnmr,
               postneonatal_mortality = merged$u5mr - merged$nnmr,
               exposure = merged$exposure, stringsAsFactors = FALSE)
  }
  extract_survey <- function(svrow) {
    cache <- file.path(MORT_CACHE, paste0(svrow$svkey, ".rds"))
    if (file.exists(cache)) return(readRDS(cache))
    path <- recode_path(svrow$no_extension)
    if (is.na(path)) stop("recode unavailable")
    br <- readRDS(path)
    if (!"v024" %in% names(br)) stop("recode lacks v024")
    out <- do.call(rbind, lapply(PERIODS, function(p)
      cbind(svkey = svrow$svkey, mortality_at(br, p))))
    out <- out[nzchar(out$regkey) & is.finite(out$postneonatal_mortality), , drop = FALSE]
    saveRDS(out, cache)
    out
  }

  message("Recomputing mortality at Period = ", paste(PERIODS, collapse = "/"),
          " months for ", nrow(registry), " surveys ...")
  rows <- list()
  for (i in seq_len(nrow(registry))) {
    svrow <- registry[i, , drop = FALSE]
    result <- tryCatch(extract_survey(svrow), error = function(e) {
      message("  skip ", svrow$svkey, ": ", conditionMessage(e)); NULL })
    if (!is.null(result) && nrow(result)) rows[[svrow$svkey]] <- result
  }
  long <- do.call(rbind, rows)
  rownames(long) <- NULL
  write.csv(long, MORT_LONG_CSV, row.names = FALSE)
  long
}

pfpr_long <- if (file.exists(PFPR_LONG_CSV)) {
  read.csv(PFPR_LONG_CSV, stringsAsFactors = FALSE)
} else {
  tryCatch(build_lagged_pfpr(), error = function(e) {
    message("Lagged PfPR extraction unavailable (", conditionMessage(e), ")."); NULL })
}
mort_long <- if (file.exists(MORT_LONG_CSV)) {
  read.csv(MORT_LONG_CSV, stringsAsFactors = FALSE)
} else {
  tryCatch(build_period_mortality(), error = function(e) {
    message("Period mortality recomputation unavailable (", conditionMessage(e), ")."); NULL })
}
if (is.null(pfpr_long)) {
  message("Skipping the timing sensitivity: no lagged prevalence available.")
  quit(save = "no", status = 0)
}
pfpr_long$k <- paste(pfpr_long$svkey, pfpr_long$regkey, sep = "|")
if (!is.null(mort_long)) mort_long$k <- paste(mort_long$svkey, mort_long$regkey, sep = "|")

## reproduction checks: lag 0 must equal the published prevalence, and
## Period = 60 must equal the published mortality
lag0 <- pfpr_long[pfpr_long$lag == 0, ]
main <- as.logical(analysis$main_sample)
check <- merge(analysis[main, c("k", "pfpr2_10")],
               setNames(lag0[, c("k", "pfpr2_10")], c("k", "refit")), by = "k")
cat(sprintf("check: lag-0 vs published PfPR      %d rows, mean abs diff %.3f pp, r = %.4f\n",
            nrow(check), mean(abs(check$pfpr2_10 - check$refit), na.rm = TRUE),
            cor(check$pfpr2_10, check$refit, use = "complete.obs")))
if (!is.null(mort_long)) {
  p60 <- mort_long[mort_long$period == 60, ]
  chk <- merge(analysis[main, c("k", "postneonatal_mortality")],
               setNames(p60[, c("k", "postneonatal_mortality")], c("k", "refit")), by = "k")
  cat(sprintf("check: Period-60 vs published rate %d rows, median abs diff %.3f, r = %.4f\n",
              nrow(chk), median(abs(chk$postneonatal_mortality - chk$refit), na.rm = TRUE),
              cor(chk$postneonatal_mortality, chk$refit, use = "complete.obs")))
}

## ---- shared fitting helpers -------------------------------------------------
# prevalence averaged over lags 0..maxlag, keyed by survey-region
window_mean <- function(maxlag) {
  sub <- pfpr_long[pfpr_long$lag <= maxlag, ]
  tapply(sub$pfpr2_10, sub$k, mean)
}
window_length <- function(maxlag) {
  sub <- pfpr_long[pfpr_long$lag <= maxlag, ]
  tapply(sub$pfpr2_10, sub$k, length)
}
single_lag <- function(lag) {
  sub <- pfpr_long[pfpr_long$lag == lag, ]
  setNames(sub$pfpr2_10, sub$k)
}

# assemble the main-sample frame with a chosen prevalence vector and, optionally,
# mortality/exposure recomputed at a chosen chmort period
make_frame <- function(prevalence, outcome, period = NULL, rows = main) {
  dd <- analysis[rows, , drop = FALSE]
  if (!is.null(period)) {
    sub <- mort_long[mort_long$period == period,
                     c("k", "postneonatal_mortality", "nnmr", "exposure")]
    idx <- match(dd$k, sub$k)
    dd$postneonatal_mortality <- sub$postneonatal_mortality[idx]
    dd$nnmr <- sub$nnmr[idx]
    dd$exposure <- sub$exposure[idx]
  }
  prev <- as.numeric(prevalence[dd$k])
  dd$pfpr2_10 <- prev
  dd$pfpr10 <- prev / 10
  dd[is.finite(dd[[outcome]]) & dd[[outcome]] > 0 &
       is.finite(dd$exposure) & dd$exposure > 0 &
       is.finite(dd$pfpr10) & dd$pfpr10 > 0, , drop = FALSE]
}

fit_pair <- function(dd, outcome) {
  linear <- model_summary_row(fit_ridge_gam(
    dd, outcome, catalog, "linear_no_interaction",
    method = "REML", preprocessing = bundle$preprocessing))
  spline <- fit_ridge_gam(
    dd, outcome, catalog, "spline_no_interaction",
    method = "REML", preprocessing = bundle$preprocessing)
  anchors <- af_from_model(spline$model, c(10, 30, 50), year_c = 0)
  list(model = spline$model, row = data.frame(
    n = nrow(dd), median_exposure = round(median(dd$exposure)),
    pct_change_per_10 = linear$pct_change_per_10,
    lo = linear$pct_change_lo, hi = linear$pct_change_hi,
    ci_width = linear$pct_change_hi - linear$pct_change_lo,
    p_value = linear$pfpr_p, af10 = 100 * anchors$af[1],
    af30 = 100 * anchors$af[2], af50 = 100 * anchors$af[3]))
}

PREV_GRID <- exp(seq(log(1), log(60), length.out = 220))
curves <- function(model, label) {
  prediction <- link_prediction(model, newdata_at_mean(model, PREV_GRID / 10, year_c = 0))
  fraction <- af_from_model(model, PREV_GRID, year_c = 0)
  list(
    rate = data.frame(label = label, prevalence = PREV_GRID,
                      rate = 1000 * exp(prediction$fit),
                      lo = 1000 * exp(prediction$fit - 1.96 * prediction$se),
                      hi = 1000 * exp(prediction$fit + 1.96 * prediction$se)),
    af = data.frame(label = label, prevalence = PREV_GRID, af = 100 * fraction$af))
}

BLUE <- "#08519c"
x_log <- ggplot2::scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50))
base_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

forest_plot <- function(res, reference_label, title, subtitle) {
  reference <- res[res$label == reference_label, ]
  res$row <- factor(res$label, levels = rev(res$label))
  ggplot2::ggplot(res, ggplot2::aes(pct_change_per_10, row)) +
    ggplot2::annotate("rect", xmin = reference$lo, xmax = reference$hi,
                      ymin = -Inf, ymax = Inf, fill = BLUE, alpha = 0.10) +
    ggplot2::geom_vline(xintercept = reference$pct_change_per_10,
                        linetype = "dashed", colour = BLUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lo, xmax = hi, colour = label == reference_label),
      orientation = "y", width = 0.25, linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(colour = label == reference_label), size = 3) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey30", "TRUE" = BLUE), guide = "none") +
    ggplot2::labs(
      x = expression("Change in post-neonatal mortality per +10 " * italic(Pf) * "PR"[2-10] * " points (%, 95% CI)"),
      y = NULL, title = title, subtitle = subtitle) +
    base_theme +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
}

overlay_plot <- function(rate_df, af_df, levels, colours, rate_label, titles, log_y = TRUE) {
  rate_df$label <- factor(rate_df$label, levels = levels)
  af_df$label <- factor(af_df$label, levels = levels)
  panel_a <- ggplot2::ggplot(rate_df, ggplot2::aes(prevalence, rate, colour = label, fill = label)) +
    ggplot2::geom_ribbon(data = rate_df[rate_df$label == levels[1], ],
                         ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    ggplot2::geom_line(linewidth = 1) + x_log +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%), log scale"),
                  y = rate_label, title = titles[1]) +
    base_theme +
    ggplot2::theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
                   legend.background = ggplot2::element_rect(
                     fill = scales::alpha("white", 0.7), colour = NA))
  if (log_y) panel_a <- panel_a + ggplot2::scale_y_log10()
  panel_b <- ggplot2::ggplot(af_df, ggplot2::aes(prevalence, af, colour = label)) +
    ggplot2::geom_line(linewidth = 1) + x_log +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, 60)) +
    ggplot2::labs(x = expression(italic(Pf) * "PR"[2-10] * " (%), log scale"),
                  y = "Attributable share of deaths\n(%, versus 1% prevalence)",
                  title = titles[2]) +
    base_theme
  panel_a | panel_b
}

save_plot <- function(plot, file, width, height) {
  ggplot2::ggsave(file.path(RESULTS_DIR, file), plot,
                  width = width, height = height, dpi = 320, bg = "white")
}
report <- function(res, caption, columns) {
  cat("\n=== ", caption, " ===\n", sep = "")
  out <- res
  for (nm in intersect(c("pct_change_per_10", "lo", "hi", "ci_width"), names(out))) {
    out[[nm]] <- round(out[[nm]], 1)
  }
  for (nm in intersect(c("af10", "af30", "af50"), names(out))) {
    out[[nm]] <- round(out[[nm]])
  }
  print(out[, columns], row.names = FALSE)
}

## ===========================================================================
## Arm 1 — exposure side: lag / average the prevalence (mortality Period = 60)
## ===========================================================================
exposure_defs <- list(
  "Survey year (primary)" = single_lag(0),
  "Lag 1 year" = single_lag(1),
  "Lag 2 years" = single_lag(2),
  "5-year window mean" = window_mean(4)
)
# a common sample keeps the comparison about timing rather than about which
# region-years survive each definition
common <- main
for (prev in exposure_defs) {
  value <- as.numeric(prev[analysis$k])
  common <- common & is.finite(value) & value > 0
}
cat(sprintf("\nexposure-side common sample: %d region-years\n", sum(common)))
window_n <- as.numeric(window_length(4)[analysis$k])
cat(sprintf("prevalence window length: mean %.2f years; %d region-years truncated below 5\n",
            mean(window_n[main], na.rm = TRUE), sum(window_n[main] < 5, na.rm = TRUE)))

exposure_rows <- list(); exposure_rate <- list(); exposure_af <- list()
neonatal_rows <- list(); neonatal_rate <- list(); neonatal_af <- list()
for (label in names(exposure_defs)) {
  prev <- exposure_defs[[label]]
  fit <- fit_pair(make_frame(prev, "postneonatal_mortality", rows = common),
                  "postneonatal_mortality")
  exposure_rows[[label]] <- cbind(label = label, fit$row)
  drawn <- curves(fit$model, label)
  exposure_rate[[label]] <- drawn$rate; exposure_af[[label]] <- drawn$af

  neo <- fit_pair(make_frame(prev, "nnmr", rows = common), "nnmr")
  neonatal_rows[[label]] <- cbind(label = label, neo$row)
  drawn <- curves(neo$model, label)
  neonatal_rate[[label]] <- drawn$rate; neonatal_af[[label]] <- drawn$af
}
exposure_res <- do.call(rbind, exposure_rows)
neonatal_res <- do.call(rbind, neonatal_rows)
write.csv(exposure_res, file.path(RESULTS_DIR, "sensitivity_prevalence_timing_summary.csv"), row.names = FALSE)
write.csv(neonatal_res, file.path(RESULTS_DIR, "sensitivity_prevalence_timing_neonatal_summary.csv"), row.names = FALSE)
report(exposure_res, "Exposure side: post-neonatal association by prevalence timing",
       c("label", "n", "pct_change_per_10", "lo", "hi", "af10", "af30", "af50"))
report(neonatal_res, "Exposure side: NEONATAL negative control",
       c("label", "n", "pct_change_per_10", "lo", "hi", "p_value"))

exposure_levels <- names(exposure_defs)
exposure_colours <- setNames(c(BLUE, "#d95f0e", "#d73027", "#238b45"), exposure_levels)
save_plot(forest_plot(exposure_res, "Survey year (primary)",
                      "Prevalence-mortality association under temporally-matched prevalence",
                      "DHS mortality covers 0-59 months pre-survey (centre ~2.5y before); band = survey-year 95% CI"),
          "sensitivity_prevalence_timing.png", 9, 4.4)
save_plot(overlay_plot(do.call(rbind, exposure_rate), do.call(rbind, exposure_af),
                       exposure_levels, exposure_colours,
                       "Post-neonatal mortality\n(per 1000, log scale)",
                       c("A  Spline dose-response", "B  Attributable fraction")),
          "sensitivity_prevalence_timing_splines.png", 11, 4.6)
save_plot(overlay_plot(do.call(rbind, neonatal_rate), do.call(rbind, neonatal_af),
                       exposure_levels, exposure_colours,
                       "Neonatal mortality (per 1000)",
                       c("A  Neonatal spline dose-response (negative control)",
                         "B  Neonatal attributable fraction"), log_y = FALSE),
          "sensitivity_prevalence_timing_neonatal.png", 11, 4.6)

## ===========================================================================
## Arm 2 — outcome side: shorten the mortality horizon (survey-year prevalence)
## ===========================================================================
if (!is.null(mort_long)) {
  horizon_defs <- c("60 months (primary, 0-5y)" = 60, "36 months (0-3y)" = 36, "24 months (0-2y)" = 24)
  horizon_rows <- list(); horizon_rate <- list(); horizon_af <- list()
  for (label in names(horizon_defs)) {
    period <- horizon_defs[[label]]
    fit <- fit_pair(make_frame(single_lag(0), "postneonatal_mortality", period = period),
                    "postneonatal_mortality")
    horizon_rows[[label]] <- cbind(label = label, fit$row)
    drawn <- curves(fit$model, label)
    horizon_rate[[label]] <- drawn$rate; horizon_af[[label]] <- drawn$af
  }
  horizon_res <- do.call(rbind, horizon_rows)
  write.csv(horizon_res, file.path(RESULTS_DIR, "sensitivity_mortality_horizon_summary.csv"), row.names = FALSE)
  report(horizon_res, "Outcome side: post-neonatal association by mortality horizon",
         c("label", "n", "median_exposure", "pct_change_per_10", "lo", "hi", "ci_width", "af30"))
  horizon_levels <- names(horizon_defs)
  save_plot(overlay_plot(do.call(rbind, horizon_rate), do.call(rbind, horizon_af),
                         horizon_levels, setNames(c(BLUE, "#d95f0e", "#d73027"), horizon_levels),
                         "Post-neonatal mortality\n(per 1000, log scale)",
                         c("A  Spline dose-response by mortality horizon", "B  Attributable fraction")),
            "sensitivity_mortality_horizon.png", 11, 4.6)

  ## =========================================================================
  ## Arm 3 — centre-matched combinations (and deliberate mismatches)
  ## =========================================================================
  # L = P/12 matches the two window centres; MAP's 2000 start caps L at 4, so the
  # Period-60 pairing is approximate (centre 2.0y versus 2.5y).
  combinations <- list(
    list(label = "A  Period 60 + survey-year prevalence (published primary)", period = 60, prev = single_lag(0)),
    list(label = "B  Period 60 + 5-year window (approximately matched)",      period = 60, prev = window_mean(4)),
    list(label = "C  Period 36 + 4-year window (centre-matched)",             period = 36, prev = window_mean(3)),
    list(label = "D  Period 24 + 3-year window (centre-matched)",             period = 24, prev = window_mean(2)),
    list(label = "E  Period 24 + survey-year prevalence (outcome fix only)",  period = 24, prev = single_lag(0)),
    list(label = "F  Period 24 + 5-year window (over-lagged)",                period = 24, prev = window_mean(4))
  )
  combination_rows <- lapply(combinations, function(cb) {
    fit <- fit_pair(make_frame(cb$prev, "postneonatal_mortality", period = cb$period),
                    "postneonatal_mortality")
    cbind(label = cb$label, fit$row)
  })
  combination_res <- do.call(rbind, combination_rows)
  write.csv(combination_res, file.path(RESULTS_DIR, "sensitivity_timing_matched_combinations.csv"), row.names = FALSE)
  report(combination_res, "Centre-matched mortality-horizon x prevalence-window combinations",
         c("label", "n", "median_exposure", "pct_change_per_10", "lo", "hi", "ci_width", "af10", "af30", "af50"))
  cat("NOTE: AIC and deviance are NOT comparable across mortality horizons (the outcome data differ).\n")
  save_plot(forest_plot(combination_res, combination_res$label[1],
                        "Mortality horizon x prevalence window: centre-matched combinations",
                        "Matched pairings (B, C, D) versus deliberate mismatches (A, E too recent; F over-lagged)"),
            "sensitivity_timing_matched_combinations.png", 11, 5)
} else {
  message("Mortality horizons skipped: no recomputed mortality available ",
          "(raw DHS Births Recodes required).")
}
cat("\nTiming sensitivity complete. Outputs in ", RESULTS_DIR, "\n", sep = "")
