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
#   Period 60 -> lags 0-5 (lags 0-4 are used here, leaving it under-matched by
#   half a year), Period 36 -> lags 0-3 (exact), Period 24 -> lags 0-2 (exact;
#   note this is THREE calendar years, because a 24-month span ending mid-year
#   overlaps three of them with weights 0.25/0.50/0.25 -- a two-year mean would
#   sit at 0.5 years and correct only half the gap).
#   Pairing a short mortality window with a long
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
#   mortality_horizon_agreement.csv                  60- versus 24-month windows
#   mortality_p60_vs_p24_by_subregion.png            the same, by sub-region
#   sensitivity_timing_burden_by_structure.csv       burden decline by horizon
#                                                    AND time structure
#
# HEADLINE CAVEAT. Arm 4 shows the 2000-2024 burden decline is governed chiefly by
# whether a PfPR x calendar-year interaction is admitted (52-63% without it,
# 1-38% with it) rather than by the mortality horizon. That interaction is not
# significant at Period 60 on the published sample but becomes so as the horizon
# shortens, so the reported decline rests on an assumption of time-stability that
# the timing analysis cannot itself settle.
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
# Exact centre-matching for Period 60 needs lags 0-5 (centre 2.5 years); lags 0-4
# give a centre of 2.0 years, so the Period-60 pairing below is under-matched by
# half a year and is labelled "approximately matched" throughout. Lag 5 is
# available for surveys from 2005 onwards but not for earlier ones (MAP surfaces
# start in 2000), and extending the window would deepen the truncation already
# affecting the earliest surveys, so the window is capped at 4.
MAX_LAG <- 4L
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


  ## ---- diagnostic: do the two mortality horizons agree? --------------------
  # A scatter of the 60-month against the 24-month estimate, one panel per
  # sub-region. If shortening the window merely added noise the points would
  # scatter symmetrically about equality; a systematic offset would mean the
  # horizon changes the LEVEL of the outcome and not just its precision.
  SUBREGION_MAP <- c(
    BEN="West Africa",BFA="West Africa",CIV="West Africa",CPV="West Africa",
    GHA="West Africa",GIN="West Africa",GMB="West Africa",GNB="West Africa",
    LBR="West Africa",MLI="West Africa",MRT="West Africa",NER="West Africa",
    NGA="West Africa",SEN="West Africa",SLE="West Africa",TGO="West Africa",
    AGO="Central Africa",CAF="Central Africa",CMR="Central Africa",COD="Central Africa",
    COG="Central Africa",GAB="Central Africa",GNQ="Central Africa",STP="Central Africa",
    TCD="Central Africa",
    BDI="East & Southern",BWA="East & Southern",COM="East & Southern",DJI="East & Southern",
    ERI="East & Southern",ETH="East & Southern",KEN="East & Southern",LSO="East & Southern",
    MDG="East & Southern",MOZ="East & Southern",MWI="East & Southern",NAM="East & Southern",
    RWA="East & Southern",SDN="East & Southern",SOM="East & Southern",SSD="East & Southern",
    SWZ="East & Southern",TZA="East & Southern",UGA="East & Southern",ZAF="East & Southern",
    ZMB="East & Southern",ZWE="East & Southern")
  horizons <- merge(
    mort_long[mort_long$period == 60, c("k", "postneonatal_mortality", "exposure")],
    mort_long[mort_long$period == 24, c("k", "postneonatal_mortality", "exposure")],
    by = "k", suffixes = c("_p60", "_p24"))
  idx <- match(horizons$k, analysis$k)
  horizons$iso3 <- analysis$iso3[idx]
  horizons <- horizons[!is.na(idx) & as.logical(analysis$main_sample)[idx] &
                         is.finite(horizons$postneonatal_mortality_p60) &
                         horizons$postneonatal_mortality_p60 > 0 &
                         is.finite(horizons$postneonatal_mortality_p24) &
                         horizons$postneonatal_mortality_p24 > 0, , drop = FALSE]
  horizons$subregion <- factor(unname(SUBREGION_MAP[horizons$iso3]),
    levels = c("West Africa", "Central Africa", "East & Southern"))
  horizons <- horizons[!is.na(horizons$subregion), , drop = FALSE]
  horizons$ratio <- horizons$postneonatal_mortality_p60 / horizons$postneonatal_mortality_p24
  horizon_res <- do.call(rbind, lapply(c(levels(horizons$subregion), "ALL"), function(sr) {
    z <- if (sr == "ALL") horizons else horizons[horizons$subregion == sr, ]
    data.frame(subregion = sr, n = nrow(z),
      median_period60 = median(z$postneonatal_mortality_p60),
      median_period24 = median(z$postneonatal_mortality_p24),
      median_ratio = median(z$ratio), pct_period60_higher = 100 * mean(z$ratio > 1),
      log_correlation = cor(log(z$postneonatal_mortality_p60), log(z$postneonatal_mortality_p24)))
  }))
  write.csv(horizon_res, file.path(RESULTS_DIR, "mortality_horizon_agreement.csv"), row.names = FALSE)
  cat("\n=== Agreement between the 60-month and 24-month mortality windows ===\n")
  print(within(horizon_res, { median_period60 <- round(median_period60, 1)
    median_period24 <- round(median_period24, 1); median_ratio <- round(median_ratio, 3)
    pct_period60_higher <- round(pct_period60_higher)
    log_correlation <- round(log_correlation, 3) }), row.names = FALSE)
  horizon_plot <- ggplot2::ggplot(horizons,
      ggplot2::aes(postneonatal_mortality_p24, postneonatal_mortality_p60)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(size = exposure_p60, colour = subregion), alpha = 0.4) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20",
                         fill = "grey70", linewidth = 0.7) +
    ggplot2::facet_wrap(~ subregion, nrow = 1) +
    ggplot2::scale_colour_manual(values = c("West Africa" = "#08519c",
      "Central Africa" = "#e08214", "East & Southern" = "#238b45"), guide = "none") +
    ggplot2::scale_size_area(max_size = 4.5, name = "Births", labels = scales::comma) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + ggplot2::coord_equal() +
    ggplot2::labs(x = "Post-neonatal mortality, 24-month window (per 1000, log scale)",
      y = "Post-neonatal mortality,\n60-month window (per 1000, log scale)",
      title = "DHS mortality over a 5-year versus a 2-year reference period",
      subtitle = "One point per survey-region; dashed line is equality") +
    base_theme +
    ggplot2::theme(legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(RESULTS_DIR, "mortality_p60_vs_p24_by_subregion.png"),
                  horizon_plot, width = 12, height = 5, dpi = 320, bg = "white")

  ## =========================================================================
  ## Arm 4 — burden consequence of the mortality horizon AND the time structure
  ## =========================================================================
  # The headline 2000-2024 decline turns out to depend far more on whether a
  # PfPR x calendar-year interaction is admitted than on the mortality horizon or
  # the dose-response shape. This matters because the interaction becomes
  # significant as the horizon shortens (ti() p = 0.16 at Period 60 on the
  # published sample, 0.0022 at Period 24), and because the AIC margin between
  # the additive and time-varying structures is small and unstable. Under an
  # interaction the attributable fraction RISES over the period, offsetting the
  # fall in prevalence and flattening the burden trajectory.
  national <- merge(
    read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"),
             stringsAsFactors = FALSE)[, c("iso3", "year", "allcause_1mo5y", "births")],
    read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE),
    by = c("iso3", "year"))
  national$region <- countrycode::countrycode(national$iso3, "iso3c", "region", warn = FALSE)
  national <- national[national$region == "Sub-Saharan Africa" &
                         is.finite(national$allcause_1mo5y) & is.finite(national$pfpr_pct) &
                         is.finite(national$births) &
                         national$year >= 2000 & national$year <= 2024, , drop = FALSE]
  year_center <- unique(bundle$year_center)[1]
  burden_series <- function(model) {
    high <- newdata_at_mean(model, national$pfpr_pct / 10, year_c = national$year - year_center)
    low <- newdata_at_mean(model, rep(AF_REFERENCE / 10, nrow(national)),
                           year_c = national$year - year_center)
    dX <- population_lpmatrix(model, high) - population_lpmatrix(model, low)
    af <- pmax(1 - exp(-as.numeric(dX %*% coef(model))), 0)
    tapply(af * national$allcause_1mo5y, national$year, sum)
  }
  af_at <- function(model, prevalence, year) {
    high <- newdata_at_mean(model, prevalence / 10, year_c = year - year_center)
    low <- newdata_at_mean(model, AF_REFERENCE / 10, year_c = year - year_center)
    dX <- population_lpmatrix(model, high) - population_lpmatrix(model, low)
    100 * max(1 - exp(-as.numeric(dX %*% coef(model))), 0)
  }
  # NOTE: every row here is REFITTED on the rows that survive this script's
  # filters, so the Period-60 additive row is a close approximation to, not an
  # exact reproduction of, the headline burden reported in Figure 4 (which uses
  # the originally fitted model bundle). Compare declines across rows, not the
  # absolute totals against the manuscript.
  structure_cases <- list(
    list(label = "Period 60, additive spline (primary structure)", period = NULL, spec = "spline_no_interaction"),
    list(label = "Period 60, spline x time interaction",   period = 60,   spec = "spline_time_interaction"),
    list(label = "Period 24, additive spline",             period = 24,   spec = "spline_no_interaction"),
    list(label = "Period 24, spline x time interaction",   period = 24,   spec = "spline_time_interaction"),
    list(label = "Period 24, linear x time interaction",   period = 24,   spec = "linear_time_interaction"),
    list(label = "Period 24, additive linear",             period = 24,   spec = "linear_no_interaction")
  )
  structure_rows <- lapply(structure_cases, function(cs) {
    dd <- make_frame(single_lag(0), "postneonatal_mortality", period = cs$period)
    model <- fit_ridge_gam(dd, "postneonatal_mortality", catalog, cs$spec,
                           method = "REML", preprocessing = bundle$preprocessing)$model
    series <- burden_series(model)
    data.frame(label = cs$label, n = nrow(dd),
               deaths_2000 = as.numeric(series[["2000"]]),
               deaths_2024 = as.numeric(series[["2024"]]),
               decline_pct = 100 * (1 - series[["2024"]] / series[["2000"]]),
               af30_2000 = af_at(model, 30, 2000), af30_2024 = af_at(model, 30, 2024))
  })
  structure_res <- do.call(rbind, structure_rows)
  write.csv(structure_res, file.path(RESULTS_DIR, "sensitivity_timing_burden_by_structure.csv"),
            row.names = FALSE)
  cat("\n=== Burden 2000-2024 by mortality horizon and time structure ===\n")
  print(within(structure_res, {
    deaths_2000 <- round(deaths_2000 / 1000); deaths_2024 <- round(deaths_2024 / 1000)
    decline_pct <- round(decline_pct); af30_2000 <- round(af30_2000, 1); af30_2024 <- round(af30_2024, 1)
  }), row.names = FALSE)
  cat("Deaths in thousands. The decline depends chiefly on whether a PfPR x year\n",
      "interaction is admitted (52-63% without, 1-38% with), not on the horizon.\n", sep = "")
} else {
  message("Mortality horizons skipped: no recomputed mortality available ",
          "(raw DHS Births Recodes required).")
}
cat("\nTiming sensitivity complete. Outputs in ", RESULTS_DIR, "\n", sep = "")
