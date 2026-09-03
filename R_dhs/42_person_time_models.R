# =============================================================================
# 42_person_time_models.R — the piecewise-exponential model on window-by-age
# person-time, paired with the MAP prevalence of each window's own year.
#
# Rows are region x 12-month window x DHS age segment (script 40). Each window is
# paired with the person-time-weighted average of the annual MAP surfaces for
# the calendar years it spans (script 41). Within an age band the model is a
# negative-binomial regression on deaths with log person-time as the offset:
#
#   log E[deaths] = log(person-months) + segment + window index
#                 + s(PfPR) + s(calendar year) + ridge covariate block
#                 + country intercept + survey intercept
#
# One model per age band, with a smooth prevalence effect. The dose-response is
# not linear: it rises steeply to about 30% prevalence and flattens above, so a
# single "% change per 10 points" is not an interpretable summary and is not
# produced. Each band's effect is reported as the hazard-ratio curve against 0%
# prevalence and as the attributable fraction 1 - 1/HR at 10, 30 and 50%.
# Every region is kept whatever its prevalence, which is why the reference is
# 0% (PERSON_TIME_AF_REFERENCE) and not the 1% floor of the region-level
# pipeline.
#
# Why one model per band. A single joint model with age-specific prevalence
# effects but SHARED covariate effects and random effects gave a strongly
# negative neonatal effect while the neonatal month on its own gives a null
# one: the shared regional level is calibrated to the older ages, where
# mortality rises with prevalence, and the neonatal effect is pushed down to
# compensate. Fitting each band separately lets every nuisance term differ by
# age and is the fully interacted version of "PfPR x age". Script 47 fits the
# joint model with band-specific nuisance terms as a check.
#
# Counts. Weighted deaths are not integers, so each cell's count is the
# design-weighted rate applied to the UNWEIGHTED person-months, rounded:
# deaths_eff = round(deaths_w / pm_w x pm_n), with pm_n as the offset. The rate
# is the weighted one; the information content is the person-months observed.
#
# Groupings: the analysis plan's four age groups (0-3 months, 3-12 months, 1-2
# years, 2-5 years), the neonatal month separated (five groups), and two
# six-band splits (<1, 1-2, 3-11 | <1, 1-3, 4-11, then 12-23, 24-35, 36-59
# months). The 4-month split is the main presentation. Windows centred before
# 2000 have no MAP surface and are dropped; windows straddling 1999/2000 or
# 2024/2025 borrow the nearest surface for that share.
#
# Fits use mgcv::bam with discretised covariates (a plain gam took eight minutes
# per fit; bam takes a second and reproduces it to four decimals).
#
# Outputs (results/dhs_rebuild)
#   person_time_model_data.csv        the modelling table
#   person_time_effects.csv           attributable fractions at 10/30/50% by age group, all groupings
#   person_time_window_effects.csv    recall (window index) effects by age group
#   person_time_sensitivity.csv       AF at 10/30/50% under: windows 1-4; window 1 only;
#                                     lag-1 prevalence; Poisson; no survey intercept
#   person_time_dose_response_curves.csv  hazard ratio against 0% on a prevalence grid
#   figure22_person_time_age_effects.png                 AF at 10/30/50% by age group
#   figure23_person_time_dose_response.png               four primary groups
#   figure26_person_time_dose_response_six_bands.png     <1, 1-2, 3-11, 12-23, 24-35, 36-59 months
#   figure28_person_time_dose_response_six_bands_4m.png  <1, 1-3, 4-11, 12-23, 24-35, 36-59 months
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

PERSON_TIME_CSV <- file.path(DERIVED_DIR, "person_time_region_window_segment.csv")
YEAR_SHARE_CSV <- file.path(DERIVED_DIR, "person_time_window_year_shares.csv")
MAP_WINDOW_CSV <- file.path(DERIVED_DIR, "map_pfpr_window_years.csv")
MODEL_DATA_CSV <- file.path(RESULTS_DIR, "person_time_model_data.csv")
ANCHORS <- c(10, 30, 50)
GRID <- seq(0, 80, by = 1)
FIRST_MAP_YEAR <- 2000L
LAST_MAP_YEAR <- 2024L
BAM_THREADS <- max(1L, min(8L, parallel::detectCores() - 2L))
AGE5 <- c("0 months", "1-2 months", "3-11 months", "12-23 months", "24-59 months")
# AGE6 (split at 3 months) and AGE6B (split at 4 months) come from 00_config.R.
PREVALENCE_TERM <- "s(pfpr10, k = 5)"

bundle <- readRDS(MODEL_BUNDLE_RDS)
year_center <- unique(bundle$year_center)[1]
catalog <- read.csv(COVARIATE_CSV, stringsAsFactors = FALSE)
catalog$included_in_main <- as.logical(catalog$included_in_main)

## ---- 1. pair each region-window with its prevalence and calendar time -----------
cells <- read.csv(PERSON_TIME_CSV, stringsAsFactors = FALSE)
shares <- read.csv(YEAR_SHARE_CSV, stringsAsFactors = FALSE)
map_years <- read.csv(MAP_WINDOW_CSV, stringsAsFactors = FALSE)
map_key <- paste(map_years$svkey, map_years$regkey, map_years$year)

pair_windows <- function(shift) {
  shares$map_year <- pmin(pmax(shares$calendar_year - shift, FIRST_MAP_YEAR), LAST_MAP_YEAR)
  shares$pfpr <- map_years$pfpr2_10[match(paste(shares$svkey, shares$regkey, shares$map_year), map_key)]
  out <- do.call(rbind, lapply(
    split(shares, list(shares$svkey, shares$regkey, shares$window), drop = TRUE),
    function(d) data.frame(
      svkey = d$svkey[1], regkey = d$regkey[1], window = d$window[1],
      pfpr = sum(d$share * d$pfpr),
      mid_year = sum(d$share * (d$calendar_year + 0.5)),
      pre_map_share = sum(d$share[d$calendar_year < FIRST_MAP_YEAR]),
      map_missing = any(is.na(d$pfpr)),
      stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}
pairing <- pair_windows(0L)
lag1 <- pair_windows(1L)
pairing$pfpr_lag1 <- lag1$pfpr[match(paste(pairing$svkey, pairing$regkey, pairing$window),
                                     paste(lag1$svkey, lag1$regkey, lag1$window))]
dropped_pre_map <- sum(pairing$pre_map_share > 0.5)
dropped_missing <- sum(pairing$map_missing & pairing$pre_map_share <= 0.5)
pairing <- pairing[pairing$pre_map_share <= 0.5 & !pairing$map_missing, , drop = FALSE]
message("Region-windows: ", nrow(pairing), " kept; ", dropped_pre_map,
        " centred before ", FIRST_MAP_YEAR, " dropped; ", dropped_missing,
        " without a MAP value dropped")

## ---- 2. the modelling table --------------------------------------------------------
analysis <- read_analysis_data()
region_cov <- analysis[, c("svkey", "regkey", "iso3",
                           grep("_analysis$", names(analysis), value = TRUE))]
model_data <- merge(cells, pairing, by = c("svkey", "regkey", "window"))
model_data <- merge(model_data, region_cov, by = c("svkey", "regkey", "iso3"))
model_data <- model_data[model_data$person_months_w > 0 & model_data$person_months_n > 0, ]
model_data$deaths_eff <- round(model_data$deaths_w / model_data$person_months_w *
                                 model_data$person_months_n)
model_data$log_pm <- log(model_data$person_months_n)
model_data$pfpr10 <- model_data$pfpr / 10
model_data$pfpr10_lag1 <- model_data$pfpr_lag1 / 10
model_data$year_c <- model_data$mid_year - year_center
model_data$window_f <- factor(model_data$window)
model_data$age_group <- factor(model_data$age_group, levels = AGE_GROUPS)
model_data$age5 <- factor(ifelse(model_data$seg_lo == 0, AGE5[1],
                          ifelse(model_data$seg_lo < 3, AGE5[2],
                          ifelse(model_data$seg_lo < 12, AGE5[3],
                          ifelse(model_data$seg_lo < 24, AGE5[4], AGE5[5])))), levels = AGE5)
model_data$age6 <- band_from_segment(model_data$seg_lo, c(1, 3, 12, 24, 36), AGE6)
model_data$age6b <- band_from_segment(model_data$seg_lo, c(1, 4, 12, 24, 36), AGE6B)
model_data$country <- factor(model_data$iso3)
model_data$survey <- factor(model_data$svkey)
# descriptive only: how much of the prevalence variation is within a region over
# its five windows against between regions
region_mean <- tapply(model_data$pfpr10, paste(model_data$svkey, model_data$regkey), mean)
model_data$pfpr10_between <- as.numeric(region_mean[paste(model_data$svkey, model_data$regkey)])
model_data$pfpr10_within <- model_data$pfpr10 - model_data$pfpr10_between

ridge <- make_ridge_matrix(model_data, catalog, bundle$preprocessing)
model_data$G <- ridge$matrix
penalty <- list(G = list(diag(ncol(ridge$matrix))))
write.csv(model_data[, setdiff(names(model_data), "G")], MODEL_DATA_CSV, row.names = FALSE)
message(sprintf(paste0(
  "Modelling table: %d cells; %d region-windows; %d regions; %d surveys; %d countries; ",
  "%.0f weighted deaths (%.0f effective) over %.2f million person-months"),
  nrow(model_data), nrow(pairing),
  length(unique(paste(model_data$svkey, model_data$regkey))),
  nlevels(model_data$survey), nlevels(model_data$country),
  sum(model_data$deaths_w), sum(model_data$deaths_eff), sum(model_data$person_months_n) / 1e6))
message("Prevalence across region-windows: ", paste(round(range(pairing$pfpr), 1), collapse = "-"),
        "%; ", sum(pairing$pfpr < 1), " region-windows below 1%; within-region SD ",
        round(10 * stats::sd(model_data$pfpr10_within), 1), " points against between-region SD ",
        round(10 * stats::sd(model_data$pfpr10_between), 1))

## ---- 3. fitting helpers ------------------------------------------------------------------
fit_group <- function(data, terms = PREVALENCE_TERM, family = mgcv::nb(), survey_re = TRUE) {
  data$segment_f <- droplevels(factor(data$segment))
  parts <- c(terms,
             if (nlevels(data$segment_f) > 1) "segment_f",
             if (length(unique(data$window)) > 1) "window_f",
             "s(year_c, k = 8)", "G", "s(country, bs = 're')",
             if (survey_re) "s(survey, bs = 're')",
             "offset(log_pm)")
  f <- as.formula(paste("deaths_eff ~", paste(parts, collapse = " + ")))
  mgcv::bam(f, family = family, method = "fREML", paraPen = penalty, data = data,
            discrete = TRUE, nthreads = BAM_THREADS)
}

# log hazard ratio of prevalence p against the reference, with its standard
# error, for an average cell: first segment of the band, first window, calendar
# year at the centre, covariates at their means, random effects at zero
log_hazard_ratio <- function(fit, data, p, variable = "pfpr10") {
  segment_levels <- levels(droplevels(factor(data$segment)))
  frame <- function(values) {
    out <- data.frame(segment_f = factor(segment_levels[1], levels = segment_levels),
                      window_f = factor("1", levels = levels(model_data$window_f)),
                      year_c = 0, log_pm = 0,
                      country = factor(levels(model_data$country)[1], levels = levels(model_data$country)),
                      survey = factor(levels(model_data$survey)[1], levels = levels(model_data$survey)))
    out <- out[rep(1, length(values)), , drop = FALSE]
    out[[variable]] <- values / 10
    out$G <- matrix(0, nrow(out), ncol(model_data$G), dimnames = list(NULL, colnames(model_data$G)))
    out
  }
  Xh <- predict(fit, frame(p), type = "lpmatrix", discrete = FALSE)
  Xl <- predict(fit, frame(rep(PERSON_TIME_AF_REFERENCE, length(p))), type = "lpmatrix", discrete = FALSE)
  re <- grep("^s\\(country\\)|^s\\(survey\\)", colnames(Xh))
  Xh[, re] <- 0; Xl[, re] <- 0
  dX <- Xh - Xl
  data.frame(pfpr = p, est = as.numeric(dX %*% coef(fit)),
             se = sqrt(rowSums((dX %*% vcov(fit)) * dX)))
}
# attributable fraction 1 - 1/HR at the anchor prevalences, interval from the
# log hazard ratio
anchor_row <- function(fit, data, variable = "pfpr10") {
  lhr <- log_hazard_ratio(fit, data, ANCHORS, variable)
  row <- data.frame(row.names = NULL)
  for (i in seq_along(ANCHORS)) {
    p <- ANCHORS[i]
    row[1, paste0("af", p)] <- 1 - exp(-lhr$est[i])
    row[1, paste0("af", p, "_lo")] <- 1 - exp(-(lhr$est[i] - 1.96 * lhr$se[i]))
    row[1, paste0("af", p, "_hi")] <- 1 - exp(-(lhr$est[i] + 1.96 * lhr$se[i]))
  }
  row
}

## ---- 4. per age group: effects, recall, sensitivities ---------------------------------------
run_grouping <- function(grouping, variable) {
  levels_g <- levels(model_data[[variable]])
  out <- list(effects = list(), windows = list(), sens = list(), fits = list())
  for (g in levels_g) {
    data <- model_data[model_data[[variable]] == g, , drop = FALSE]
    message(sprintf("  %s | %-13s %6d cells, %6.0f deaths", grouping, g, nrow(data),
                    sum(data$deaths_w)))
    fit <- fit_group(data)
    out$fits[[g]] <- fit
    out$effects[[g]] <- cbind(
      data.frame(grouping = grouping, age_group = g, cells = nrow(data),
                 deaths = sum(data$deaths_w), nb_theta = fit$family$getTheta(TRUE),
                 prevalence_edf = sum(fit$edf[grep("s\\(pfpr10\\)", names(coef(fit)))]),
                 stringsAsFactors = FALSE),
      anchor_row(fit, data))
    tab <- summary(fit)$p.table
    wr <- grep("^window_f", rownames(tab))
    out$windows[[g]] <- data.frame(
      grouping = grouping, age_group = g,
      window = c(1L, as.integer(sub("^window_f", "", rownames(tab)[wr]))),
      rate_ratio_vs_window1 = c(1, exp(tab[wr, "Estimate"])),
      lo = c(NA, exp(tab[wr, "Estimate"] - 1.96 * tab[wr, "Std. Error"])),
      hi = c(NA, exp(tab[wr, "Estimate"] + 1.96 * tab[wr, "Std. Error"])),
      stringsAsFactors = FALSE)
    sens_fits <- list(
      "primary (NB, windows 1-5, lag 0)" = list(fit = fit, data = data, variable = "pfpr10"),
      "windows 1-4 only" = list(data = data[data$window <= 4, ], variable = "pfpr10"),
      "window 1 only (12 months before interview)" = list(data = data[data$window == 1, ], variable = "pfpr10"),
      "prevalence lagged one year" = list(data = data[is.finite(data$pfpr10_lag1), ], variable = "pfpr10_lag1",
                                          terms = "s(pfpr10_lag1, k = 5)"),
      "Poisson likelihood" = list(data = data, variable = "pfpr10", family = poisson()),
      "no survey intercept" = list(data = data, variable = "pfpr10", survey_re = FALSE))
    for (name in names(sens_fits)) {
      spec <- sens_fits[[name]]
      sfit <- if (!is.null(spec$fit)) spec$fit else
        fit_group(spec$data, terms = spec$terms %||% PREVALENCE_TERM,
                  family = spec$family %||% mgcv::nb(), survey_re = spec$survey_re %||% TRUE)
      out$sens[[paste(g, name)]] <- cbind(
        data.frame(grouping = grouping, age_group = g, sensitivity = name,
                   cells = nrow(spec$data), stringsAsFactors = FALSE),
        anchor_row(sfit, spec$data, spec$variable))
    }
  }
  out
}

message("\nPrimary grouping (four age groups)")
primary <- run_grouping("four groups", "age_group")
message("\nNeonatal month separated (five groups)")
split5 <- run_grouping("neonatal split", "age5")
message("\nSix bands, split at 3 months")
bands6 <- run_grouping("six bands", "age6")
message("\nSix bands, split at 4 months")
bands6b <- run_grouping("six bands (4-month split)", "age6b")

## ---- 5. tables -------------------------------------------------------------------------------
collect <- function(field) rbind(do.call(rbind, primary[[field]]),
                                 do.call(rbind, split5[[field]]),
                                 do.call(rbind, bands6[[field]]),
                                 do.call(rbind, bands6b[[field]]))
effects <- collect("effects"); windows <- collect("windows"); sens <- collect("sens")
rownames(effects) <- rownames(windows) <- rownames(sens) <- NULL

# post-neonatal (1-59 months) attributable fraction from the neonatal split,
# weighting the four groups from 1 month up by their share of deaths
pn <- effects[effects$grouping == "neonatal split" & effects$age_group != AGE5[1], ]
pn_weights <- pn$deaths / sum(pn$deaths)
postneonatal <- data.frame(
  grouping = "neonatal split", age_group = "1-59 months (death-share weighted)",
  cells = sum(pn$cells), deaths = sum(pn$deaths), stringsAsFactors = FALSE)
for (p in ANCHORS) {
  postneonatal[[paste0("af", p)]] <- sum(pn_weights * pn[[paste0("af", p)]])
  postneonatal[[paste0("af", p, "_lo")]] <- sum(pn_weights * pn[[paste0("af", p, "_lo")]])
  postneonatal[[paste0("af", p, "_hi")]] <- sum(pn_weights * pn[[paste0("af", p, "_hi")]])
}
effects <- merge(effects, postneonatal, all = TRUE, sort = FALSE)

write.csv(effects, file.path(RESULTS_DIR, "person_time_effects.csv"), row.names = FALSE)
write.csv(windows, file.path(RESULTS_DIR, "person_time_window_effects.csv"), row.names = FALSE)
write.csv(sens, file.path(RESULTS_DIR, "person_time_sensitivity.csv"), row.names = FALSE)

fmt_ci <- function(m, lo, hi, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), m, lo, hi)
show_af <- function(x) data.frame(
  x[, intersect(c("grouping", "age_group", "sensitivity"), names(x))], deaths = round(x$deaths),
  af10 = fmt_ci(100 * x$af10, 100 * x$af10_lo, 100 * x$af10_hi),
  af30 = fmt_ci(100 * x$af30, 100 * x$af30_lo, 100 * x$af30_hi),
  af50 = fmt_ci(100 * x$af50, 100 * x$af50_lo, 100 * x$af50_hi))
message("\nAttributable fraction of all-cause mortality (%) against ", PERSON_TIME_AF_REFERENCE,
        "% prevalence, by age group:")
print(cbind(show_af(effects), theta = round(effects$nb_theta, 1), edf = round(effects$prevalence_edf, 2)),
      row.names = FALSE)
message("\nRecorded mortality by window relative to the year before interview:")
print(transform(windows[windows$grouping == "four groups", ],
                rate_ratio_vs_window1 = round(rate_ratio_vs_window1, 3),
                lo = round(lo, 3), hi = round(hi, 3)), row.names = FALSE)
message("\nSensitivities (AF at 10/30/50%), six bands split at 4 months:")
s6 <- sens[sens$grouping == "six bands (4-month split)", ]
s6$deaths <- NA
print(show_af(s6)[, c("age_group", "sensitivity", "af10", "af30", "af50")], row.names = FALSE, right = FALSE)

## ---- 6. figures ------------------------------------------------------------------------------
show <- effects[!grepl("^six bands", effects$grouping) & !grepl("weighted", effects$age_group), ]
show_long <- do.call(rbind, lapply(ANCHORS, function(p) data.frame(
  grouping = show$grouping, age_group = show$age_group, anchor = paste0(p, "% prevalence"),
  af = show[[paste0("af", p)]], lo = show[[paste0("af", p, "_lo")]], hi = show[[paste0("af", p, "_hi")]],
  stringsAsFactors = FALSE)))
show_long$age_group <- factor(show_long$age_group, levels = c(AGE_GROUPS, AGE5))
show_long$grouping <- factor(show_long$grouping, levels = c("four groups", "neonatal split"),
                             labels = c("Primary: four age groups", "Neonatal month separated"))
plot_effects <- ggplot2::ggplot(show_long, ggplot2::aes(age_group, 100 * af, colour = anchor)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.15,
                         position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::geom_point(size = 2.6, position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::facet_wrap(~grouping, scales = "free_x") +
  ggplot2::scale_colour_manual(values = c("#7FB3C8", "#1D6F8B", "#0B3C4F"), name = NULL) +
  ggplot2::labs(x = NULL, y = sprintf("Attributable fraction of all-cause mortality (%%)\nagainst %d%% prevalence",
                                      PERSON_TIME_AF_REFERENCE),
                title = "Malaria-attributable fraction by age, person-time model",
                subtitle = paste("One negative-binomial model per age group with a smooth prevalence effect;",
                                 "country and survey intercepts,\nridge covariate block, window and calendar-year terms; 95% CIs")) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1), legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure22_person_time_age_effects.png"),
                plot_effects, width = 10, height = 4.8, dpi = 200)

curve_for <- function(fit, data, group) {
  lhr <- log_hazard_ratio(fit, data, GRID)
  data.frame(age_group = group, pfpr = GRID, hr = exp(lhr$est), lo = exp(lhr$est - 1.96 * lhr$se),
             hi = exp(lhr$est + 1.96 * lhr$se), stringsAsFactors = FALSE)
}
dose_response_figure <- function(result, variable, levels, file, title) {
  curves <- do.call(rbind, lapply(levels, function(g)
    curve_for(result$fits[[g]], model_data[model_data[[variable]] == g, ], g)))
  curves$age_group <- factor(curves$age_group, levels = levels)
  observed <- pairing$pfpr
  plot <- ggplot2::ggplot(curves, ggplot2::aes(pfpr, hr)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#1D6F8B", alpha = 0.15) +
    ggplot2::geom_line(colour = "#1D6F8B", linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_rug(data = data.frame(pfpr = observed[observed <= max(GRID)]),
                      ggplot2::aes(x = pfpr), inherit.aes = FALSE, alpha = 0.05, sides = "b") +
    ggplot2::facet_wrap(~age_group, nrow = 1) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "MAP PfPR2-10 in the window's own year (%)",
                  y = sprintf("Mortality hazard ratio\nversus %d%% prevalence (log scale)", PERSON_TIME_AF_REFERENCE),
                  title = title,
                  subtitle = "Person-time model, one fit per age band; bands are 95% CIs; rug shows region-window prevalence") +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(file, plot, width = 11, height = 3.8, dpi = 200)
  curves
}
curves4 <- dose_response_figure(primary, "age_group", AGE_GROUPS,
                                file.path(RESULTS_DIR, "figure23_person_time_dose_response.png"),
                                "Dose-response by age group, smooth prevalence effect")
curves6 <- dose_response_figure(bands6, "age6", AGE6,
                                file.path(RESULTS_DIR, "figure26_person_time_dose_response_six_bands.png"),
                                "Dose-response by age band, smooth prevalence effect")
curves6b <- dose_response_figure(bands6b, "age6b", AGE6B,
                                 file.path(RESULTS_DIR, "figure28_person_time_dose_response_six_bands_4m.png"),
                                 "Dose-response by age band (split at 4 months), smooth prevalence effect")
write.csv(rbind(cbind(grouping = "four groups", curves4), cbind(grouping = "six bands", curves6),
                cbind(grouping = "six bands (4-month split)", curves6b)),
          file.path(RESULTS_DIR, "person_time_dose_response_curves.csv"), row.names = FALSE)
saveRDS(list(smooth_fits = primary$fits,
             split_smooth_fits = split5$fits,
             bands6_smooth_fits = bands6$fits,
             bands6b_smooth_fits = bands6b$fits,
             year_center = year_center, af_reference = PERSON_TIME_AF_REFERENCE),
        file.path(DERIVED_DIR, "person_time_model_bundle.rds"))
message("\nWrote person_time_model_data.csv, person_time_effects.csv, person_time_window_effects.csv, ",
        "person_time_sensitivity.csv, person_time_dose_response_curves.csv and figures 22, 23, 26, 28")
