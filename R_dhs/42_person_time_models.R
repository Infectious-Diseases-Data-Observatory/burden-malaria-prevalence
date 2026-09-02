# =============================================================================
# 42_person_time_models.R — the piecewise-exponential model on window-by-age
# person-time, paired with the MAP prevalence of each window's own year.
#
# Rows are region x 12-month window x DHS age segment (script 40). Each window is
# paired with the person-time-weighted average of the annual MAP surfaces for
# the calendar years it spans (script 41). Within an age group the model is a
# negative-binomial regression on deaths with log person-time as the offset:
#
#   log E[deaths] = log(person-months) + segment + window index
#                 + PfPR + s(calendar year) + ridge covariate block
#                 + country intercept + survey intercept
#
# One model per age group. A single joint model with age-specific prevalence
# slopes but SHARED covariate effects and random effects gave a strongly
# negative neonatal slope (-6% per 10 points) while the neonatal month on its
# own gives a null one: the shared regional level is calibrated to the older
# ages, where mortality rises with prevalence, and the neonatal slope is pushed
# down to compensate. Fitting each age group separately lets every nuisance term
# differ by age and is the fully interacted version of "PfPR x age". The joint
# model is kept as a sensitivity row so the artefact is on record.
#
# Counts. Weighted deaths are not integers, so each cell's count is the
# design-weighted rate applied to the UNWEIGHTED person-months, rounded:
# deaths_eff = round(deaths_w / pm_w x pm_n), with pm_n as the offset. The rate
# is the weighted one; the information content is the person-months observed.
#
# Choices fixed by the analysis plan: prevalence effects by age group (0-3
# months, 3-12 months, 1-2 years, 2-5 years); every region kept whatever its
# prevalence; regional DHS covariates as survey-level constants (the national
# series too, for now). The neonatal month is also separated out, because the
# negative-control logic applies to it and not to months 1-2. Windows centred
# before 2000 have no MAP surface and are dropped; windows straddling 1999/2000
# or 2024/2025 borrow the nearest surface for that share.
#
# Fits use mgcv::bam with discretised covariates (a plain gam took eight minutes
# per fit; bam takes a second and reproduces it to four decimals). fREML is the
# only criterion discretisation supports, so the within-group ladders are
# compared by mgcv's AIC under fREML.
#
# Outputs (results/dhs_rebuild)
#   person_time_model_data.csv        the modelling table
#   person_time_model_comparison.csv  AIC ladder within each age group
#   person_time_effects.csv           hazard ratios and attributable fractions by age group
#   person_time_within_between.csv    between-region against within-region prevalence effects
#   person_time_window_effects.csv    recall (window index) effects by age group
#   person_time_sensitivity.csv       windows 1-4; window 1 only; lag-1 prevalence; Poisson;
#                                     no survey intercept; joint model with shared nuisance terms
#   figure22_person_time_age_effects.png
#   figure23_person_time_dose_response.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "ggplot2"))

PERSON_TIME_CSV <- file.path(DERIVED_DIR, "person_time_region_window_segment.csv")
YEAR_SHARE_CSV <- file.path(DERIVED_DIR, "person_time_window_year_shares.csv")
MAP_WINDOW_CSV <- file.path(DERIVED_DIR, "map_pfpr_window_years.csv")
MODEL_DATA_CSV <- file.path(RESULTS_DIR, "person_time_model_data.csv")
ANCHORS <- c(10, 30, 50)
FIRST_MAP_YEAR <- 2000L
LAST_MAP_YEAR <- 2024L
BAM_THREADS <- max(1L, min(8L, parallel::detectCores() - 2L))
AGE5 <- c("0 months", "1-2 months", "3-11 months", "12-23 months", "24-59 months")

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
model_data$country <- factor(model_data$iso3)
model_data$survey <- factor(model_data$svkey)
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
fit_group <- function(data, terms, family = mgcv::nb(), survey_re = TRUE) {
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
slope <- function(fit, term = "pfpr10") {
  tab <- summary(fit)$p.table
  c(beta = tab[term, "Estimate"], se = tab[term, "Std. Error"], p = tab[term, 4])
}
pct <- function(b) 100 * (exp(b) - 1)
af <- function(beta, prevalence) 1 - exp(-beta * (prevalence - AF_REFERENCE) / 10)

STRUCTURES <- list(
  linear = list(label = "Linear prevalence effect", terms = "pfpr10"),
  smooth = list(label = "Smooth prevalence effect", terms = "s(pfpr10, k = 5)"),
  linear_time = list(label = "Linear effect changing with calendar year",
                     terms = "pfpr10 + pfpr10:year_c"),
  within_between = list(label = "Between-region and within-region prevalence",
                        terms = "pfpr10_between + pfpr10_within")
)

## ---- 4. per age group: ladder, effects, decomposition, recall, sensitivities ------------
run_grouping <- function(grouping, variable) {
  levels_g <- levels(model_data[[variable]])
  out <- list(comparison = list(), effects = list(), wb = list(), windows = list(),
              sens = list(), fits = list())
  for (g in levels_g) {
    data <- model_data[model_data[[variable]] == g, , drop = FALSE]
    message(sprintf("  %s | %-13s %6d cells, %6.0f deaths", grouping, g, nrow(data),
                    sum(data$deaths_w)))
    fits <- lapply(STRUCTURES, function(s) fit_group(data, s$terms))
    out$fits[[g]] <- fits
    aic <- vapply(fits, stats::AIC, numeric(1))
    out$comparison[[g]] <- data.frame(
      grouping = grouping, age_group = g, structure = names(STRUCTURES),
      label = vapply(STRUCTURES, `[[`, "", "label"),
      edf = vapply(fits, function(f) sum(f$edf), numeric(1)), AIC = aic,
      dAIC = aic - min(aic), stringsAsFactors = FALSE)
    s <- slope(fits$linear)
    row <- data.frame(grouping = grouping, age_group = g, cells = nrow(data),
                      deaths = sum(data$deaths_w), beta_per_10 = s[["beta"]],
                      se = s[["se"]], p_value = s[["p"]],
                      hr_per_10 = exp(s[["beta"]]),
                      hr_lo = exp(s[["beta"]] - 1.96 * s[["se"]]),
                      hr_hi = exp(s[["beta"]] + 1.96 * s[["se"]]),
                      nb_theta = fits$linear$family$getTheta(TRUE),
                      time_interaction_pct_per_year = pct(slope(fits$linear_time, "pfpr10:year_c")[["beta"]]),
                      time_interaction_p = slope(fits$linear_time, "pfpr10:year_c")[["p"]],
                      stringsAsFactors = FALSE)
    for (p in ANCHORS) {
      row[[paste0("af", p)]] <- af(s[["beta"]], p)
      row[[paste0("af", p, "_lo")]] <- af(s[["beta"]] - 1.96 * s[["se"]], p)
      row[[paste0("af", p, "_hi")]] <- af(s[["beta"]] + 1.96 * s[["se"]], p)
    }
    out$effects[[g]] <- row
    for (component in c("pfpr10_between", "pfpr10_within")) {
      w <- slope(fits$within_between, component)
      out$wb[[paste(g, component)]] <- data.frame(
        grouping = grouping, age_group = g,
        component = ifelse(component == "pfpr10_between", "between regions",
                           "within region over time"),
        pct_per_10 = pct(w[["beta"]]), pct_lo = pct(w[["beta"]] - 1.96 * w[["se"]]),
        pct_hi = pct(w[["beta"]] + 1.96 * w[["se"]]), p_value = w[["p"]],
        stringsAsFactors = FALSE)
    }
    tab <- summary(fits$linear)$p.table
    wr <- grep("^window_f", rownames(tab))
    out$windows[[g]] <- data.frame(
      grouping = grouping, age_group = g,
      window = c(1L, as.integer(sub("^window_f", "", rownames(tab)[wr]))),
      rate_ratio_vs_window1 = c(1, exp(tab[wr, "Estimate"])),
      lo = c(NA, exp(tab[wr, "Estimate"] - 1.96 * tab[wr, "Std. Error"])),
      hi = c(NA, exp(tab[wr, "Estimate"] + 1.96 * tab[wr, "Std. Error"])),
      stringsAsFactors = FALSE)
    sens_fits <- list(
      "primary (NB, windows 1-5, lag 0)" = fits$linear,
      "windows 1-4 only" = fit_group(data[data$window <= 4, ], "pfpr10"),
      "window 1 only (12 months before interview)" = fit_group(data[data$window == 1, ], "pfpr10"),
      "prevalence lagged one year" = fit_group(data[is.finite(data$pfpr10_lag1), ], "pfpr10_lag1"),
      "Poisson likelihood" = fit_group(data, "pfpr10", family = poisson()),
      "no survey intercept" = fit_group(data, "pfpr10", survey_re = FALSE))
    for (name in names(sens_fits)) {
      term <- if (grepl("lagged", name)) "pfpr10_lag1" else "pfpr10"
      w <- slope(sens_fits[[name]], term)
      out$sens[[paste(g, name)]] <- data.frame(
        grouping = grouping, age_group = g, sensitivity = name,
        pct_per_10 = pct(w[["beta"]]), pct_lo = pct(w[["beta"]] - 1.96 * w[["se"]]),
        pct_hi = pct(w[["beta"]] + 1.96 * w[["se"]]), stringsAsFactors = FALSE)
    }
  }
  out
}

message("\nPrimary grouping (four age groups)")
primary <- run_grouping("four groups", "age_group")
message("\nNeonatal month separated (five groups)")
split5 <- run_grouping("neonatal split", "age5")

## ---- 5. the joint model with shared nuisance terms, for the record -----------------------
message("\nJoint model with shared covariate effects and random effects")
model_data$segment_f <- factor(model_data$segment)
joint <- mgcv::bam(
  deaths_eff ~ pfpr10:age_group + segment_f + window_f + s(year_c, k = 8) + G +
    s(country, bs = "re") + s(survey, bs = "re") + offset(log_pm),
  family = mgcv::nb(), method = "fREML", paraPen = penalty, data = model_data,
  discrete = TRUE, nthreads = BAM_THREADS)
jt <- summary(joint)$p.table
jr <- grep("^pfpr10:age_group", rownames(jt))
joint_rows <- data.frame(
  grouping = "four groups", age_group = sub("^pfpr10:age_group", "", rownames(jt)[jr]),
  sensitivity = "joint model, shared random effects and covariates (artefact)",
  pct_per_10 = pct(jt[jr, "Estimate"]),
  pct_lo = pct(jt[jr, "Estimate"] - 1.96 * jt[jr, "Std. Error"]),
  pct_hi = pct(jt[jr, "Estimate"] + 1.96 * jt[jr, "Std. Error"]), stringsAsFactors = FALSE)

## ---- 6. tables -------------------------------------------------------------------------------
collect <- function(field) rbind(do.call(rbind, primary[[field]]), do.call(rbind, split5[[field]]))
comparison <- collect("comparison"); effects <- collect("effects"); wb <- collect("wb")
windows <- collect("windows"); sens <- rbind(collect("sens"), joint_rows)
rownames(comparison) <- rownames(effects) <- rownames(wb) <- rownames(windows) <- rownames(sens) <- NULL

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

write.csv(comparison, file.path(RESULTS_DIR, "person_time_model_comparison.csv"), row.names = FALSE)
write.csv(effects, file.path(RESULTS_DIR, "person_time_effects.csv"), row.names = FALSE)
write.csv(wb, file.path(RESULTS_DIR, "person_time_within_between.csv"), row.names = FALSE)
write.csv(windows, file.path(RESULTS_DIR, "person_time_window_effects.csv"), row.names = FALSE)
write.csv(sens, file.path(RESULTS_DIR, "person_time_sensitivity.csv"), row.names = FALSE)

fmt_ci <- function(m, lo, hi, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), m, lo, hi)
message("\nHazard ratio per +10 PfPR points and attributable fraction by age group:")
show <- effects[!is.na(effects$hr_per_10), ]
print(data.frame(grouping = show$grouping, age_group = show$age_group, deaths = round(show$deaths),
                 pct_per_10 = fmt_ci(pct(show$beta_per_10), pct(show$beta_per_10 - 1.96 * show$se),
                                     pct(show$beta_per_10 + 1.96 * show$se)),
                 p = signif(show$p_value, 2),
                 af10 = fmt_ci(100 * show$af10, 100 * show$af10_lo, 100 * show$af10_hi),
                 af30 = fmt_ci(100 * show$af30, 100 * show$af30_lo, 100 * show$af30_hi),
                 af50 = fmt_ci(100 * show$af50, 100 * show$af50_lo, 100 * show$af50_hi),
                 time_pct_per_year = round(show$time_interaction_pct_per_year, 2),
                 time_p = signif(show$time_interaction_p, 2)), row.names = FALSE)
message("Post-neonatal (1-59 months), death-share weighted, AF at 10/30/50%: ", paste(sprintf("%.1f (%.1f to %.1f)", 100 * unlist(postneonatal[paste0("af", ANCHORS)]),
                                             100 * unlist(postneonatal[paste0("af", ANCHORS, "_lo")]),
                                             100 * unlist(postneonatal[paste0("af", ANCHORS, "_hi")])), collapse = "; "))
message("\nModel ladder within each age group (dAIC):")
print(reshape(comparison[, c("grouping", "age_group", "structure", "dAIC")],
              idvar = c("grouping", "age_group"), timevar = "structure", direction = "wide"),
      row.names = FALSE, digits = 3)
message("\nBetween-region against within-region prevalence effects (% per +10 points):")
print(transform(wb, pct_per_10 = round(pct_per_10, 1), pct_lo = round(pct_lo, 1),
                pct_hi = round(pct_hi, 1), p_value = signif(p_value, 2)), row.names = FALSE)
message("\nRecorded mortality by window relative to the year before interview:")
print(transform(windows[windows$grouping == "four groups", ],
                rate_ratio_vs_window1 = round(rate_ratio_vs_window1, 3),
                lo = round(lo, 3), hi = round(hi, 3)), row.names = FALSE)
message("\nSensitivities (% per +10 points):")
print(transform(sens, pct_per_10 = round(pct_per_10, 1), pct_lo = round(pct_lo, 1),
                pct_hi = round(pct_hi, 1)), row.names = FALSE)

## ---- 7. figures ------------------------------------------------------------------------------
show$age_group <- factor(show$age_group, levels = c(AGE_GROUPS, AGE5))
show$grouping <- factor(show$grouping, levels = c("four groups", "neonatal split"),
                        labels = c("Primary: four age groups", "Neonatal month separated"))
plot_effects <- ggplot2::ggplot(show, ggplot2::aes(age_group, pct(beta_per_10))) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct(beta_per_10 - 1.96 * se),
                                      ymax = pct(beta_per_10 + 1.96 * se)),
                         width = 0.15, colour = "#1D6F8B") +
  ggplot2::geom_point(size = 3, colour = "#1D6F8B") +
  ggplot2::facet_wrap(~grouping, scales = "free_x") +
  ggplot2::labs(x = NULL, y = "Change in mortality hazard per +10 PfPR2-10 points (%)",
                title = "Prevalence effect by age, person-time model",
                subtitle = paste("One negative-binomial model per age group on deaths by region,",
                                 "12-month window and age segment;\ncountry and survey intercepts,",
                                 "ridge covariate block, window and calendar-year terms; 95% CIs")) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure22_person_time_age_effects.png"),
                plot_effects, width = 10, height = 4.4, dpi = 200)

curve_for <- function(fit, data, group) {
  p <- seq(0, 60, by = 1)
  segment_levels <- levels(droplevels(factor(data$segment)))
  frame <- function(pf) {
    out <- data.frame(pfpr10 = pf / 10,
                      segment_f = factor(segment_levels[1], levels = segment_levels),
                      window_f = factor("1", levels = levels(model_data$window_f)),
                      year_c = 0, log_pm = 0,
                      country = factor(levels(model_data$country)[1], levels = levels(model_data$country)),
                      survey = factor(levels(model_data$survey)[1], levels = levels(model_data$survey)))
    out$G <- matrix(0, nrow(out), ncol(model_data$G), dimnames = list(NULL, colnames(model_data$G)))
    out
  }
  Xh <- predict(fit, frame(p), type = "lpmatrix", discrete = FALSE)
  Xl <- predict(fit, frame(rep(AF_REFERENCE, length(p))), type = "lpmatrix", discrete = FALSE)
  re <- grep("^s\\(country\\)|^s\\(survey\\)", colnames(Xh))
  Xh[, re] <- 0; Xl[, re] <- 0
  dX <- Xh - Xl
  est <- as.numeric(dX %*% coef(fit)); se <- sqrt(rowSums((dX %*% vcov(fit)) * dX))
  data.frame(age_group = group, pfpr = p, hr = exp(est), lo = exp(est - 1.96 * se),
             hi = exp(est + 1.96 * se), stringsAsFactors = FALSE)
}
curves <- do.call(rbind, lapply(AGE_GROUPS, function(g)
  curve_for(primary$fits[[g]]$smooth, model_data[model_data$age_group == g, ], g)))
curves$age_group <- factor(curves$age_group, levels = AGE_GROUPS)
observed <- pairing$pfpr
plot_curves <- ggplot2::ggplot(curves, ggplot2::aes(pfpr, hr)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#1D6F8B", alpha = 0.15) +
  ggplot2::geom_line(colour = "#1D6F8B", linewidth = 0.9) +
  ggplot2::geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.4) +
  ggplot2::geom_rug(data = data.frame(pfpr = observed[observed <= 60]),
                    ggplot2::aes(x = pfpr), inherit.aes = FALSE, alpha = 0.05, sides = "b") +
  ggplot2::facet_wrap(~age_group, nrow = 1) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = "MAP PfPR2-10 in the window's own year (%)",
                y = "Mortality hazard ratio\nversus 1% prevalence (log scale)",
                title = "Dose-response by age group, smooth prevalence effect",
                subtitle = "Person-time model, one fit per age group; bands are 95% CIs; rug shows region-window prevalence") +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(file.path(RESULTS_DIR, "figure23_person_time_dose_response.png"),
                plot_curves, width = 11, height = 3.8, dpi = 200)

saveRDS(list(primary_fits = lapply(primary$fits, `[[`, "linear"),
             split_fits = lapply(split5$fits, `[[`, "linear"),
             smooth_fits = lapply(primary$fits, `[[`, "smooth"),
             joint = joint, year_center = year_center),
        file.path(DERIVED_DIR, "person_time_model_bundle.rds"))
message("\nWrote person_time_model_data.csv, person_time_model_comparison.csv, ",
        "person_time_effects.csv, person_time_within_between.csv, ",
        "person_time_window_effects.csv, person_time_sensitivity.csv, ",
        "figure22_person_time_age_effects.png, figure23_person_time_dose_response.png")
