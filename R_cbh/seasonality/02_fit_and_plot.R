#!/usr/bin/env Rscript
# =============================================================================
# Sahel mortality by calendar month: plot the data and fit a simple model.
#
# Reads the aggregate cells from 01_build_cells.R. For each age band (<1, 1-5,
# 6-11 and 12-23 completed months) separately:
#
#   Data. The pooled design-weighted death rate per 1,000 child-years in each
#   calendar month (sum of weighted deaths over sum of weighted person-months
#   across all Sahelian regions), with a 95% interval from a leave-one-survey-out
#   jackknife.
#
#   Model. A negative-binomial regression of the cell death count on calendar
#   month as a factor, with age in completed months as a factor (so the within-
#   band age composition, which shifts with birth seasonality, cannot masquerade
#   as a calendar effect), a survey random intercept, and log person-months as
#   the offset. Cells are survey x calendar month x age in months. The count is
#   the weighted rate applied to the unweighted person-months, as in the earlier
#   person-time work, so the information content is the months actually
#   observed. Month effects are reported as rate ratios against the annual
#   geometric mean, with a Wald test of the eleven month contrasts and the
#   peak-to-trough ratio. A cyclic cubic spline in month (bs = "cc", December
#   joining January) is fitted alongside for a smooth curve.
#
#   Reference. Births by calendar month, as a share of the annual mean, are drawn
#   as a dashed line. Because a death month is derived as birth month plus age,
#   any heaping of reported ages at death (6, 12, 18 months) copies the birth
#   seasonality onto the death months; a death pattern that merely tracks the
#   birth pattern should be read with that in mind.
#
# Outputs (results/cbh/seasonality_sahel_v1)
#   month_rates_and_effects.csv   pooled rates, jackknife intervals, model rate ratios, births
#   cyclic_curves.csv             the cyclic-spline curve by band
#   seasonality_tests.csv         Wald test, AIC, peak and trough by band
#   sahel_mortality_by_calendar_month.png
#   REPORT.md
# =============================================================================
source("R_cbh/load_pipeline.R")
source("R_cbh/seasonality/settings.R")
for (p in c("mgcv", "ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) stop("Package ", p, " is required.")

st <- cbh_seasonality_settings()
cells <- cbh_read_csv(file.path(st$out, "cells_survey_month_age.csv"))
births <- cbh_read_csv(file.path(st$out, "births_survey_month.csv"))
regions <- cbh_read_csv(file.path(st$out, "sahel_regions.csv"))
BANDS <- st$bands$band
MONTHS <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
cells$band <- factor(cells$band, levels = BANDS)

## ---- births by month, relative to the annual mean --------------------------------------------
births_month <- tapply(births$births_w, births$month, sum)
births_rel <- data.frame(month = as.integer(names(births_month)),
                         births_rel = as.numeric(births_month) / mean(births_month))

## ---- per band ---------------------------------------------------------------------------------
jackknife_rates <- function(d) {
  rate <- function(x) 12000 * tapply(x$deaths_w, factor(x$month, levels = 1:12), sum) /
    tapply(x$person_months_w, factor(x$month, levels = 1:12), sum)
  pooled <- rate(d)
  surveys <- unique(d$svkey)
  leave_out <- sapply(surveys, function(s) rate(d[d$svkey != s, ]))
  n <- length(surveys)
  se <- sqrt((n - 1) / n * rowSums((leave_out - rowMeans(leave_out))^2))
  data.frame(month = 1:12, rate = as.numeric(pooled), rate_lo = as.numeric(pooled - 1.96 * se),
             rate_hi = as.numeric(pooled + 1.96 * se), annual_rate = 12000 * sum(d$deaths_w) / sum(d$person_months_w))
}

fit_band <- function(band, exclude_ages = integer(0), variant = "all ages in band") {
  d <- cells[cells$band == band & cells$person_months_n > 0 & !cells$age_months %in% exclude_ages, ]
  d$deaths_eff <- ifelse(d$person_months_w > 0, round(d$deaths_w / d$person_months_w * d$person_months_n), 0)
  d$log_pm <- log(d$person_months_n)
  d$month_f <- factor(d$month, levels = 1:12)
  d$age_f <- factor(d$age_months)
  d$survey <- factor(d$svkey)
  age_term <- if (nlevels(d$age_f) > 1) "age_f" else NULL
  base <- paste(c(age_term, "s(survey, bs = 're')", "offset(log_pm)"), collapse = " + ")
  m0 <- mgcv::gam(as.formula(paste("deaths_eff ~", base)), family = mgcv::nb(), data = d, method = "REML")
  m1 <- mgcv::gam(as.formula(paste("deaths_eff ~ month_f +", base)), family = mgcv::nb(), data = d, method = "REML")
  m2 <- mgcv::gam(as.formula(paste("deaths_eff ~ s(month, bs = 'cc', k = 6) +", base)), family = mgcv::nb(),
                  data = d, method = "REML", knots = list(month = c(0.5, 12.5)))
  # month effects centred on the annual geometric mean
  idx <- grep("^month_f", names(coef(m1)))
  b <- coef(m1)[idx]; V <- vcov(m1)[idx, idx]
  A <- rbind(0, diag(11))                       # 12 months from 11 contrasts (January reference)
  C <- (diag(12) - matrix(1 / 12, 12, 12)) %*% A
  centred <- as.numeric(C %*% b); se <- sqrt(diag(C %*% V %*% t(C)))
  wald <- as.numeric(t(b) %*% solve(V) %*% b)
  peak <- which.max(centred); trough <- which.min(centred)
  cpt <- (A[peak, ] - A[trough, ])
  ratio <- exp(sum(cpt * b)); ratio_se <- sqrt(as.numeric(t(cpt) %*% V %*% cpt))
  # cyclic curve, centred over the year, random effect and age at reference removed
  grid <- data.frame(month = seq(0.5, 12.5, by = 0.1), log_pm = 0, survey = levels(d$survey)[1])
  if (!is.null(age_term)) grid$age_f <- factor(levels(d$age_f)[1], levels = levels(d$age_f))
  X <- predict(m2, grid, type = "lpmatrix")
  X[, grep("^s\\(survey\\)", colnames(X))] <- 0
  Xc <- sweep(X, 2, colMeans(X))
  curve <- data.frame(band = band, month = grid$month, log_rr = as.numeric(Xc %*% coef(m2)),
                      se = sqrt(rowSums((Xc %*% vcov(m2)) * Xc)))
  effects <- cbind(band = band, jackknife_rates(d),
                   data.frame(rr = exp(centred), rr_lo = exp(centred - 1.96 * se), rr_hi = exp(centred + 1.96 * se)))
  effects <- merge(effects, births_rel, by = "month")
  tests <- data.frame(band = band, variant = variant, excluded_ages = paste(exclude_ages, collapse = ","),
                      surveys = nlevels(d$survey), cells = nrow(d), deaths_n = sum(d$deaths_n),
                      deaths_w = sum(d$deaths_w), person_years_w = sum(d$person_months_w) / 12,
                      annual_rate_per_1000 = effects$annual_rate[1],
                      wald_chi2 = wald, df = 11, p_value = pchisq(wald, 11, lower.tail = FALSE),
                      aic_no_month = AIC(m0), aic_month = AIC(m1), aic_cyclic = AIC(m2),
                      peak_month = month.abb[peak], trough_month = month.abb[trough],
                      peak_trough_ratio = ratio, ratio_lo = exp(log(ratio) - 1.96 * ratio_se),
                      ratio_hi = exp(log(ratio) + 1.96 * ratio_se), cyclic_edf = sum(m2$edf[grep("month", names(m2$edf))]),
                      nb_theta = m1$family$getTheta(TRUE), stringsAsFactors = FALSE)
  # survey-level relative rates for the plot (surveys with enough deaths)
  s <- aggregate(cbind(deaths_w, person_months_w) ~ svkey + month, data = d, FUN = sum)
  totals <- aggregate(cbind(deaths_w, person_months_w) ~ svkey, data = d, FUN = sum)
  s$rel <- (s$deaths_w / s$person_months_w) / (totals$deaths_w / totals$person_months_w)[match(s$svkey, totals$svkey)]
  s <- s[s$svkey %in% totals$svkey[totals$deaths_w >= st$min_deaths_for_survey_line], ]
  s <- s[is.finite(s$rel) & s$rel > 0, ]
  s$band <- band
  list(effects = effects, curve = curve, tests = tests, survey_lines = s)
}
results <- lapply(BANDS, fit_band)
effects <- do.call(rbind, lapply(results, `[[`, "effects"))
curves <- do.call(rbind, lapply(results, `[[`, "curve"))
tests <- do.call(rbind, lapply(results, `[[`, "tests"))
survey_lines <- do.call(rbind, lapply(results, `[[`, "survey_lines"))
# Robustness to age-at-death heaping. Reported ages heap at 6, 12 and 18 months, and a
# heaped death's derived month is its birth month plus a round number, so it carries the
# birth-month pattern. Refit the month-factor model on the unheaped ages only.
robust <- list(fit_band("6-11 months", exclude_ages = 6L, variant = "excluding age 6 months"),
               fit_band("12-23 months", exclude_ages = c(12L, 18L), variant = "excluding ages 12 and 18 months"))
tests_robust <- do.call(rbind, lapply(robust, `[[`, "tests"))
effects$band <- factor(effects$band, levels = BANDS); curves$band <- factor(curves$band, levels = BANDS)
survey_lines$band <- factor(survey_lines$band, levels = BANDS)
effects <- effects[order(effects$band, effects$month), ]
cbh_atomic_csv(effects, file.path(st$out, "month_rates_and_effects.csv"))
cbh_atomic_csv(curves, file.path(st$out, "cyclic_curves.csv"))
cbh_atomic_csv(rbind(tests, tests_robust), file.path(st$out, "seasonality_tests.csv"))

## ---- pooled monthly pattern by regional prevalence stratum (descriptive) ----------------------
region_cells <- cbh_read_csv(file.path(st$out, "cells_survey_region_month_band.csv"))
region_cells$mean_pfpr <- regions$mean_pfpr_pct[match(paste(region_cells$svkey, region_cells$regkey),
                                                       paste(regions$svkey, regions$regkey))]
STRATA <- c("PfPR < 10%", "PfPR 10-25%", "PfPR >= 25%")
region_cells$stratum <- cut(region_cells$mean_pfpr, c(-Inf, 10, 25, Inf), labels = STRATA)
strata <- do.call(rbind, lapply(BANDS, function(b) do.call(rbind, lapply(STRATA, function(s) {
  z <- region_cells[region_cells$band == b & region_cells$stratum %in% s, ]
  if (!nrow(z)) return(NULL)
  rate <- tapply(z$deaths_w, factor(z$month, levels = 1:12), sum) /
    tapply(z$person_months_w, factor(z$month, levels = 1:12), sum)
  rr <- as.numeric(rate) / (sum(z$deaths_w) / sum(z$person_months_w))
  data.frame(band = b, stratum = s, regions = length(unique(paste(z$svkey, z$regkey))),
             deaths_w = sum(z$deaths_w), person_years_w = sum(z$person_months_w) / 12,
             month = 1:12, rate_ratio = rr, stringsAsFactors = FALSE)
}))))
cbh_atomic_csv(strata, file.path(st$out, "month_rates_by_prevalence_stratum.csv"))
strata_summary <- do.call(rbind, lapply(split(strata, list(strata$band, strata$stratum), drop = TRUE), function(z) data.frame(
  band = z$band[1], stratum = z$stratum[1], regions = z$regions[1], deaths_w = z$deaths_w[1],
  peak = month.abb[which.max(z$rate_ratio)], trough = month.abb[which.min(z$rate_ratio)],
  peak_trough_ratio = max(z$rate_ratio) / min(z$rate_ratio),
  aug_oct_mean_rr = mean(z$rate_ratio[8:10]), stringsAsFactors = FALSE)))
strata_summary <- strata_summary[order(factor(strata_summary$band, levels = BANDS), factor(strata_summary$stratum, levels = STRATA)), ]

## ---- console -------------------------------------------------------------------------------------
message(sprintf("Sahel: %d surveys, %d regions, %s countries", length(unique(regions$svkey)), nrow(regions),
                paste(sort(unique(regions$iso3)), collapse = " ")))
message("\nSeasonality by age band (month-factor model against none):")
print(transform(tests, deaths = round(deaths_w), person_years = round(person_years_w),
                annual_rate = round(annual_rate_per_1000, 1), wald = round(wald_chi2, 1), p = signif(p_value, 2),
                dAIC = round(aic_no_month - aic_month, 1),
                peak_to_trough = sprintf("%.2f (%.2f-%.2f)", peak_trough_ratio, ratio_lo, ratio_hi))[
                  , c("band", "surveys", "deaths", "person_years", "annual_rate", "wald", "p", "dAIC",
                      "peak_month", "trough_month", "peak_to_trough")], row.names = FALSE)
message("\nRobustness to age heaping (month-factor model on unheaped ages only):")
print(transform(tests_robust, deaths = round(deaths_w), wald = round(wald_chi2, 1), p = signif(p_value, 2),
                peak_to_trough = sprintf("%.2f (%.2f-%.2f)", peak_trough_ratio, ratio_lo, ratio_hi))[
                  , c("band", "variant", "deaths", "wald", "p", "peak_month", "trough_month", "peak_to_trough")], row.names = FALSE)
message("\nRate ratio by calendar month against the annual mean:")
rr <- reshape(transform(effects, rr = sprintf("%.2f", rr))[, c("band", "month", "rr")],
              idvar = "band", timevar = "month", direction = "wide")
names(rr) <- c("band", month.abb); print(rr, row.names = FALSE)
message("\nPooled monthly pattern by regional mean PfPR (descriptive, no age adjustment):")
print(transform(strata_summary, deaths = round(deaths_w), ratio = round(peak_trough_ratio, 2),
                aug_oct = round(aug_oct_mean_rr, 2))[, c("band", "stratum", "regions", "deaths", "peak", "trough", "ratio", "aug_oct")],
      row.names = FALSE)
message("\nBirths by month, relative to the annual mean: ",
        paste(month.abb, sprintf("%.2f", births_rel$births_rel), collapse = "; "))

## ---- figure ---------------------------------------------------------------------------------------
palette <- c("#7570B3", "#1B9E77", "#66A61E", "#E6AB02"); names(palette) <- BANDS
labels <- data.frame(band = factor(tests$band, levels = BANDS),
                     text = sprintf("peak %s, trough %s\nratio %.2f (%.2f-%.2f)\np %s", tests$peak_month, tests$trough_month,
                                    tests$peak_trough_ratio, tests$ratio_lo, tests$ratio_hi,
                                    ifelse(tests$p_value < 0.001, "< 0.001", sprintf("= %.3f", tests$p_value))))
top <- ggplot2::ggplot(effects, ggplot2::aes(month, rate)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = annual_rate), colour = "grey60", linetype = "dashed") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = rate_lo, ymax = rate_hi, colour = band), width = 0.25, show.legend = FALSE) +
  ggplot2::geom_point(ggplot2::aes(colour = band), size = 2.4, show.legend = FALSE) +
  ggplot2::facet_wrap(~ band, nrow = 1, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = 1:12, labels = MONTHS) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::labs(x = NULL, y = "Deaths per 1,000 child-years",
                title = "Sahelian survey regions: child mortality by calendar month",
                subtitle = sprintf(paste0("Pooled design-weighted rate across %d regions in %d surveys (%s).\n",
                                          "Bars are 95%% leave-one-survey-out jackknife intervals; the dashed line is the annual rate."),
                                   nrow(regions), length(unique(regions$svkey)), paste(sort(unique(regions$iso3)), collapse = ", "))) +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
bottom <- ggplot2::ggplot(effects, ggplot2::aes(month, rr)) +
  ggplot2::geom_line(data = survey_lines, ggplot2::aes(month, rel, group = svkey), colour = "grey80", linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = 1, colour = "grey55") +
  ggplot2::geom_line(data = effects, ggplot2::aes(month, births_rel), colour = "grey40", linetype = "dashed", linewidth = 0.6) +
  ggplot2::geom_ribbon(data = curves, ggplot2::aes(month, exp(log_rr), ymin = exp(log_rr - 1.96 * se), ymax = exp(log_rr + 1.96 * se), fill = band),
                       alpha = 0.18, inherit.aes = FALSE, show.legend = FALSE) +
  ggplot2::geom_line(data = curves, ggplot2::aes(month, exp(log_rr), colour = band), linewidth = 0.9, show.legend = FALSE) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = rr_lo, ymax = rr_hi, colour = band), width = 0.25, show.legend = FALSE) +
  ggplot2::geom_point(ggplot2::aes(colour = band), size = 2.2, show.legend = FALSE) +
  ggplot2::geom_text(data = labels, ggplot2::aes(x = 0.6, y = Inf, label = text), hjust = 0, vjust = 1.15, size = 2.7, colour = "grey25") +
  ggplot2::facet_wrap(~ band, nrow = 1) +
  ggplot2::scale_x_continuous(breaks = 1:12, labels = MONTHS, limits = c(0.5, 12.5)) +
  ggplot2::scale_y_log10(breaks = c(0.5, 0.7, 1, 1.5, 2, 3)) +
  ggplot2::coord_cartesian(ylim = c(0.45, 3.2)) +
  ggplot2::scale_colour_manual(values = palette) + ggplot2::scale_fill_manual(values = palette) +
  ggplot2::labs(x = "Calendar month", y = "Rate ratio against the annual mean (log scale)",
                title = "Model: negative-binomial regression on calendar month",
                subtitle = paste0("Points: month effects adjusted for age in months, with a survey random intercept and 95% CI; curve: cyclic spline in month.\n",
                                  "Thin grey lines: surveys with at least ", st$min_deaths_for_survey_line, " deaths in the band; ",
                                  "dashed: births by calendar month relative to their mean.")) +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(st$out, "sahel_mortality_by_calendar_month.png"), top / bottom,
                width = 14, height = 8.5, dpi = 200, bg = "white")

## ---- report -----------------------------------------------------------------------------------------
fmt <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
excl <- cbh_read_csv(file.path(st$out, "exclusions.csv"))
excl_tab <- aggregate(children ~ status, data = excl, FUN = sum); excl_tab <- excl_tab[order(-excl_tab$children), ]
report <- c(
  "# Sahelian child mortality by calendar month",
  "",
  sprintf("Regions with a boundary centroid at or above %g degrees north and west of %g degrees east, excluding %s: %d survey regions in %d surveys (%s).",
          st$latitude_min, st$longitude_max, paste(st$exclude_countries, collapse = ", "), nrow(regions),
          length(unique(regions$svkey)), paste(sort(unique(regions$iso3)), collapse = ", ")),
  "",
  "Month of death is derived as month of birth plus age at death in completed months (exact for deaths reported in days, within about half a month for deaths reported in months; deaths reported as one year cannot be placed and those children are excluded). Children with a DHS-imputed birth month are excluded. Exposure is every calendar month in the 60 months before the interview month in which the child was alive and aged 0 to 23 completed months.",
  "",
  "## Seasonality by age band",
  "",
  "| Band | Surveys | Deaths (weighted) | Child-years (weighted) | Annual rate per 1,000 | Wald chi-square (11 df) | p | AIC gain from month | Peak | Trough | Peak / trough (95% CI) |",
  "|---|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %.2f (%.2f to %.2f) |", tests$band, tests$surveys, fmt(tests$deaths_w, 0),
          fmt(tests$person_years_w, 0), fmt(tests$annual_rate_per_1000, 1), fmt(tests$wald_chi2, 1),
          ifelse(tests$p_value < 0.001, "< 0.001", fmt(tests$p_value, 3)), fmt(tests$aic_no_month - tests$aic_month, 1),
          tests$peak_month, tests$trough_month, tests$peak_trough_ratio, tests$ratio_lo, tests$ratio_hi),
  "",
  "## Robustness to age-at-death heaping",
  "",
  "Reported ages at death heap at 6, 12 and 18 months. A heaped death's derived month is its birth month plus a round number of months, so it carries the birth-month pattern rather than a death-month pattern. The month-factor model refitted on the unheaped ages only:",
  "",
  "| Band | Ages used | Deaths (weighted) | Wald chi-square (11 df) | p | Peak | Trough | Peak / trough (95% CI) |", "|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %s | %s | %s | %s | %.2f (%.2f to %.2f) |", tests_robust$band, tests_robust$variant, fmt(tests_robust$deaths_w, 0),
          fmt(tests_robust$wald_chi2, 1), ifelse(tests_robust$p_value < 0.001, "< 0.001", fmt(tests_robust$p_value, 3)),
          tests_robust$peak_month, tests_robust$trough_month, tests_robust$peak_trough_ratio, tests_robust$ratio_lo, tests_robust$ratio_hi),
  "",
  "## Rate ratio by calendar month against the annual mean",
  "",
  paste("| Band |", paste(month.abb, collapse = " | "), "|"),
  paste("|---|", paste(rep("---", 12), collapse = "|"), "|"),
  vapply(BANDS, function(b) { e <- effects[effects$band == b, ]; paste("|", b, "|", paste(sprintf("%.2f", e$rr), collapse = " | "), "|") }, character(1)),
  paste("| Births |", paste(sprintf("%.2f", births_rel$births_rel), collapse = " | "), "|"),
  "",
  "## Pooled monthly pattern by regional mean MAP prevalence",
  "",
  "Descriptive pooled rates without age adjustment, by tercile-like strata of each region's mean PfPR over the exposure years. Aug-Oct is the mean rate ratio over August to October.",
  "",
  "| Band | Stratum | Regions | Deaths (weighted) | Peak | Trough | Peak / trough | Aug-Oct |", "|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %s | %s | %s | %.2f | %.2f |", strata_summary$band, strata_summary$stratum, strata_summary$regions,
          fmt(strata_summary$deaths_w, 0), strata_summary$peak, strata_summary$trough, strata_summary$peak_trough_ratio,
          strata_summary$aug_oct_mean_rr),
  "",
  sprintf("Deaths reported as one year of age cannot be placed in a calendar month: %s children with such a report are excluded, against %s weighted deaths used in the 12-23 month band, so that band rests on the roughly half of its deaths that were reported in months.",
          fmt(sum(excl$children[excl$status == "death_month_unresolved_year_unit"]), 0),
          fmt(tests$deaths_w[tests$band == "12-23 months"], 0)),
  "",
  "## Children by status",
  "",
  "| Status | Children |", "|---|---|",
  sprintf("| %s | %s |", excl_tab$status, fmt(excl_tab$children, 0)),
  "",
  "Figure: `sahel_mortality_by_calendar_month.png`. Tables: `month_rates_and_effects.csv`, `cyclic_curves.csv`, `seasonality_tests.csv`, `month_rates_by_prevalence_stratum.csv`, `sahel_regions.csv`, `exclusions.csv`.")
writeLines(report, file.path(st$out, "REPORT.md"))
message("\nWrote month_rates_and_effects.csv, cyclic_curves.csv, seasonality_tests.csv, sahel_mortality_by_calendar_month.png and REPORT.md")
