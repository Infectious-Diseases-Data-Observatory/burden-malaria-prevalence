#!/usr/bin/env Rscript
# Report for the linear-PfPR sensitivity (records with PfPR 1-30%): hazard ratios over equal steps of
# the range, the linear slope against the average slope of the splines, curves against the v7 primary,
# AIC of linear against spline on the same records, and the national burden 2000-2024 with the linear
# model (log hazard linear in PfPR, extrapolated below 1% and above 30%). Reads saved components and
# aggregates only.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
suppressPackageStartupMessages({library(data.table); library(mgcv); library(ggplot2)})
base <- cbh_primary_settings("regional_mics")
id <- "linear_pfpr_1_30_dhsmics_map_gamma2_v2"; out <- file.path("results/cbh", id); private <- file.path("data/derived_cbh/models", id)
ages <- cbh_config()$age_bands$age_band
# Display label of the fitted primary (these reports were first written against v7).
pv <- sub(".*_(v[0-9]+)$", "\\1", base$id)
labels <- c(v7 = sprintf("Primary %s spline, all records", pv), spline = "Spline, PfPR 1–30% records", linear = "Linear, PfPR 1–30% records")
v7 <- readRDS(file.path(base$private, "pfpr_components.rds")); names(v7) <- ages[as.integer(sub(".*_", "", names(v7)))]
new <- readRDS(file.path(private, "pfpr_components.rds"))
pick <- function(m) { z <- new[vapply(new, function(x) x$model == m, NA)]; setNames(z, vapply(z, `[[`, "", "age_band")) }
comp <- list(v7 = v7, spline = pick("spline"), linear = pick("linear"))
stopifnot(all(vapply(comp, function(z) setequal(names(z), ages), NA)))
# Log hazard ratio for PfPR `from` -> `to` (vectors recycled) with its SE.
lhr <- function(piece, from, to) {
  n <- max(length(from), length(to)); from <- rep_len(from, n); to <- rep_len(to, n)
  if (is.null(piece$smooth)) return(list(est = unname(piece$coef) * (to - from), se = abs(to - from) * sqrt(piece$covariance[1, 1])))
  L <- PredictMat(piece$smooth, data.frame(pfpr_pct = to)) - PredictMat(piece$smooth, data.frame(pfpr_pct = from))
  list(est = drop(L %*% piece$coef), se = sqrt(pmax(rowSums((L %*% piece$covariance) * L), 0)))
}

# Hazard ratios for reductions across the range (the paper's X% -> Y% convention) and the slope.
steps <- list(c(30, 20), c(20, 10), c(10, 1), c(30, 1))
con <- rbindlist(lapply(names(comp), function(m) rbindlist(lapply(ages, function(a) rbindlist(lapply(steps, function(x) {
  k <- lhr(comp[[m]][[a]], x[1], x[2])
  data.table(model = m, age_band = a, contrast = sprintf("%g%% to %g%%", x[1], x[2]), log_hr = k$est, se = k$se,
    hr = exp(k$est), lower_95 = exp(k$est - 1.96 * k$se), upper_95 = exp(k$est + 1.96 * k$se))
}))))))
cbh_atomic_csv(as.data.frame(con), file.path(out, "contrasts.csv"))
# HR per 10-point increase: the linear coefficient, or the average slope of a spline from 1% to 30%.
slope <- con[contrast == "30% to 1%", .(model, age_band, log_hr_10 = -log_hr * 10 / 29, se_10 = se * 10 / 29)]
slope[, `:=`(hr_10 = exp(log_hr_10), lower_95 = exp(log_hr_10 - 1.96 * se_10), upper_95 = exp(log_hr_10 + 1.96 * se_10))]
cbh_atomic_csv(as.data.frame(slope), file.path(out, "slope_per_10_points.csv"))
af <- rbindlist(lapply(names(comp), function(m) rbindlist(lapply(ages, function(a) {
  k <- lhr(comp[[m]][[a]], c(10, 20, 30), 0)
  data.table(model = m, age_band = a, pfpr_pct = c(10, 20, 30), attributable_fraction = 1 - exp(k$est))
}))))
cbh_atomic_csv(as.data.frame(af), file.path(out, "attributable_fraction_by_age.csv"))

# Curves relative to PfPR 20%: the primary over 0-40%, the restricted models over their 1-30% range.
curves <- rbindlist(lapply(names(comp), function(m) rbindlist(lapply(ages, function(a) {
  g <- if (m == "v7") seq(0, 40, by = .25) else seq(1, 30, by = .25); k <- lhr(comp[[m]][[a]], 20, g)
  data.table(model = m, age_band = a, pfpr_pct = g, log_hr = k$est, lower = k$est - 1.96 * k$se, upper = k$est + 1.96 * k$se)
}))))
cbh_atomic_csv(as.data.frame(curves), file.path(out, "pfpr_curves.csv"))

# Fit: linear against spline on the same records.
fd <- fread(file.path(out, "fit_diagnostics.csv"))
fit <- dcast(fd, age_band ~ model, value.var = c("aic", "fREML", "loglik", "version", "children", "deaths"))[match(ages, age_band)]
stopifnot(all(fit$children_linear == fit$children_spline), all(fit$deaths_linear == fit$deaths_spline))
ss <- fread(file.path(out, "smooth_summaries.csv"))
fit <- merge(fit, ss[model == "spline" & term == "s(pfpr_pct)", .(age_band, spline_pfpr_edf = edf, spline_pfpr_p = `p-value`)], by = "age_band")[match(ages, age_band)]
fit[, `:=`(linear_p = vapply(comp$linear[age_band], `[[`, 0, "p_value"), delta_aic = aic_linear - aic_spline, delta_loglik = loglik_linear - loglik_spline)]
cbh_atomic_csv(as.data.frame(fit), file.path(out, "fit_comparison.csv"))
sel <- fread(file.path(out, "selection_by_age.csv"))
selt <- sel[, .(records = sum(records), deaths = sum(deaths)), by = subset]
sel_r <- sel[subset == "PfPR 1-30%", .(surveys = max(surveys), countries = max(countries), regions = max(regions))]

# National burden 2000-2024 (same IHME all-cause deaths and national PfPR as the primary).
age_in <- fread(file.path(base$out, "annual_comparison/country_age_estimates_2000_2024.csv"))
cty <- fread(file.path(base$out, "annual_comparison/country_estimates_2000_2024.csv"))
bur <- rbindlist(lapply(c("v7", "linear"), function(m) rbindlist(lapply(ages, function(a) {
  z <- age_in[age_band == a]; k <- lhr(comp[[m]][[a]], z$pfpr_pct, 0); k1 <- lhr(comp[[m]][[a]], z$pfpr_pct, pmin(z$pfpr_pct, 1))
  z[, .(model = m, iso3, year, age_band, attributable_deaths = allcause_deaths * (1 - exp(k$est)), attributable_deaths_floor1 = allcause_deaths * (1 - exp(k1$est)))]
}))))
chk <- merge(bur[model == "v7"], age_in[, .(iso3, year, age_band, saved = attributable_deaths)], by = c("iso3", "year", "age_band"))
stopifnot(max(abs(chk$attributable_deaths - chk$saved)) < 1e-6 * max(chk$saved))
cy <- merge(bur[, .(model_deaths = sum(attributable_deaths), model_deaths_floor1 = sum(attributable_deaths_floor1)), by = .(model, iso3, year)],
  cty[, .(iso3, year, country, pfpr_pct, under5_person_years, ihme_malaria_deaths, who_cacode_deaths)], by = c("iso3", "year"))
cy[, pfpr_group := factor(fifelse(pfpr_pct < 1, "<1%", fifelse(pfpr_pct <= 30, "1–30%", ">30%")), levels = c("<1%", "1–30%", ">30%"))]
cbh_atomic_csv(as.data.frame(cy), file.path(out, "country_year_burden.csv"))
yt <- cy[, .(model_deaths = sum(model_deaths), model_deaths_floor1 = sum(model_deaths_floor1), ihme = sum(ihme_malaria_deaths), unigme = sum(who_cacode_deaths), py = sum(under5_person_years)), by = .(model, year)]
yt[, `:=`(rate = 1000 * model_deaths / py, ratio_ihme = model_deaths / ihme, ratio_unigme = model_deaths / unigme)]
setorder(yt, model, year); cbh_atomic_csv(as.data.frame(yt), file.path(out, "year_totals.csv"))
grp <- dcast(cy[year == 2024, .(countries = uniqueN(iso3), deaths = sum(model_deaths), ihme = sum(ihme_malaria_deaths)), by = .(model, pfpr_group)],
  pfpr_group + countries + ihme ~ model, value.var = "deaths")[order(pfpr_group)]
cbh_atomic_csv(as.data.frame(grp), file.path(out, "burden_2024_by_national_pfpr.csv"))
byage <- dcast(bur[year == 2024, .(deaths = sum(attributable_deaths)), by = .(model, age_band)], age_band ~ model, value.var = "deaths")[match(ages, age_band)]

# Figure.
band <- function(a) factor(ifelse(a == "<1", "<1 month", paste(a, "months")), levels = ifelse(ages == "<1", "<1 month", paste(ages, "months")))
curves[, `:=`(band = band(age_band), label = factor(labels[model], levels = labels))]
cols <- setNames(c("#222222", "#2A7AB0", "#C34D26"), labels)
p <- ggplot(curves, aes(pfpr_pct, log_hr, colour = label, fill = label)) +
  annotate("rect", xmin = c(0, 30), xmax = c(1, 40), ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, colour = NA) + geom_line(linewidth = .7) +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40), expand = c(0, 0)) +
  facet_wrap(~band, nrow = 2) + labs(x = expression(italic(Pf)*PR["2–10"]~"(%) at band entry"), y = "Log hazard ratio relative to PfPR 20%") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), panel.spacing.x = unit(1, "lines"))
ggsave(file.path(out, "sfig_pfpr_linear_1_30.png"), p, width = 13, height = 7.5, dpi = 300, device = ragg::agg_png, bg = "white")

f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
ci <- function(h, l, u) sprintf("%.2f (%.2f–%.2f)", h, l, u)
tot <- function(m, y) yt[model == m & year == y]
pc <- function(m, a, b) 100 * (tot(m, b)$rate / tot(m, a)$rate - 1)
writeLines(c("# Supplementary figure: linear PfPR on 1–30%", "",
  sprintf(paste("Association between PfPR[2–10] and all-cause mortality by age band when the log hazard is forced to be linear in PfPR over 1–30%%.",
    "Red: linear term fitted to the %s child-band records (%s deaths) with PfPR 1–30%% at band entry; blue: penalised cubic spline fitted to the same records;",
    "black: the primary spline fitted to all %s records. At 1–23 months the spline is penalised to a straight line and lies under the linear fit. Log hazard ratios relative to PfPR 20%% with 95%% intervals conditional on the smoothing parameters;",
    "shading marks PfPR outside 1–30%%. All models share the primary adjustment (17 covariates, calendar-year spline, survey, country and survey-region random intercepts)."),
    f0(selt[subset == "PfPR 1-30%"]$records), f0(selt[subset == "PfPR 1-30%"]$deaths), f0(selt[subset == "v7 primary"]$records))), file.path(out, "CAPTION.md"))

w <- dcast(con, age_band + contrast ~ model, value.var = c("hr", "lower_95", "upper_95"))
w <- w[order(match(contrast, sprintf("%g%% to %g%%", sapply(steps, `[`, 1), sapply(steps, `[`, 2))), match(age_band, ages))]
sw <- dcast(slope, age_band ~ model, value.var = c("hr_10", "lower_95", "upper_95"))[match(ages, age_band)]
afw <- dcast(af[pfpr_pct == 20], age_band ~ model, value.var = "attributable_fraction")[match(ages, age_band)]
lines <- c("# Linear PfPR[2–10] on the 1–30% range", "",
  sprintf(paste("Records restricted to band entries with PfPR 1–30%%: %s of the %s v7 child-band records (%.0f%%) and %s of %s deaths (%.0f%%), in %d surveys, %d countries and %d survey-regions.",
    "Each age band is fitted twice on these records with the v7 primary specification (17 covariates, `s(calendar_year)` on the primary reference knots, survey, country and survey-region random intercepts, `offset(log(band_years))`, gamma = 2):",
    "**linear**, with `pfpr_pct` as a linear term, and **spline**, with `s(pfpr_pct, cr, k = 5)` on knots equally spaced over 1–30%%.",
    "Fitted by bam (fREML) on the collapsed survey-region × entry-month cells, whose binomial likelihood equals the child-level likelihood; covariates keep the v7 scaling.",
    "The primary curves are the v7 fits to all records."),
    f0(selt[subset == "PfPR 1-30%"]$records), f0(selt[subset == "v7 primary"]$records), 100 * selt[subset == "PfPR 1-30%"]$records / selt[subset == "v7 primary"]$records,
    f0(selt[subset == "PfPR 1-30%"]$deaths), f0(selt[subset == "v7 primary"]$deaths), 100 * selt[subset == "PfPR 1-30%"]$deaths / selt[subset == "v7 primary"]$deaths,
    sel_r$surveys, sel_r$countries, sel_r$regions), "",
  sprintf("All 14 fits converged with positive-definite smoothing Hessians%s.", if (file.exists(file.path(out, "restarts.csv"))) " (strict restarts: see [restarts](restarts.csv))" else ""), "",
  "## Slope: hazard ratio per 10-point increase in PfPR", "",
  "For the splines, the average slope from 1% to 30% (the 30%→1% log hazard ratio × 10/29).", "",
  "| Age (months) | Linear, 1–30% | Spline, 1–30% (average) | Primary v7 (average over 1–30%) |", "|---|---:|---:|---:|",
  sw[, sprintf("| %s | %s | %s | %s |", age_band, ci(hr_10_linear, lower_95_linear, upper_95_linear), ci(hr_10_spline, lower_95_spline, upper_95_spline), ci(hr_10_v7, lower_95_v7, upper_95_v7))], "",
  "## Hazard ratios across the range (95% intervals conditional on smoothing parameters)", "",
  "Equal-width steps show where the splines depart from linearity; the linear model gives the same hazard ratio for every 10-point step.", "",
  "| Age (months) | Contrast | Linear, 1–30% | Spline, 1–30% | Primary v7 |", "|---|---|---:|---:|---:|",
  w[, sprintf("| %s | %s | %s | %s | %s |", age_band, contrast, ci(hr_linear, lower_95_linear, upper_95_linear), ci(hr_spline, lower_95_spline, upper_95_spline), ci(hr_v7, lower_95_v7, upper_95_v7))], "",
  "## Fit on the 1–30% records: linear against spline", "",
  "| Age (months) | Spline PfPR EDF | ΔAIC (linear − spline) | Δ log-likelihood | Linear term p | Spline term p |", "|---|---:|---:|---:|---:|---:|",
  fit[, sprintf("| %s | %.2f | %.1f | %.1f | %s | %s |", age_band, spline_pfpr_edf, delta_aic, delta_loglik, format.pval(linear_p, digits = 2, eps = 1e-4), format.pval(spline_pfpr_p, digits = 2, eps = 1e-4))], "",
  "AIC is mgcv's conditional AIC with the corrected degrees of freedom (lower is better); the log-likelihood is the binomial kernel. A spline EDF near 1 means the penalised spline is itself close to linear.", "",
  "## Burden across the 42 countries", "",
  "The linear model's attributable deaths use 1 − exp(−βP) with national PfPR P, so the log hazard is extrapolated linearly below 1% (to the 0% reference) and above 30%.", "",
  sprintf("Attributable fraction at PfPR 20%% (primary → linear): %s.", paste(sprintf("%s %.1f%% → %.1f%%", afw$age_band, 100 * afw$v7, 100 * afw$linear), collapse = "; ")), "",
  "| Year | Primary v7 | Linear, 1–30% | IHME | UN IGME |", "|---|---:|---:|---:|---:|",
  vapply(c(2000, 2015, 2024), function(y) sprintf("| %d | %s | %s | %s | %s |", y, f0(tot("v7", y)$model_deaths), f0(tot("linear", y)$model_deaths), f0(tot("v7", y)$ihme), f0(tot("v7", y)$unigme)), ""), "",
  sprintf("2024: linear %.2f × IHME and %.2f × UN IGME (primary %.2f and %.2f). Rate change 2000–2015 %.1f%% and 2015–2024 %.1f%% (primary %.1f%% and %.1f%%).",
    tot("linear", 2024)$ratio_ihme, tot("linear", 2024)$ratio_unigme, tot("v7", 2024)$ratio_ihme, tot("v7", 2024)$ratio_unigme,
    pc("linear", 2000, 2015), pc("linear", 2015, 2024), pc("v7", 2000, 2015), pc("v7", 2015, 2024)), "",
  "With a 1% counterfactual floor (deaths attributable to PfPR above 1%, as in `reference_floor_dhsmics_map_gamma2_v2`), which removes the extrapolated 0–1% segment from both models:", "",
  "| Year | Primary v7 | Linear, 1–30% |", "|---|---:|---:|",
  vapply(c(2000, 2015, 2024), function(y) sprintf("| %d | %s | %s |", y, f0(tot("v7", y)$model_deaths_floor1), f0(tot("linear", y)$model_deaths_floor1)), ""), "",
  "2024 deaths by national PfPR of the country:", "",
  "| National PfPR | Countries | Primary v7 | Linear, 1–30% | IHME |", "|---|---:|---:|---:|---:|",
  grp[, sprintf("| %s | %d | %s | %s | %s |", pfpr_group, countries, f0(v7), f0(linear), f0(ihme))], "",
  "2024 deaths by age band:", "", "| Age (months) | Primary v7 | Linear, 1–30% |", "|---|---:|---:|",
  byage[, sprintf("| %s | %s | %s |", age_band, f0(v7), f0(linear))], "",
  "![PfPR curves](sfig_pfpr_linear_1_30.png)", "",
  "Files: `selection_by_age.csv`, `fit_diagnostics.csv`, `fit_comparison.csv`, `smooth_summaries.csv`, `contrasts.csv`, `slope_per_10_points.csv`, `attributable_fraction_by_age.csv`, `pfpr_curves.csv`, `year_totals.csv`, `country_year_burden.csv`, `burden_2024_by_national_pfpr.csv`, `CAPTION.md`. Reproduce: `Rscript R_cbh/sensitivity/linear_pfpr_1_30/01_fit.R` then `02_report.R`.")
writeLines(gsub("\\bv7\\b", pv, lines), file.path(out, "REPORT.md"))
inputs <- c(file.path(base$private, "pfpr_components.rds"), file.path(private, "pfpr_components.rds"),
  file.path(base$out, c("annual_comparison/country_age_estimates_2000_2024.csv", "annual_comparison/country_estimates_2000_2024.csv")),
  file.path(out, c("fit_diagnostics.csv", "smooth_summaries.csv", "selection_by_age.csv")), "R_cbh/sensitivity/linear_pfpr_1_30/02_report.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "report_provenance.csv"))
message("Linear-PfPR report written: ", file.path(out, "REPORT.md"))
