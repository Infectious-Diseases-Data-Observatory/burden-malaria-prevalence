#!/usr/bin/env Rscript
# Report for the tensor-product sensitivity: model fit (tensor and ti versus separate splines),
# PfPR hazard ratios by calendar year, curves, PfPR support by period and the 2000-2024 burden
# with year-specific hazard ratios, and the support/identification checks of 03_checks.R when present.
# Calendar time is evaluated at mid-year, held at the last observed entry month (December 2023,
# which is also the last calendar-year knot) so that 2024 is not extrapolated. Reads aggregates only.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
suppressPackageStartupMessages({library(data.table); library(mgcv); library(ggplot2)})
base <- cbh_primary_settings("regional_mics")
id <- "tensor_pfpr_year_dhsmics_map_gamma2_v1"; out <- file.path("results/cbh", id); private <- file.path("data/derived_cbh/models", id)
ages <- cbh_config()$age_bands$age_band
band_label <- setNames(ifelse(ages == "<1", "<1 month", paste(ages, "months")), ages)
comp <- readRDS(file.path(private, "pfpr_year_components.rds"))
get_comp <- function(m, a) comp[[sprintf("%s_age_%d", m, match(a, ages))]]
models <- c("separate", "tensor", "ti")
model_label <- c(separate = "Separate splines", tensor = "Tensor product te(PfPR, year)", ti = "Separate + ti(PfPR, year)")

# log HR for PfPR `from` -> `to` at calendar time `year` (all vectors recycled), with SE.
lhr <- function(z, from, to, year) {
  n <- max(length(from), length(to), length(year)); from <- rep_len(from, n); to <- rep_len(to, n); year <- rep_len(year, n)
  L <- matrix(0, n, length(z$ix))
  for (k in seq_along(z$smooths)) {
    s <- z$smooths[[k]]; cols <- match(z$index[[k]], z$ix)
    L[, cols] <- L[, cols] + PredictMat(s, data.frame(pfpr_pct = to, calendar_year = year)) - PredictMat(s, data.frame(pfpr_pct = from, calendar_year = year))
  }
  list(est = drop(L %*% z$coef), se = sqrt(pmax(rowSums((L %*% z$Vp) * L), 0)))
}

# Model fit.
dg <- fread(file.path(out, "fit_diagnostics.csv")); ss <- fread(file.path(out, "smooth_summaries.csv"))
fit <- dcast(dg, age_band ~ model, value.var = c("aic", "aic_df", "loglik", "edf_total"))[match(ages, age_band)]
fit[, `:=`(delta_aic_tensor = aic_tensor - aic_separate, delta_aic_ti = aic_ti - aic_separate)]
ti_term <- ss[model == "ti" & grepl("^ti\\(", term), .(age_band, ti_edf = edf, ti_p = `p-value`)]
te_term <- ss[model == "tensor" & grepl("^te\\(", term), .(age_band, te_edf = edf)]
sep_terms <- ss[model == "separate" & term %in% c("s(pfpr_pct)", "s(calendar_year)"), .(sep_edf = sum(edf)), by = age_band]
svy <- dcast(ss[term == "s(survey)"], age_band ~ model, value.var = "edf"); setnames(svy, c("separate", "tensor", "ti"), paste0("survey_edf_", c("separate", "tensor", "ti")))
fit <- Reduce(function(x, y) merge(x, y, by = "age_band"), list(fit, ti_term, te_term, sep_terms, svy))[match(ages, age_band)]
cbh_atomic_csv(as.data.frame(fit), file.path(out, "model_fit_comparison.csv"))

# Hazard ratios by year.
years <- c(2005, 2010, 2015, 2020)
con <- rbindlist(lapply(models, function(m) rbindlist(lapply(ages, function(a) rbindlist(lapply(list(c(40, 20), c(20, 0)), function(x) {
  k <- lhr(get_comp(m, a), x[1], x[2], years + .5)
  data.table(model = m, age_band = a, contrast = sprintf("%g%% to %g%%", x[1], x[2]), year = years, hr = exp(k$est),
             lower_95 = exp(k$est - 1.96 * k$se), upper_95 = exp(k$est + 1.96 * k$se))
}))))))
cbh_atomic_csv(as.data.frame(con), file.path(out, "hazard_ratios_by_year.csv"))
last_entry <- 2023 + 11 / 12; eval_year <- function(y) pmin(y + .5, last_entry)
grid_y <- c(seq(2000.5, 2023.75, by = .25), last_entry)
hr20 <- rbindlist(lapply(models, function(m) rbindlist(lapply(ages, function(a) {
  k <- lhr(get_comp(m, a), 20, 0, grid_y)
  data.table(model = m, age_band = a, calendar_year = grid_y, hr = exp(k$est), lower_95 = exp(k$est - 1.96 * k$se), upper_95 = exp(k$est + 1.96 * k$se))
}))))
cbh_atomic_csv(as.data.frame(hr20), file.path(out, "hr_20_to_0_by_year.csv"))
grid_p <- seq(0, 70, by = .5)
# Curves relative to PfPR 20% in the same year, flagged by whether PfPR lies within that period's
# observed range (children-weighted 2.5th-97.5th percentile of band entries; pfpr_support_by_period.csv).
support <- fread(file.path(out, "pfpr_support_by_period.csv"))
period_of <- c(`2005` = "2005-09", `2010` = "2010-14", `2015` = "2015-19", `2020` = "2020-23")
curves <- rbindlist(lapply(c("separate", "tensor"), function(m) rbindlist(lapply(ages, function(a) rbindlist(lapply(years, function(y) {
  k <- lhr(get_comp(m, a), 20, grid_p, y + .5); sp <- support[age_band == a & period == period_of[as.character(y)]]
  data.table(model = m, age_band = a, year = y, pfpr_pct = grid_p, log_hr = k$est, lower_95 = k$est - 1.96 * k$se, upper_95 = k$est + 1.96 * k$se,
             within_period_support = grid_p >= sp$pfpr_p025 & grid_p <= sp$pfpr_p975)
}))))))
cbh_atomic_csv(as.data.frame(curves), file.path(out, "pfpr_curves_by_year.csv"))

# Burden 2000-2024: each model's year-specific HR (PfPR 0 vs current) at mid-year.
age_in <- fread(file.path(base$out, "annual_comparison/country_age_estimates_2000_2024.csv"))
cty <- fread(file.path(base$out, "annual_comparison/country_estimates_2000_2024.csv"))
bur <- rbindlist(lapply(models, function(m) rbindlist(lapply(ages, function(a) {
  z <- age_in[age_band == a]; k <- lhr(get_comp(m, a), z$pfpr_pct, 0, eval_year(z$year))
  z[, .(model = m, iso3, year, age_band, deaths = allcause_deaths * (1 - exp(k$est)))]
}))))
bur <- rbind(bur, age_in[, .(model = "v7", iso3, year, age_band, deaths = attributable_deaths)])
yt <- merge(bur[, .(deaths = sum(deaths)), by = .(model, year)],
  cty[, .(ihme = sum(ihme_malaria_deaths), unigme = sum(who_cacode_deaths), py = sum(under5_person_years)), by = year], by = "year")
yt[, `:=`(rate = 1000 * deaths / py, ratio_ihme = deaths / ihme)]; setorder(yt, model, year)
cbh_atomic_csv(as.data.frame(yt), file.path(out, "year_totals.csv"))
byage <- bur[year == 2024, .(deaths = sum(deaths)), by = .(model, age_band)]
cbh_atomic_csv(as.data.frame(byage), file.path(out, "deaths_by_age_2024.csv"))
# 2024 under alternative evaluation times for the time-varying models: extrapolated to mid-2024, and
# the HR of mid-2019 (centre of the last well-populated window).
alt <- rbindlist(lapply(c("tensor", "ti"), function(m) rbindlist(lapply(c(`mid-2024 (extrapolated)` = 2024.5, `mid-2019` = 2019.5), function(yy) {
  z <- age_in[year == 2024]; d <- sum(sapply(ages, function(a) { w <- z[age_band == a]; k <- lhr(get_comp(m, a), w$pfpr_pct, 0, yy); sum(w$allcause_deaths * (1 - exp(k$est))) }))
  data.table(model = m, evaluated_at = yy, deaths_2024 = d)
}))))
cbh_atomic_csv(as.data.frame(alt), file.path(out, "burden_2024_alternative_times.csv"))

# Figures.
theme_t <- theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
  panel.border = element_rect(fill = NA, colour = "grey85"))
curves[, `:=`(band = factor(band_label[age_band], levels = band_label), year_lab = factor(paste0("Mid-", year), levels = paste0("Mid-", years)))]
te <- curves[model == "tensor"]; se <- curves[model == "separate"]
# Solid within the period's observed PfPR range, dashed outside (segments split so each is one line).
te[, seg := rleid(within_period_support), by = .(age_band, year)]
p1 <- ggplot() + geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
  geom_ribbon(data = te, aes(pfpr_pct, ymin = lower_95, ymax = upper_95), fill = "#2171b5", alpha = .15) +
  geom_line(data = se, aes(pfpr_pct, log_hr, colour = "Separate splines (same in every year)"), linewidth = .5, linetype = "22") +
  geom_line(data = te, aes(pfpr_pct, log_hr, colour = "Tensor product", linetype = within_period_support, group = interaction(seg, within_period_support)), linewidth = .8) +
  scale_colour_manual(values = c(`Tensor product` = "#2171b5", `Separate splines (same in every year)` = "#222222"), name = NULL) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "31"), breaks = c("TRUE", "FALSE"),
    labels = c(`TRUE` = "Within the period's observed PfPR range", `FALSE` = "Outside it"), name = NULL) +
  facet_grid(band ~ year_lab, scales = "free_y") + coord_cartesian(xlim = c(0, 70)) +
  labs(x = expression(italic(Pf)*PR["2–10"]~"(%) at band entry"), y = "Log hazard ratio relative to PfPR 20% in the same year") +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2, override.aes = list(colour = "#2171b5"))) +
  theme_t + theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(out, "sfig_pfpr_curves_by_year_tensor.png"), p1, width = 12, height = 15, dpi = 300, device = ragg::agg_png, bg = "white")
hr20[, band := factor(band_label[age_band], levels = band_label)]; hr20[, label := factor(model_label[model], levels = model_label)]
p2 <- ggplot(hr20[model != "ti"], aes(calendar_year, hr, colour = label, fill = label)) + geom_hline(yintercept = 1, colour = "grey60", linewidth = .3) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = .15, colour = NA) + geom_line(linewidth = .75) +
  scale_colour_manual(values = c("#222222", "#2171b5"), name = NULL) + scale_fill_manual(values = c("#222222", "#2171b5"), name = NULL) +
  scale_y_log10() + facet_wrap(~band, nrow = 2) + labs(x = "Calendar year of band entry", y = "Hazard ratio, PfPR 0% vs 20% (log scale)") + theme_t
ggsave(file.path(out, "sfig_hr_20_to_0_by_year.png"), p2, width = 13, height = 7.5, dpi = 300, device = ragg::agg_png, bg = "white")

# Report.
f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
cw <- dcast(con[contrast == "20% to 0%" & model %in% c("separate", "tensor")], age_band ~ model + year, value.var = "hr")[match(ages, age_band)]
cw4 <- dcast(con[contrast == "40% to 20%" & model %in% c("separate", "tensor")], age_band ~ model + year, value.var = "hr")[match(ages, age_band)]
row_hr <- function(w) w[, sprintf("| %s | %.2f | %.2f | %.2f | %.2f | %.2f |", age_band, separate_2005, tensor_2005, tensor_2010, tensor_2015, tensor_2020)]
tot <- function(m, y) yt[model == m & year == y]$deaths; pc <- function(m, a, b) 100 * (yt[model == m & year == b]$rate / yt[model == m & year == a]$rate - 1)
ba <- dcast(byage, age_band ~ model, value.var = "deaths")[match(ages, age_band)]
has_checks <- file.exists(file.path(out, "interaction_checks.csv"))
checks_block <- character()
if (has_checks) {
  ic <- fread(file.path(out, "interaction_checks.csv")); bn <- fread(file.path(out, "binned_pfpr_by_period.csv")); sp <- fread(file.path(out, "support_pfpr_class_by_period.csv"))
  low <- dcast(sp[pclass %in% c("<2", "2-5")][, .(d = sum(deaths)), by = .(age_band, period)], age_band ~ period, value.var = "d")
  b5 <- dcast(bn[contrast == "<5 vs 20-40%"], age_band ~ period, value.var = c("log_hr", "deaths_in_class"))
  checks_block <- c("## Support and identification checks (24–59 months; [03_checks.R](../../../R_cbh/sensitivity/tensor_pfpr_year/03_checks.R))", "",
    "Deaths in cells with PfPR below 5%, by period of band entry:", "", paste0("| Age (months) | ", paste(names(low)[-1], collapse = " | "), " |"), paste0("|---|", strrep("---:|", ncol(low) - 1)),
    apply(low, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "",
    "Binned adjusted log hazard ratio for PfPR <5% versus 20–40% by period (survey and country random intercepts, calendar spline, 17 covariates), with the deaths in the <5% class:", "",
    paste0("| Age (months) | ", paste(sort(unique(bn$period)), collapse = " | "), " |"), paste0("|---|", strrep("---:|", uniqueN(bn$period))),
    sapply(unique(bn$age_band), function(a) paste0("| ", a, " | ", paste(bn[age_band == a & contrast == "<5 vs 20-40%"][order(period), sprintf("%.2f (%d)", log_hr, deaths_in_class)], collapse = " | "), " |")), "",
    "Interaction tests (ΔAIC against the matching model without the interaction; HR for PfPR 0% vs 20% at mid-2005 and mid-2020):", "",
    "| Age (months) | Check | ΔAIC | ti EDF | ti p | survey RE EDF | HR 2005 | HR 2020 |", "|---|---|---:|---:|---:|---:|---:|---:|",
    ic[, sprintf("| %s | %s | %.1f | %s | %s | %s | %.2f | %.2f |", age_band, check, delta_aic, ifelse(is.na(ti_p), "–", sprintf("%.2f", ti_edf)), ifelse(is.na(ti_p), "–", format.pval(ti_p, digits = 2, eps = 1e-4)),
      ifelse(survey_edf == 0, "(fixed effects)", sprintf("%.1f", survey_edf)), hr_2005, hr_2020)], "",
    "Full tables: [support_pfpr_class_by_period.csv](support_pfpr_class_by_period.csv), [binned_pfpr_by_period.csv](binned_pfpr_by_period.csv), [interaction_checks.csv](interaction_checks.csv).", "")
}
alt <- fread(file.path(out, "burden_2024_alternative_times.csv"))
lines <- c("# Tensor product of PfPR and calendar time (no region random intercept)", "",
  "Seven age-band models without the survey-region random intercept, fitted with bam (fREML, no discretisation) to the collapsed survey-region × entry-month cells (binomial deaths out of children entering the band; the same likelihood as the child-level model). Three specifications per band, otherwise identical (17 covariates, survey and country random intercepts, band-width offset, primary reference knots, gamma = 2):", "",
  "- **Separate splines:** `s(pfpr_pct, cr, k=5) + s(calendar_year, cr, k=6)` (as the primary).",
  "- **Tensor product:** `te(pfpr_pct, calendar_year, bs = c(\"cr\",\"cr\"), k = c(5, 6))`, which lets the PfPR curve change with calendar time.",
  "- **Interaction test:** separate splines plus `ti(pfpr_pct, calendar_year)`, the pure interaction. The function space contains the separate-spline model, but the penalised fits are not nested: smoothing parameters, including the survey random-effect variance, are re-estimated.", "",
  if (file.exists(file.path(out, "restarts.csv"))) "Strict restarts: see [restarts.csv](restarts.csv)." else "All 21 fits converged with positive-definite smoothing Hessians; none needed a restart.", "",
  "## Conclusion", "",
  "Keep the separate-spline (time-constant) PfPR curves as the primary specification and report the tensor product as an exploratory sensitivity. Allowing the PfPR effect to change with calendar time improves AIC in the full-period fits only at <1 and 1–5 months, where the curves barely change. At 24–59 months the fitted low-PfPR gradient steepens over time, and the binned data show the same direction, but the evidence is mixed: in the full-period fits AIC does not favour the interaction and the log-likelihood does not improve, because the interaction takes over between-survey heterogeneity (the survey random-effect EDF falls); it is favoured when the random-effect variances are held at the separate-spline values and at 36–47 months within 2005–2019; it disappears when identified within surveys (survey fixed effects); and it rests on a few dozen deaths in PfPR <5% areas at each end of the period. It therefore cannot be separated from between-survey (country-period) differences. If the steepening were real it would raise the 2024 attributable deaths and shrink the post-2015 decline (below).", "",
  "## Model fit", "",
  "| Age (months) | ΔAIC tensor − separate | ΔAIC ti − separate | ti EDF | ti p-value | survey RE EDF: separate / tensor / ti | EDF of PfPR/time terms: separate, tensor |", "|---|---:|---:|---:|---:|---:|---:|",
  fit[, sprintf("| %s | %.1f | %.1f | %.2f | %s | %.1f / %.1f / %.1f | %.1f, %.1f |", age_band, delta_aic_tensor, delta_aic_ti, ti_edf, format.pval(ti_p, digits = 2, eps = 1e-4),
    survey_edf_separate, survey_edf_tensor, survey_edf_ti, sep_edf, te_edf)], "",
  "AIC is mgcv's (binomial likelihood with the smoothing-corrected degrees of freedom, `aic_df` in [fit_diagnostics.csv](fit_diagnostics.csv); the `loglik` column there is the kernel log-likelihood without the binomial coefficient). Negative ΔAIC favours the more flexible model; |ΔAIC| below about 3 (6–11, 36–47 and 48–59 months) is no discernible difference, and its sign depends on the degrees-of-freedom definition. The ti p-value is mgcv's approximate test conditional on the estimated smoothing parameters; in the older bands it coexists with no gain in log-likelihood because the survey random effect is shrunk when the interaction enters, so it should not be read as stand-alone evidence of effect modification.", "",
  "## Hazard ratio for PfPR 20% → 0% by calendar year of band entry (evaluated at mid-year)", "",
  "| Age (months) | Separate (all years) | Tensor mid-2005 | mid-2010 | mid-2015 | mid-2020 |", "|---|---:|---:|---:|---:|---:|", row_hr(cw), "",
  "## Hazard ratio for PfPR 40% → 20% by calendar year (mid-year)", "", "| Age (months) | Separate (all years) | Tensor mid-2005 | mid-2010 | mid-2015 | mid-2020 |", "|---|---:|---:|---:|---:|---:|", row_hr(cw4), "",
  checks_block,
  "## Burden across the 42 countries (IHME all-cause inputs, year-specific hazard ratios)", "",
  "Calendar time is evaluated at mid-year, held at December 2023 (the last observed band entry and last calendar-year knot) for 2024.", "",
  "| Model | 2000 | 2015 | 2024 | 2024 × IHME | Rate change 2000–2015 | Rate change 2015–2024 |", "|---|---:|---:|---:|---:|---:|---:|",
  sapply(c("v7", "separate", "tensor", "ti"), function(m) sprintf("| %s | %s | %s | %s | %.2f | %.1f%% | %.1f%% |",
    c(v7 = "Primary v7 (with region RE, child-level)", separate = "Separate splines", tensor = "Tensor product", ti = "Separate + ti (see note)")[m], f0(tot(m, 2000)), f0(tot(m, 2015)), f0(tot(m, 2024)), yt[model == m & year == 2024]$ratio_ihme, pc(m, 2000, 2015), pc(m, 2015, 2024))), "",
  sprintf("2024 under other evaluation times: tensor %s at mid-2024 (extrapolated) and %s with the mid-2019 hazard ratios; ti %s and %s.",
    f0(alt[model == "tensor" & evaluated_at == 2024.5]$deaths_2024), f0(alt[model == "tensor" & evaluated_at == 2019.5]$deaths_2024),
    f0(alt[model == "ti" & evaluated_at == 2024.5]$deaths_2024), f0(alt[model == "ti" & evaluated_at == 2019.5]$deaths_2024)), "",
  "2024 attributable deaths by age band:", "", "| Age (months) | Separate | Tensor | ti |", "|---|---:|---:|---:|", ba[, sprintf("| %s | %s | %s | %s |", age_band, f0(separate), f0(tensor), f0(ti))], "",
  "Note: the ti model's neonatal interaction (about 4.6 EDF) oscillates over time around a hazard ratio close to 1; applied to the large neonatal all-cause totals this moves tens of thousands of deaths, so the ti burden row is unstable and is shown only for completeness. The tensor product, which spans essentially the same function space, does not show this. The time-varying models' lower 2000 and higher 2024 totals, and hence their smaller declines, follow directly from the assumed steepening of the PfPR effect and are not separate findings.", "",
  "PfPR curves of the tensor product by year of band entry (relative to PfPR 20% in the same year, 95% intervals conditional on smoothing parameters; solid within the period's children-weighted 2.5th–97.5th PfPR percentiles, dashed outside), with the time-constant separate-spline curve for reference:", "", "![PfPR curves by year](sfig_pfpr_curves_by_year_tensor.png)", "", "![HR 20% to 0% by year](sfig_hr_20_to_0_by_year.png)", "",
  "Files: `model_fit_comparison.csv`, `hazard_ratios_by_year.csv`, `hr_20_to_0_by_year.csv`, `pfpr_curves_by_year.csv`, `pfpr_support_by_period.csv`, `year_totals.csv`, `deaths_by_age_2024.csv`, `burden_2024_alternative_times.csv`, `fit_diagnostics.csv`, `smooth_summaries.csv`, `model_formulas.txt`. Reproduce: `Rscript R_cbh/sensitivity/tensor_pfpr_year/01_fit.R`, then `03_checks.R`, then `02_report.R`.")
writeLines(lines, file.path(out, "REPORT.md"))
inputs <- c(file.path(private, "pfpr_year_components.rds"), file.path(out, c("fit_diagnostics.csv", "smooth_summaries.csv", "pfpr_support_by_period.csv")),
  if (has_checks) file.path(out, c("interaction_checks.csv", "binned_pfpr_by_period.csv", "support_pfpr_class_by_period.csv")),
  file.path(base$out, c("annual_comparison/country_age_estimates_2000_2024.csv", "annual_comparison/country_estimates_2000_2024.csv")), "R_cbh/sensitivity/tensor_pfpr_year/02_report.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "report_provenance.csv"))
print(fit[, .(age_band, delta_aic_tensor, delta_aic_ti, ti_edf, ti_p)]); message("Tensor-product report written: ", file.path(out, "REPORT.md"))
