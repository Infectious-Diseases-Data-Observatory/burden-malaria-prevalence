#!/usr/bin/env Rscript
# Comparison of the primary without the survey-region random intercept with the primary (written for v7, now v9):
# PfPR hazard ratios, attributable fractions, curves, random-effect EDFs, AIC and the national
# burden 2000-2024 against IHME and UN IGME. Reads saved fits and aggregates only.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/sensitivity/no_region_re/settings.R")
suppressPackageStartupMessages({library(data.table); library(mgcv); library(ggplot2)})
st <- cbh_no_region_settings(); base <- st$base; out <- st$out
ages <- cbh_config()$age_bands$age_band
# Display label of the fitted primary (these reports were first written against v7).
pv <- sub(".*_(v[0-9]+)$", "\\1", base$id)
labels <- c(v7 = sprintf("Primary %s (survey, country and region random intercepts)", pv), no_region = "Without the region random intercept")
comp <- list(v7 = readRDS(file.path(base$private, "pfpr_components.rds")), no_region = readRDS(file.path(st$private, "pfpr_components.rds")))
names(comp$v7) <- ages[as.integer(sub(".*_", "", names(comp$v7)))]; names(comp$no_region) <- ages[as.integer(sub(".*_", "", names(comp$no_region)))]
stopifnot(setequal(names(comp$v7), ages), setequal(names(comp$no_region), ages))
lhr <- function(piece, from, to) {
  n <- max(length(from), length(to)); from <- rep_len(from, n); to <- rep_len(to, n)
  L <- PredictMat(piece$smooth, data.frame(pfpr_pct = to)) - PredictMat(piece$smooth, data.frame(pfpr_pct = from))
  list(est = drop(L %*% piece$coef), se = sqrt(pmax(rowSums((L %*% piece$covariance) * L), 0)))
}

# Contrasts and attributable fractions.
con <- rbindlist(lapply(names(comp), function(m) rbindlist(lapply(ages, function(a) rbindlist(lapply(list(c(40, 20), c(20, 0)), function(x) {
  k <- lhr(comp[[m]][[a]], x[1], x[2])
  data.table(model = m, age_band = a, contrast = sprintf("%g%% to %g%%", x[1], x[2]), hr = exp(k$est), lower_95 = exp(k$est - 1.96 * k$se), upper_95 = exp(k$est + 1.96 * k$se))
}))))))
cbh_atomic_csv(as.data.frame(con), file.path(out, "comparison_contrasts.csv"))
af <- rbindlist(lapply(names(comp), function(m) rbindlist(lapply(ages, function(a) {
  k <- lhr(comp[[m]][[a]], c(10, 20, 30, 40), 0)
  data.table(model = m, age_band = a, pfpr_pct = c(10, 20, 30, 40), attributable_fraction = 1 - exp(k$est))
}))))
cbh_atomic_csv(as.data.frame(af), file.path(out, "attributable_fraction_by_age.csv"))
grid <- seq(0, 70, by = .5)
curves <- rbindlist(lapply(names(comp), function(m) rbindlist(lapply(ages, function(a) {
  k <- lhr(comp[[m]][[a]], rep(20, length(grid)), grid)
  data.table(model = m, age_band = a, pfpr_pct = grid, log_hr = k$est, lower = k$est - 1.96 * k$se, upper = k$est + 1.96 * k$se)
}))))
cbh_atomic_csv(as.data.frame(curves), file.path(out, "comparison_pfpr_curves.csv"))

# Random-effect EDFs, fREML and AIC.
edf <- rbind(cbind(model = "v7", fread(file.path(base$out, "smooth_summaries.csv"))[, .(age_band, term, edf)]),
             cbind(model = "no_region", fread(file.path(out, "smooth_summaries.csv"))[, .(age_band, term, edf)]))
cbh_atomic_csv(as.data.frame(edf), file.path(out, "comparison_edf.csv"))
v7m <- fread(file.path(base$out, "fit_manifest.csv"))[series == "map_full"][match(ages, age_band)]
v7_fit_stats <- rbindlist(lapply(seq_len(nrow(v7m)), function(i) {
  stopifnot(identical(cbh_file_hash(v7m$model_file[i]), v7m$md5[i]))
  f <- readRDS(v7m$model_file[i])$fit; z <- data.table(age_band = v7m$age_band[i], aic = AIC(f), fREML = unname(f$gcv.ubre), edf_total = sum(f$edf)); rm(f); gc(FALSE); z
}))
nr <- fread(file.path(out, "fit_diagnostics.csv"))[, .(age_band, aic, fREML, edf_total, converged)]
fitcmp <- merge(v7_fit_stats, nr, by = "age_band", suffixes = c("_v7", "_no_region"))[match(ages, age_band)]
fitcmp[, `:=`(delta_aic = aic_no_region - aic_v7, delta_fREML = fREML_no_region - fREML_v7)]
cbh_atomic_csv(as.data.frame(fitcmp), file.path(out, "comparison_fit_statistics.csv"))

# National burden 2000-2024 with each model's curves (same IHME inputs and national PfPR).
age_in <- fread(file.path(base$out, "annual_comparison/country_age_estimates_2000_2024.csv"))
cty <- fread(file.path(base$out, "annual_comparison/country_estimates_2000_2024.csv"))
bur <- rbindlist(lapply(names(comp), function(m) rbindlist(lapply(ages, function(a) {
  z <- age_in[age_band == a]; k <- lhr(comp[[m]][[a]], z$pfpr_pct, 0)
  z[, .(model = m, iso3, year, age_band, attributable_deaths = allcause_deaths * (1 - exp(k$est)))]
}))))
chk <- merge(bur[model == "v7"], age_in[, .(iso3, year, age_band, saved = attributable_deaths)], by = c("iso3", "year", "age_band"))
stopifnot(max(abs(chk$attributable_deaths - chk$saved)) < 1e-6 * max(chk$saved))
cy <- merge(bur[, .(model_deaths = sum(attributable_deaths)), by = .(model, iso3, year)],
  cty[, .(iso3, year, country, under5_person_years, ihme_malaria_deaths, who_cacode_deaths)], by = c("iso3", "year"))
cbh_atomic_csv(as.data.frame(cy), file.path(out, "country_year_burden.csv"))
yt <- cy[, .(model_deaths = sum(model_deaths), ihme = sum(ihme_malaria_deaths), unigme = sum(who_cacode_deaths), py = sum(under5_person_years)), by = .(model, year)]
yt[, `:=`(rate = 1000 * model_deaths / py, ratio_ihme = model_deaths / ihme, ratio_unigme = model_deaths / unigme)]
setorder(yt, model, year); cbh_atomic_csv(as.data.frame(yt), file.path(out, "year_totals.csv"))
byage <- bur[year == 2024, .(deaths = sum(attributable_deaths)), by = .(model, age_band)]

# Figure: PfPR curves, v7 against no region random intercept.
curves[, `:=`(band = factor(ifelse(age_band == "<1", "<1 month", paste(age_band, "months")), levels = ifelse(ages == "<1", "<1 month", paste(ages, "months"))),
              label = factor(labels[model], levels = labels))]
cols <- setNames(c("#222222", "#C34D26"), labels)
p <- ggplot(curves, aes(pfpr_pct, log_hr, colour = label, fill = label)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, colour = NA) + geom_line(linewidth = .7) +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
  facet_wrap(~band, nrow = 2) + labs(x = expression(italic(Pf)*PR["2–10"]~"(%) at band entry"), y = "Log hazard ratio relative to PfPR 20%") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(out, "sfig_pfpr_splines_no_region_re.png"), p, width = 13, height = 7.5, dpi = 300, device = ragg::agg_png, bg = "white")

f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
w <- dcast(con, age_band + contrast ~ model, value.var = c("hr", "lower_95", "upper_95"))
w <- w[order(contrast != "40% to 20%", match(age_band, ages))]
ci <- function(h, l, u) sprintf("%.2f (%.2f–%.2f)", h, l, u)
t24 <- yt[year == 2024]; t15 <- yt[year == 2015]; t00 <- yt[year == 2000]
pc <- function(m, a, b) 100 * (yt[model == m & year == b]$rate / yt[model == m & year == a]$rate - 1)
e <- dcast(edf[term %in% c("s(pfpr_pct)", "s(calendar_year)", "s(survey)", "s(country)", "s(region)")], age_band + term ~ model, value.var = "edf")
e <- e[order(match(term, c("s(pfpr_pct)", "s(calendar_year)", "s(survey)", "s(country)", "s(region)")), match(age_band, ages))]
afw <- dcast(af[pfpr_pct == 20], age_band ~ model, value.var = "attributable_fraction")[match(ages, age_band)]
lines <- c("# Primary model without the survey-region random intercept", "",
  sprintf("Same prepared v7 sample (%s child-band records, %s deaths), formula, reference knots, gamma = 2 and bam settings as `%s`, with `s(region, bs = \"re\")` removed; survey and country random intercepts kept. All seven fits converged; the 1-5 month fit first stopped with a non-positive-definite smoothing Hessian (PfPR smooth penalised to a straight line, fREML 1.3 worse) and was replaced by the strict restart under the primary rule ([restarts](restarts.csv)).",
    f0(base$expected_records), f0(base$expected_deaths), base$id), "",
  "## Hazard ratios (95% intervals conditional on smoothing parameters)", "",
  "| Age (months) | Contrast | v7 | Without region RE |", "|---|---|---:|---:|",
  w[, sprintf("| %s | %s | %s | %s |", age_band, contrast, ci(hr_v7, lower_95_v7, upper_95_v7), ci(hr_no_region, lower_95_no_region, upper_95_no_region))], "",
  sprintf("Largest absolute change: %.3f in the 40%%→20%% hazard ratio and %.3f in the 20%%→0%% hazard ratio.",
    max(abs(w[contrast == "40% to 20%"]$hr_no_region - w[contrast == "40% to 20%"]$hr_v7)), max(abs(w[contrast == "20% to 0%"]$hr_no_region - w[contrast == "20% to 0%"]$hr_v7))), "",
  "Attributable fraction at PfPR 20% (v7 → without region RE): ", paste(sprintf("%s %.1f%% → %.1f%%", afw$age_band, 100 * afw$v7, 100 * afw$no_region), collapse = "; "), "",
  "## Burden across the 42 countries", "",
  "| Year | v7 | Without region RE | IHME | UN IGME |", "|---|---:|---:|---:|---:|",
  sprintf("| %d | %s | %s | %s | %s |", c(2000, 2015, 2024), f0(c(t00[model == "v7"]$model_deaths, t15[model == "v7"]$model_deaths, t24[model == "v7"]$model_deaths)),
    f0(c(t00[model == "no_region"]$model_deaths, t15[model == "no_region"]$model_deaths, t24[model == "no_region"]$model_deaths)),
    f0(c(t00$ihme[1], t15$ihme[1], t24$ihme[1])), f0(c(t00$unigme[1], t15$unigme[1], t24$unigme[1]))), "",
  sprintf("2024: %.2f × IHME and %.2f × UN IGME without the region random intercept (v7 %.2f and %.2f). Rate change 2000–2015 %.1f%% and 2015–2024 %.1f%% (v7 %.1f%% and %.1f%%).",
    t24[model == "no_region"]$ratio_ihme, t24[model == "no_region"]$ratio_unigme, t24[model == "v7"]$ratio_ihme, t24[model == "v7"]$ratio_unigme,
    pc("no_region", 2000, 2015), pc("no_region", 2015, 2024), pc("v7", 2000, 2015), pc("v7", 2015, 2024)), "",
  "## Effective degrees of freedom", "", "| Term | Age (months) | v7 | Without region RE |", "|---|---|---:|---:|",
  e[, sprintf("| %s | %s | %.2f | %s |", term, age_band, v7, ifelse(is.na(no_region), "–", sprintf("%.2f", no_region)))], "",
  "## Fit statistics", "", "| Age (months) | ΔAIC (without − v7) | ΔfREML (without − v7) |", "|---|---:|---:|",
  fitcmp[, sprintf("| %s | %.1f | %.2f |", age_band, delta_aic, delta_fREML)], "",
  "AIC is mgcv's conditional AIC; fREML is the restricted likelihood criterion minimised by bam (lower is better for both). Both models share the fixed effects, so the fREML difference compares the random-effect structures directly.", "",
  "![PfPR curves](sfig_pfpr_splines_no_region_re.png)", "",
  "Files: `comparison_contrasts.csv`, `attributable_fraction_by_age.csv`, `comparison_pfpr_curves.csv`, `comparison_edf.csv`, `comparison_fit_statistics.csv`, `year_totals.csv`, `country_year_burden.csv`. Reproduce: `Rscript R_cbh/sensitivity/no_region_re/01_fit.R` then `02_report.R`.")
writeLines(gsub("\\bv7\\b", pv, lines), file.path(out, "REPORT.md"))
inputs <- c(file.path(base$private, "pfpr_components.rds"), file.path(st$private, "pfpr_components.rds"), file.path(base$out, c("smooth_summaries.csv", "fit_manifest.csv",
  "annual_comparison/country_age_estimates_2000_2024.csv", "annual_comparison/country_estimates_2000_2024.csv")), file.path(out, c("smooth_summaries.csv", "fit_diagnostics.csv")),
  "R_cbh/sensitivity/no_region_re/02_report.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "report_provenance.csv"))
print(byage); message("No-region-RE report written: ", file.path(out, "REPORT.md"))
