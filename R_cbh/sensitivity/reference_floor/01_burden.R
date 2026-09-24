#!/usr/bin/env Rscript
# Burden with a prevalence floor (added 24 September 2026). The primary attributes to malaria
# the fall in each band's hazard from the current national PfPR to PfPR = 0, which lies below
# the observed exposure support (2.5th percentile about 1.7%). Here the counterfactual is
# min(P, floor): deaths attributable to transmission above the floor. floor = 1% is the
# headline; 0% (the primary, reproduced exactly), 2% and 5% are shown alongside.
# Uses the saved v7 spline components and annual inputs only: no refitting, no downloads.
# Deaths below the floor are not estimated here (the WMR case-fatality method would need a
# case source); the report states this.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
st <- cbh_primary_settings("regional_mics")
out <- "results/cbh/reference_floor_dhsmics_map_gamma2_v1"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
floors <- c(0, 1, 2, 5); headline <- 1
paths <- c(components = file.path(st$private, "pfpr_components.rds"),
  age = file.path(st$out, "annual_comparison/country_age_estimates_2000_2024.csv"),
  country = file.path(st$out, "annual_comparison/country_estimates_2000_2024.csv"),
  manifest = file.path(st$out, "fit_manifest.csv"),
  ihme_inputs_2024 = "results/cbh/age_band_separate_v1/country_burden_2024/ihme_disjoint_age_inputs.csv")
stopifnot(all(file.exists(paths)))
components <- readRDS(paths[["components"]]); manifest <- fread(paths[["manifest"]])
for (piece in components) stopifnot(identical(piece$model_md5, manifest$md5[match(piece$fit_id, manifest$fit_id)]))
ages <- cbh_config()$age_bands$age_band
age <- fread(paths[["age"]]); country <- fread(paths[["country"]])
stopifnot(setequal(age$age_band, ages), all(age[, .N, by = .(iso3, year)]$N == 7L))

# log HR for moving from P to min(P, floor), exact from each band's spline, with its SE.
contrast <- function(band, p, floor) {
  piece <- components[[paste0("map_full_age_", match(band, ages))]]
  L <- PredictMat(piece$smooth, data.frame(pfpr_pct = pmin(p, floor))) - PredictMat(piece$smooth, data.frame(pfpr_pct = p))
  list(est = drop(L %*% piece$coef), se = sqrt(pmax(rowSums((L %*% piece$covariance) * L), 0)))
}
rows <- rbindlist(lapply(floors, function(fl) rbindlist(lapply(ages, function(g) {
  z <- age[age_band == g]; k <- contrast(g, z$pfpr_pct, fl)
  z[, .(iso3, year, age_band, floor_pct = fl, pfpr_pct, allcause_deaths, age_person_years,
        log_hr = k$est, log_hr_se = k$se, attributable_deaths = allcause_deaths * (1 - exp(k$est)))]
}))))
# The 0% floor must reproduce the saved primary attributable deaths.
chk <- merge(rows[floor_pct == 0], age[, .(iso3, year, age_band, saved = attributable_deaths)], by = c("iso3", "year", "age_band"))
stopifnot(max(abs(chk$attributable_deaths - chk$saved)) < 1e-6 * max(chk$saved))
cbh_atomic_csv(as.data.frame(rows), file.path(out, "country_year_age.csv"))

cy <- rows[, .(attributable_deaths = sum(attributable_deaths), allcause_deaths = sum(allcause_deaths)), by = .(iso3, year, floor_pct)]
cy <- merge(cy, country[, .(iso3, year, country, pfpr_pct, under5_person_years, ihme_malaria_deaths, who_cacode_deaths)], by = c("iso3", "year"))
stopifnot(nrow(cy) == 42L * 25L * length(floors))
cbh_atomic_csv(as.data.frame(cy), file.path(out, "country_year.csv"))

yt <- cy[, .(model_deaths = sum(attributable_deaths), ihme_malaria_deaths = sum(ihme_malaria_deaths),
             unigme_cacode_deaths = sum(who_cacode_deaths), person_years = sum(under5_person_years)), by = .(floor_pct, year)]
yt[, `:=`(model_rate_per1000 = 1000 * model_deaths / person_years, ihme_rate_per1000 = 1000 * ihme_malaria_deaths / person_years,
          unigme_rate_per1000 = 1000 * unigme_cacode_deaths / person_years,
          ratio_ihme = model_deaths / ihme_malaria_deaths, ratio_unigme = model_deaths / unigme_cacode_deaths)]
setorder(yt, floor_pct, year)
cbh_atomic_csv(as.data.frame(yt), file.path(out, "year_totals.csv"))

pct <- function(a, b) 100 * (b / a - 1)
summ <- yt[year %in% c(2000, 2015, 2024)]
chg <- yt[, .(deaths_2000_2024 = pct(model_deaths[year == 2000], model_deaths[year == 2024]),
              deaths_2015_2024 = pct(model_deaths[year == 2015], model_deaths[year == 2024]),
              rate_2000_2015 = pct(model_rate_per1000[year == 2000], model_rate_per1000[year == 2015]),
              rate_2015_2024 = pct(model_rate_per1000[year == 2015], model_rate_per1000[year == 2024]),
              rate_2000_2024 = pct(model_rate_per1000[year == 2000], model_rate_per1000[year == 2024])), by = floor_pct]
cbh_atomic_csv(as.data.frame(chg), file.path(out, "trend_changes.csv"))

# 2024 country comparison (Table 1 analogue) and agreement statistics under each floor.
c24 <- cy[year == 2024]
c24[, `:=`(model_rate = 1000 * attributable_deaths / under5_person_years, ihme_rate = 1000 * ihme_malaria_deaths / under5_person_years)]
agree <- c24[, .(countries_above_ihme = sum(attributable_deaths > ihme_malaria_deaths), median_ratio_ihme = median(attributable_deaths / ihme_malaria_deaths),
                 spearman_rate_ihme = cor(model_rate, ihme_rate, method = "spearman")), by = floor_pct]
cbh_atomic_csv(as.data.frame(agree), file.path(out, "agreement_2024.csv"))
top <- c24[floor_pct == headline][order(-(attributable_deaths - ihme_malaria_deaths))]
cbh_atomic_csv(as.data.frame(top[, .(country, pfpr_pct, model_deaths = attributable_deaths, ihme_malaria_deaths, unigme_cacode_deaths = who_cacode_deaths)]),
  file.path(out, "country_comparison_2024_floor1.csv"))
byage <- rows[year == 2024, .(attributable_deaths = sum(attributable_deaths)), by = .(floor_pct, age_band)]
cbh_atomic_csv(as.data.frame(byage), file.path(out, "deaths_by_age_2024.csv"))

# Figure 2 analogue: band attributable fraction at PfPR 10/20/30/40% against the floor.
af <- rbindlist(lapply(floors, function(fl) rbindlist(lapply(ages, function(g) {
  p <- c(10, 20, 30, 40); k <- contrast(g, p, fl)
  data.table(floor_pct = fl, age_band = g, pfpr_pct = p, attributable_fraction = 1 - exp(k$est),
             af_lower_95 = 1 - exp(k$est + 1.96 * k$se), af_upper_95 = 1 - exp(k$est - 1.96 * k$se))
}))))
cbh_atomic_csv(as.data.frame(af), file.path(out, "attributable_fraction_by_age.csv"))

# Figure 4 analogue: malaria-caused probability of dying before 5 in 2024 against the floor
# (same life-table construction as R_cbh/reporting/13_under5_death_probability.R).
edges <- c(0, 28 / 365.25 * 12, 6, 12, 24, 36, 48, 60); src <- fread(paths[["ihme_inputs_2024"]])
q5 <- rbindlist(lapply(floors, function(fl) rbindlist(lapply(sort(unique(age$iso3)), function(iso) {
  z <- rows[floor_pct == fl & year == 2024 & iso3 == iso][match(ages, age_band)]
  H <- z$allcause_deaths / z$age_person_years * diff(edges) / 12
  s <- src[iso3 == iso]
  H[1] <- s[source_age == "0 to 6 days (early neonatal)"]$rate_per100000 / 1e5 * 7 / 365.25 + s[source_age == "7 to 27 days (late neonatal)"]$rate_per100000 / 1e5 * 21 / 365.25
  data.table(floor_pct = fl, iso3 = iso, q5_allcause = -expm1(-sum(H)), q5_malaria = -expm1(-sum(H)) + expm1(-sum(H * exp(z$log_hr))))
}))))
cbh_atomic_csv(as.data.frame(q5), file.path(out, "under5_probability_2024.csv"))

f <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")
h <- yt[floor_pct == headline]; p0 <- yt[floor_pct == 0]
a1 <- agree[floor_pct == headline]; c1 <- chg[floor_pct == headline]; c0 <- chg[floor_pct == 0]
ihme_rate <- function(y0, y1) pct(h$ihme_rate_per1000[h$year == y0], h$ihme_rate_per1000[h$year == y1])
ig_rate <- function(y0, y1) pct(h$unigme_rate_per1000[h$year == y0], h$unigme_rate_per1000[h$year == y1])
tab <- function(y) sprintf("| %d | %s (%s) | %s | %s | %.2f (%.2f) | %.2f | %.2f |", y, f(h$model_deaths[h$year == y]), f(p0$model_deaths[p0$year == y]),
  f(h$ihme_malaria_deaths[h$year == y]), f(h$unigme_cacode_deaths[h$year == y]), h$model_rate_per1000[h$year == y], p0$model_rate_per1000[p0$year == y],
  h$ihme_rate_per1000[h$year == y], h$unigme_rate_per1000[h$year == y])
fl_rows <- yt[year == 2024, sprintf("| %g%% | %s | %.2f | %.2f | %.2f |", floor_pct, f(model_deaths), ratio_ihme, ratio_unigme, model_rate_per1000)]
af1 <- dcast(af[floor_pct == headline], age_band ~ pfpr_pct, value.var = "attributable_fraction")[match(ages, age_band)]
af0 <- dcast(af[floor_pct == 0], age_band ~ pfpr_pct, value.var = "attributable_fraction")[match(ages, age_band)]
by1 <- byage[floor_pct == headline][match(ages, age_band)]; by0 <- byage[floor_pct == 0][match(ages, age_band)]
q1 <- q5[floor_pct == headline]; q0 <- q5[floor_pct == 0]
lines <- c("# Burden with a 1% prevalence floor", "",
  "Counterfactual PfPR[2-10] = min(current, 1%) instead of 0: deaths attributable to *P. falciparum* transmission above 1%, from the saved v7 spline components (no refitting). The 0% floor reproduces the primary exactly. Point estimates only, as in the primary. Deaths occurring at or below 1% transmission are not included; see the note at the end.", "",
  "## Totals across the 42 countries (primary 0% floor in brackets)", "",
  "| Year | PfPR-ACM, 1% floor (0%) | IHME | UN IGME | Rate per 1,000 child-years, 1% floor (0%) | IHME rate | UN IGME rate |", "|---|---:|---:|---:|---:|---:|---:|",
  tab(2000), tab(2015), tab(2024), "",
  sprintf("With the 1%% floor the 2024 total is %s, %.0f%% above IHME and %.0f%% above UN IGME (primary: %.0f%% and %.0f%%). The rate falls %.1f%% from 2000 to 2015 and %.1f%% from 2015 to 2024 (primary %.1f%% and %.1f%%; IHME %.1f%% and %.1f%%; UN IGME %.1f%% and %.1f%%). From 2000 to 2024 the rate falls %.1f%% (primary %.1f%%; IHME %.1f%%; UN IGME %.1f%%).",
    f(h$model_deaths[h$year == 2024]), 100 * (h$ratio_ihme[h$year == 2024] - 1), 100 * (h$ratio_unigme[h$year == 2024] - 1),
    100 * (p0$ratio_ihme[p0$year == 2024] - 1), 100 * (p0$ratio_unigme[p0$year == 2024] - 1),
    -c1$rate_2000_2015, -c1$rate_2015_2024, -c0$rate_2000_2015, -c0$rate_2015_2024, -ihme_rate(2000, 2015), -ihme_rate(2015, 2024), -ig_rate(2000, 2015), -ig_rate(2015, 2024),
    -c1$rate_2000_2024, -c0$rate_2000_2024, -ihme_rate(2000, 2024), -ig_rate(2000, 2024)), "",
  sprintf("Countries in 2024: the model exceeds IHME in %d of 42 (primary %d); median country ratio %.2f (primary %.2f); Spearman correlation of rates with IHME %.2f (primary %.2f).",
    a1$countries_above_ihme, agree[floor_pct == 0]$countries_above_ihme, a1$median_ratio_ihme, agree[floor_pct == 0]$median_ratio_ihme, a1$spearman_rate_ihme, agree[floor_pct == 0]$spearman_rate_ihme), "",
  "## 2024 total by floor", "", "| Floor | Deaths | × IHME | × UN IGME | Rate per 1,000 child-years |", "|---|---:|---:|---:|---:|", fl_rows, "",
  "## 2024 deaths by age band, 1% floor (0%)", "", "| Age band | Deaths |", "|---|---:|",
  sprintf("| %s | %s (%s) |", ages, f(by1$attributable_deaths), f(by0$attributable_deaths)), "",
  "## Attributable fraction within each band, 1% floor (0%)", "", "| Age band | PfPR 10% | 20% | 30% | 40% |", "|---|---:|---:|---:|---:|",
  sprintf("| %s | %.1f%% (%.1f%%) | %.1f%% (%.1f%%) | %.1f%% (%.1f%%) | %.1f%% (%.1f%%) |", ages, 100 * af1$`10`, 100 * af0$`10`, 100 * af1$`20`, 100 * af0$`20`,
    100 * af1$`30`, 100 * af0$`30`, 100 * af1$`40`, 100 * af0$`40`), "",
  sprintf("## Malaria-caused probability of dying before 5 in 2024 (Figure 4 analogue)\n\nMedian across the 42 countries %.1f per 1,000 live births with the 1%% floor (primary %.1f); largest %s %.1f (primary %.1f).",
    1000 * median(q1$q5_malaria), 1000 * median(q0$q5_malaria), q1[which.max(q5_malaria)]$iso3, 1000 * max(q1$q5_malaria), 1000 * q0[iso3 == q1[which.max(q5_malaria)]$iso3]$q5_malaria), "",
  "## Not included", "",
  "Malaria deaths that would still occur at 1% transmission are excluded by construction. The World Malaria Report's low-transmission method would estimate them as 0.256% of *P. falciparum* cases, with an under-5 share from its quadratic in all-age mortality per 1,000 at risk. That needs a case source (MAP's prevalence-to-incidence curve at 1%, or reported cases), which is not in this repository. Under that method, 10-50 *P. falciparum* cases per 1,000 at risk correspond to roughly 0.04-0.26 under-5 deaths per 1,000 child-years (assuming under-5s are 16% of the population at risk).", "",
  "Files: `year_totals.csv`, `trend_changes.csv`, `country_year.csv`, `country_year_age.csv`, `agreement_2024.csv`, `country_comparison_2024_floor1.csv`, `deaths_by_age_2024.csv`, `attributable_fraction_by_age.csv`, `under5_probability_2024.csv`. Reproduce: `Rscript R_cbh/sensitivity/reference_floor/01_burden.R`.")
writeLines(lines, file.path(out, "REPORT.md"))
inputs <- c(unname(paths), "R_cbh/sensitivity/reference_floor/01_burden.R")
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")), file.path(out, "provenance.csv"))
message("Reference-floor burden written: ", out)
