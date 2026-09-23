#!/usr/bin/env Rscript
# Figure 5 inputs and annual point estimates. No model refitting or downloads.
# Run from project root; --audit-only checks inputs before raster extraction.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mgcv))
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional_mics"))
out <- file.path(settings$out, "annual_comparison")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
years <- 2000:2024
ages <- cbh_config()$age_bands$age_band
allcause_path <- list.files("data", pattern = "2026-09-09 10-58-22[.]csv$", full.names = TRUE)
stopifnot(length(allcause_path) == 1)
malaria_path <- "data/ihme_malaria_u5_deaths_by_age_country_year.csv"
who_path <- "data/external/figure5/unicef_cacode_u5_malaria_2000_2024.csv"
component_path <- file.path(settings$private, "pfpr_components.rds")
baseline_path <- file.path(settings$out, "burden/country_totals.csv")
pop_path <- "data/pop/gpw_v4_population_density_rev11_2020_2.5m.tif"
boundary_path <- "data/africa_admin0.rds"
map_paths <- setNames(ifelse(years == 2024, "data/pfpr_2to10_africa_2024.tif",
  sprintf("data/map_annual/pfpr2_10_%d.tif", years)), years)
inputs <- c(allcause_path, malaria_path, who_path, "data/external/figure5/source.json",
  component_path, baseline_path, pop_path, boundary_path, map_paths,
  file.path(settings$out, c("fit_manifest.csv", "fit_diagnostics.csv")))
stopifnot(all(file.exists(inputs)))
previous <- fread(baseline_path)
countries <- unique(previous[status == "estimated", .(iso3, country)])
setorder(countries, iso3)
stopifnot(nrow(countries) == 42, all(previous[status == "estimated", .N, by = year]$N == 42))
components <- readRDS(component_path)
manifest <- fread(file.path(settings$out, "fit_manifest.csv"))
diagnostics <- fread(file.path(settings$out, "fit_diagnostics.csv"))
stopifnot(length(components) == 7, all(diagnostics$converged), all(diagnostics$input_verified))
for (piece in components) stopifnot(identical(piece$model_md5,
  manifest$md5[match(piece$fit_id, manifest$fit_id)]))

read_ihme <- function(path, cause) {
  z <- fread(path)[Year %in% years & Sex == "Both" & Condition == cause & Measure == "Deaths"]
  z[, iso3 := countrycode::countrycode(Location, "country.name", "iso3c", warn = FALSE)]
  stopifnot(!anyNA(z$iso3), all(is.finite(z$Value)), all(z$Value >= 0))
  cbh_unique(as.data.frame(z), c("iso3", "Year", "Age", "Unit"), path)
  z
}
allcause <- read_ihme(allcause_path, "All causes")
malaria <- read_ihme(malaria_path, "Malaria")[Age == "Under 5"]
who <- fread(who_path)
stopifnot(all(who$INDICATOR == "DEATHS"), all(who$CAUSE_OF_DEATH == "MALARIA"),
  all(who$SEX == "_T"), all(who$AGE_GROUP == "Y0T4"), all(who$UNIT_MEASURE == "D"),
  all(who$SERIES_NAME == "CA CODE"))
who <- who[REF_AREA %in% countries$iso3 & TIME_PERIOD %in% years,
  .(iso3 = REF_AREA, year = TIME_PERIOD, who_cacode_deaths = OBS_VALUE,
    who_cacode_lower90 = LOWER_BOUND, who_cacode_upper90 = UPPER_BOUND)]
cbh_unique(as.data.frame(who), c("iso3", "year"), "WHO/UNICEF CA-CODE")
stopifnot(nrow(who) == 42 * length(years), all(is.finite(who$who_cacode_deaths)),
  all(who$who_cacode_deaths >= 0))

wide <- dcast(allcause, iso3 + Location + Year + Age ~ Unit, value.var = "Value")
wide[, implied_person_years := Number / `Rate (per 100,000)` * 1e5]
stopifnot(all(is.finite(wide$implied_person_years)), all(wide$implied_person_years > 0))
leaf_ages <- c("0 to 6 days (early neonatal)", "7 to 27 days (late neonatal)",
  "1 to 5 months", "6-11 months", "12 to 23 months", "2 to 4")
leaf_sums <- wide[Age %in% leaf_ages, .(leaf_deaths = sum(Number),
  leaf_person_years = sum(implied_person_years), age_rows = .N), by = .(iso3, Year)]
denominators <- merge(wide[Age == "Under 5", .(iso3, Year,
  allcause_deaths = Number, under5_person_years = implied_person_years)], leaf_sums,
  by = c("iso3", "Year"))
denominators[, `:=`(death_sum_relative_error = leaf_deaths / allcause_deaths - 1,
  person_year_sum_relative_error = leaf_person_years / under5_person_years - 1)]
stopifnot(all(denominators$age_rows == 6),
  max(abs(denominators$death_sum_relative_error)) < 1e-5)
m <- dcast(malaria, iso3 + Year ~ Unit, value.var = "Value")
m[, malaria_native_person_years := fifelse(`Rate (per 100,000)` > 0,
  Number / `Rate (per 100,000)` * 1e5, NA_real_)]
setnames(m, c("Year", "Number", "Rate (per 100,000)"),
  c("year", "ihme_malaria_deaths", "ihme_native_rate_per100000"))
setnames(denominators, "Year", "year")
denominators <- merge(denominators, m, by = c("iso3", "year"), all.x = TRUE)
denominators[, malaria_denominator_relative_difference := malaria_native_person_years / under5_person_years - 1]
cbh_atomic_csv(denominators, file.path(out, "denominator_audit.csv"))

# Check every required country-year, never allow a changing annual country set.
grid <- CJ(iso3 = countries$iso3, year = years)
coverage <- merge(grid, denominators, by = c("iso3", "year"), all.x = TRUE)
coverage <- merge(coverage, who, by = c("iso3", "year"), all.x = TRUE)
coverage[, complete := is.finite(allcause_deaths) & is.finite(under5_person_years) &
  is.finite(ihme_malaria_deaths) & is.finite(who_cacode_deaths)]
cbh_atomic_csv(coverage, file.path(out, "country_year_input_audit.csv"))
stopifnot(nrow(coverage) == 42 * length(years), all(coverage$complete))
map_inventory <- rbindlist(lapply(years, function(y) {
  p <- terra::rast(map_paths[as.character(y)])
  stopifnot(terra::nlyr(p) == 1, terra::is.lonlat(p))
  data.table(year = y, path = map_paths[as.character(y)],
    rows = terra::nrow(p), columns = terra::ncol(p), readable = TRUE)
}))
cbh_atomic_csv(map_inventory, file.path(out, "map_input_audit.csv"))
cbh_atomic_csv(countries, file.path(out, "included_countries.csv"))
metadata <- rbindlist(lapply(list(allcause, malaria), function(z)
  unique(z[, .(Condition, Data_Suite, Data_Type, Data_Type_Level)])))
cbh_atomic_csv(metadata, file.path(out, "ihme_release_metadata.csv"))
message("Input audit passed: 42 countries x ", length(years), " years for all three death series; all MAP rasters readable.")
if ("--audit-only" %in% commandArgs(TRUE)) quit(save = "no", status = 0)

# Fixed GPW 2020 spatial weights, as in the primary burden calculation.
# Density must be multiplied by cell area before national population weighting.
b <- sf::st_make_valid(readRDS(boundary_path))
b <- b[b$iso %in% countries$iso3, ]
stopifnot(setequal(b$iso, countries$iso3), !anyDuplicated(b$iso))
nat <- rbindlist(lapply(years, function(y) {
  message("National MAP extraction: ", y)
  p <- terra::rast(map_paths[as.character(y)])
  density <- terra::resample(terra::crop(terra::rast(pop_path), p), p, method = "bilinear")
  weights <- terra::ifel(density >= 0, density, NA) * terra::cellSize(p, unit = "km")
  covered <- terra::mask(weights, p)
  layers <- c(covered, p * covered, weights)
  names(layers) <- c("covered", "numerator", "all_pop")
  agg <- exactextractr::exact_extract(layers, b, "sum", progress = FALSE)
  data.table(iso3 = b$iso, year = y, pfpr_pct = 100 * agg$sum.numerator / agg$sum.covered,
    map_population_coverage = agg$sum.covered / agg$sum.all_pop,
    population_weight_year = 2020L)
}))
stopifnot(nrow(nat) == nrow(grid), all(is.finite(nat$pfpr_pct)),
  all(nat$pfpr_pct >= 0 & nat$pfpr_pct <= 100))
cbh_atomic_csv(nat, file.path(out, "national_pfpr_2000_2024.csv"))

bands <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  iso <- grid$iso3[i]; y <- grid$year[i]
  z <- wide[iso3 == iso & Year == y & Age %in% leaf_ages]
  z <- z[match(leaf_ages, Age)]
  stopifnot(nrow(z) == 6, !anyNA(z$Age))
  data.table(iso3 = iso, year = y, age_band = ages,
    allcause_deaths = c(sum(z$Number[1:2]), z$Number[3:5], rep(z$Number[6] / 3, 3)),
    age_person_years = c(sum(z$implied_person_years[1:2]), z$implied_person_years[3:5],
      rep(z$implied_person_years[6] / 3, 3)))
}))
bands <- merge(bands, nat, by = c("iso3", "year"))
bands <- rbindlist(lapply(seq_along(ages), function(i) {
  z <- bands[age_band == ages[i]]
  piece <- components[[paste0("map_full_age_", i)]]
  L <- PredictMat(piece$smooth, data.frame(pfpr_pct = rep(0, nrow(z)))) -
    PredictMat(piece$smooth, data.frame(pfpr_pct = z$pfpr_pct))
  z[, log_hr_zero_vs_current := drop(L %*% piece$coef)]
  z[, attributable_fraction := -expm1(log_hr_zero_vs_current)]
  z[, attributable_deaths := allcause_deaths * attributable_fraction]
  z[, `:=`(zero_below_observed_support = piece$support[1] > 0,
    current_outside_observed_support = pfpr_pct < piece$support[1] | pfpr_pct > piece$support[4])]
  z
}))
stopifnot(all(is.finite(bands$attributable_deaths)))
cbh_atomic_csv(bands, file.path(out, "country_age_estimates_2000_2024.csv"))
totals <- bands[, .(model_deaths = sum(attributable_deaths),
  any_current_outside_support = any(current_outside_observed_support)), by = .(iso3, year)]
totals <- merge(totals, coverage, by = c("iso3", "year"))
totals <- merge(totals, nat, by = c("iso3", "year"))
totals <- merge(totals, countries, by = "iso3")
for (s in c("model", "ihme_malaria", "who_cacode"))
  totals[, (paste0(s, "_rate_per100000")) := get(paste0(s, "_deaths")) / under5_person_years * 1e5]
setorder(totals, year, iso3)
cbh_atomic_csv(totals, file.path(out, "country_estimates_2000_2024.csv"))

# Numerical regression check against all 126 existing primary country estimates.
check <- merge(totals, previous[status == "estimated", .(iso3, year,
  previous_deaths = attributable_under5_deaths, previous_pfpr = pfpr_pct,
  previous_ihme = ihme_malaria_deaths)], by = c("iso3", "year"))
check[, `:=`(death_difference = model_deaths - previous_deaths,
  pfpr_difference = pfpr_pct - previous_pfpr,
  ihme_difference = ihme_malaria_deaths - previous_ihme)]
stopifnot(nrow(check) == 126, max(abs(check$death_difference)) < .01,
  max(abs(check$pfpr_difference)) < 1e-6, max(abs(check$ihme_difference)) < 1e-6)
cbh_atomic_csv(check[, .(iso3, year, death_difference, pfpr_difference, ihme_difference)],
  file.path(out, "agreement_with_existing_primary.csv"))
annual <- totals[, .(countries = .N, under5_person_years = sum(under5_person_years),
  model_deaths = sum(model_deaths), ihme_malaria_deaths = sum(ihme_malaria_deaths),
  who_cacode_deaths = sum(who_cacode_deaths)), by = year]
for (s in c("model", "ihme_malaria", "who_cacode"))
  annual[, (paste0(s, "_rate_per100000")) := get(paste0(s, "_deaths")) / under5_person_years * 1e5]
stopifnot(nrow(annual) == length(years), all(annual$countries == 42))
cbh_atomic_csv(annual, file.path(out, "annual_totals_2000_2024.csv"))
inputs <- unique(c(inputs, "R_cbh/burden/04_annual_comparison.R", "R_cbh/primary/settings.R"))
cbh_atomic_csv(data.frame(file = inputs, md5 = vapply(inputs, cbh_file_hash, "")),
  file.path(out, "input_provenance.csv"))
print(annual[year %in% c(2000, 2004, 2005, 2015, 2024)])
