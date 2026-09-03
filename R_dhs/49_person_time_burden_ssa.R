# =============================================================================
# 49_person_time_burden_ssa.R — malaria-attributable under-5 deaths for every
# sub-Saharan African country from the age-band person-time models, 2005-2025,
# under two all-cause denominators (UN IGME 2025 and IHME/GBD), against IHME's
# own malaria estimates by age block.
#
# Two tiers. Countries with a DHS/MIS in the panel (35) follow script 48: the
# latest survey's admin-1 regions carry MAP prevalence by year, the survey's
# regional death shares by age band allocate the national all-cause death counts
# (three age blocks: neonatal, 1-11 months, 1-4 years; DHS shares within blocks),
# and each band's fitted dose-response gives the attributable fraction. Countries
# without a panel survey (8) use the national population-weighted MAP prevalence
# and the pooled DHS age shares, so their attributable fraction ignores
# within-country heterogeneity.
#
# Denominators. The same attributable fractions are applied to two sets of
# national all-cause death counts: UN IGME 2025 (neonatal, 1-11 months, 1-4 years;
# ends 2024, carried into 2025) and IHME/GBD (early + late neonatal, 1-4 years,
# under 5; 1-11 months by subtraction; runs to 2025). The two sources disagree
# materially for some countries, and the comparison of the denominators is
# itself an output (allcause_igme_vs_ihme_by_country.csv, figure 34).
#
# IHME's malaria deaths come from the same GBD export by age (late neonatal,
# 1-4 years, under 5; 1-11 months by subtraction; IHME assigns no malaria deaths
# to the first week), so the model is compared with IHME block by block
# (figure 35) as well as in total (figure 33).
#
# The attributable fractions come from the Stan fits of script 43: 1,000
# posterior draws per band, the same draws reused for every country and both
# denominators, so country estimates are positively correlated and the
# sub-Saharan total's interval is the sum's, not a root-sum-of-squares.
# Intervals cover the dose-response only. MAP ends in 2024
# and is carried into 2025. The comparison year is 2024, the last with actual
# inputs on every side.
#
# Inputs (data/): the two IHME Data Explorer exports of 2026-09-03 copied to
#   ihme_allcause_u5_deaths_by_age_country_year.csv and
#   ihme_malaria_u5_deaths_by_age_country_year.csv (45 and 44 countries, 2000-2025).
#
# Outputs (results/dhs_rebuild)
#   person_time_burden_ssa_timeseries.csv       country x denominator x year, with IHME malaria by block
#   person_time_burden_ssa_2024_comparison.csv  country table: both denominators against IHME
#   person_time_burden_ssa_totals.csv           sub-Saharan totals by denominator and year
#   allcause_igme_vs_ihme_by_country.csv        all-cause under-5 deaths by block: IGME against IHME
#   figure33_person_time_burden_ssa.png         totals against IHME, both denominators
#   figure34_allcause_igme_vs_ihme.png          the denominators compared
#   figure35_burden_by_age_block_vs_ihme.png    model against IHME by age block
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("terra", "sf", "ggplot2", "countrycode", "patchwork"))

YEARS <- 2005:2025
LAST_DATA_YEAR <- 2024L
COMPARE_YEAR <- 2024L
N_DRAWS <- 1000L
CACHE <- file.path(DATA_DIR, "map_country_series_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
set.seed(20260828)
BLOCKS <- c("neonatal", "months_1_11", "years_1_4")
BLOCK_LABELS <- c(neonatal = "Neonatal", months_1_11 = "1-11 months", years_1_4 = "1-4 years")
IHME_ALLCAUSE_CSV <- file.path(DATA_DIR, "ihme_allcause_u5_deaths_by_age_country_year.csv")
IHME_MALARIA_CSV <- file.path(DATA_DIR, "ihme_malaria_u5_deaths_by_age_country_year.csv")

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
registry <- registry[registry$svkey %in% map$svkey, ]
latest <- do.call(rbind, lapply(split(registry, registry$iso3), function(z) z[which.max(z$year), ]))
person_time <- read.csv(file.path(DERIVED_DIR, "person_time_region_window_segment.csv"), stringsAsFactors = FALSE)
person_time$band <- band_from_segment(person_time$seg_lo, c(1, 4, 12, 24, 36), AGE6B)
national_pfpr <- read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE)
comparison_prev <- read.csv(file.path(RESULTS_DIR, "latest_country_burden_comparison.csv"), stringsAsFactors = FALSE)

# the country set: script 10's burden countries plus any panel DHS country it lacks
countries <- sort(union(comparison_prev$iso3, latest$iso3))
countries <- countries[countries %in% national_pfpr$iso3]
tier <- ifelse(countries %in% latest$iso3, "regional (latest DHS)", "national (no panel DHS)")
names(tier) <- countries
message(length(countries), " countries: ", sum(tier == "regional (latest DHS)"), " regional, ",
        sum(tier != "regional (latest DHS)"), " national")

## ---- all-cause death counts by age block: UN IGME 2025 and IHME/GBD ------------------------
fill_years <- function(w) {
  # restrict to the analysis years; a source that stops early carries its last year forward
  w <- w[w$year %in% YEARS, , drop = FALSE]
  missing <- setdiff(YEARS, w$year)
  if (length(missing)) {
    last <- w[w$year == max(w$year), , drop = FALSE]
    w <- rbind(w, do.call(rbind, lapply(missing, function(y) { l <- last; l$year <- y; l })))
  }
  w[order(w$year), , drop = FALSE]
}

IGME_BLOCKS <- c(neonatal = "Neonatal deaths", months_1_11 = "Deaths age 1-11 months",
                 years_1_4 = "Child deaths age 1 to 4")
IGME_EXTRACT <- file.path(DERIVED_DIR, "igme_2025_ssa.csv")
if (!file.exists(IGME_EXTRACT)) {
  csv <- file.path(DATA_DIR, "igme", "UN IGME 2025.csv")
  if (!file.exists(csv)) unzip(file.path(DATA_DIR, "igme", "UN-IGME-2025.zip"), exdir = file.path(DATA_DIR, "igme"))
  x <- read.csv(csv, stringsAsFactors = FALSE)
  est <- x[x$REF_AREA %in% countries & x$Sex == "Total" & x$Wealth.Quintile == "Total" &
             x$Series.Name == "UN IGME estimate" & x$Indicator %in% c(IGME_BLOCKS, "Under-five deaths"), ]
  est$year <- floor(est$REF_DATE)
  est <- est[, c("REF_AREA", "year", "Indicator", "OBS_VALUE", "LOWER_BOUND", "UPPER_BOUND")]
  names(est)[1] <- "iso3"
  write.csv(est, IGME_EXTRACT, row.names = FALSE)
  rm(x)
}
igme <- read.csv(IGME_EXTRACT, stringsAsFactors = FALSE)
igme$block <- setNames(c(names(IGME_BLOCKS), "under5"), c(IGME_BLOCKS, "Under-five deaths"))[igme$Indicator]
igme_blocks <- function(iso3) {
  z <- igme[igme$iso3 == iso3 & igme$block %in% BLOCKS, c("year", "block", "OBS_VALUE")]
  if (!nrow(z)) return(NULL)
  w <- reshape(z, idvar = "year", timevar = "block", direction = "wide")
  names(w) <- sub("^OBS_VALUE\\.", "", names(w))
  fill_years(w[, c("year", BLOCKS)])
}

# IHME exports carry early neonatal (0-6 days), late neonatal (7-27 days), 1-4 years
# and under 5; 1-11 months is the remainder. The malaria export has no early
# neonatal rows: IHME assigns no malaria deaths to the first week.
read_ihme_blocks <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[x$Measure == "Deaths" & x$Unit == "Number" & x$Sex == "Both", ]
  x$iso3 <- countrycode::countrycode(x$Location, "country.name", "iso3c", warn = FALSE)
  x <- x[!is.na(x$iso3), c("iso3", "Year", "Age", "Value", "Lower", "Upper")]
  w <- reshape(x[, c("iso3", "Year", "Age", "Value")], idvar = c("iso3", "Year"), timevar = "Age", direction = "wide")
  names(w) <- sub("^Value\\.", "", names(w))
  column <- function(name) if (name %in% names(w)) w[[name]] else 0
  out <- data.frame(iso3 = w$iso3, year = w$Year,
                    early_neonatal = column("0 to 6 days (early neonatal)"),
                    late_neonatal = column("7 to 27 days (late neonatal)"),
                    years_1_4 = column("1 to 4"), under5 = column("Under 5"), stringsAsFactors = FALSE)
  out$neonatal <- out$early_neonatal + out$late_neonatal
  out$months_1_11 <- pmax(out$under5 - out$neonatal - out$years_1_4, 0)
  bounds <- x[x$Age == "Under 5", c("iso3", "Year", "Lower", "Upper")]
  names(bounds) <- c("iso3", "year", "under5_lo", "under5_hi")
  out <- merge(out, bounds, by = c("iso3", "year"), all.x = TRUE)
  out[order(out$iso3, out$year), ]
}
ihme_allcause <- read_ihme_blocks(IHME_ALLCAUSE_CSV)
ihme_malaria <- read_ihme_blocks(IHME_MALARIA_CSV)
ihme_blocks <- function(iso3) {
  z <- ihme_allcause[ihme_allcause$iso3 == iso3, c("year", BLOCKS)]
  if (!nrow(z)) return(NULL)
  fill_years(z)
}
DENOMINATORS <- list(IGME = igme_blocks, IHME = ihme_blocks)
DENOMINATOR_LABELS <- c(IGME = "UN IGME all-cause deaths", IHME = "IHME/GBD all-cause deaths")

## ---- MAP prevalence: regional for DHS countries, national otherwise ---------------------
raster_files <- list.files(MAP_RASTER_DIR, pattern = "pfpr2_10_[0-9]{4}[.]tif$", full.names = TRUE)
raster_years <- as.integer(regmatches(basename(raster_files), regexpr("[0-9]{4}", basename(raster_files))))
rasters <- setNames(lapply(raster_files, terra::rast), as.character(raster_years))
weights_full <- terra::resample(terra::crop(terra::rast(GPW_TIF), rasters[[1]]), rasters[[1]], method = "bilinear")
regional_prevalence <- function(svkey) {
  cache <- file.path(CACHE, paste0(svkey, "_", min(YEARS), "_", LAST_DATA_YEAR, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  survey <- registry[registry$svkey == svkey, , drop = FALSE][1, ]
  boundary <- sf::st_make_valid(readRDS(file.path(BOUNDARY_DIR, paste0(survey$SurveyId, ".rds"))))
  polygons <- terra::makeValid(terra::vect(boundary))
  years <- YEARS[YEARS <= LAST_DATA_YEAR]
  extent <- terra::ext(polygons) + 0.5
  stack <- terra::crop(terra::rast(rasters[as.character(years)]), extent)
  w <- terra::mask(terra::crop(weights_full, extent), stack[[terra::nlyr(stack)]])
  denominator <- terra::extract(w, polygons, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)[[1]]
  numerator <- terra::extract(stack * w, polygons, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)
  out <- do.call(rbind, lapply(seq_along(years), function(j) data.frame(
    svkey = svkey, regkey = rkey(boundary$DHSREGEN), year = years[j],
    pfpr = 100 * numerator[[j]] / denominator, stringsAsFactors = FALSE)))
  out <- out[nzchar(out$regkey) & is.finite(out$pfpr), ]
  saveRDS(out, cache)
  out
}
national_matrix <- function(iso3) {
  p <- national_pfpr[national_pfpr$iso3 == iso3 & national_pfpr$year %in% YEARS, ]
  p <- p[order(p$year), ]
  matrix(p$pfpr_pct, nrow = 1, dimnames = list("national", p$year))
}
carry_forward <- function(m) {
  missing <- YEARS[!YEARS %in% colnames(m)]
  if (length(missing)) {
    extra <- m[, rep(as.character(LAST_DATA_YEAR), length(missing)), drop = FALSE]
    colnames(extra) <- missing
    m <- cbind(m, extra)
  }
  m[, as.character(YEARS), drop = FALSE]
}
prevalence_matrix <- function(iso3, regions, svkey = NULL) {
  nat <- carry_forward(national_matrix(iso3))
  if (is.null(svkey)) return(nat)
  p <- regional_prevalence(svkey)
  m <- tapply(p$pfpr, list(p$regkey, p$year), mean)
  m <- m[match(regions, rownames(m)), , drop = FALSE]
  rownames(m) <- regions
  m <- carry_forward(m)
  # a region with no usable MAP pixels (no population weight) follows the national series
  empty <- rowSums(is.na(m)) == ncol(m)
  if (any(empty)) {
    message("  ", iso3, ": ", sum(empty), " region(s) without MAP coverage take the national series: ",
            paste(regions[empty], collapse = ", "))
    m[empty, ] <- matrix(nat, nrow = sum(empty), ncol = ncol(m), byrow = TRUE)
  }
  # any remaining gap is filled from the nearest year
  for (i in seq_len(nrow(m))) if (anyNA(m[i, ])) m[i, ] <- stats::approx(seq_len(ncol(m)), m[i, ], xout = seq_len(ncol(m)), rule = 2)$y
  m
}

## ---- attributable-fraction draws from the Stan fits, once per band ------------------------------
# af_draws(band, prevalence): N_DRAWS x length(prevalence) draws of 1 - 1/HR against
# the 0% reference, from the posterior of script 43's brms fit for that band
af_draws <- person_time_stan_af_draws(N_DRAWS)

## ---- age bands within the three death blocks ------------------------------------------------------------
band_block <- setNames(c("neonatal", "months_1_11", "months_1_11", "years_1_4", "years_1_4", "years_1_4"), AGE6B)
infant_bands <- AGE6B[band_block == "months_1_11"]; child_bands <- AGE6B[band_block == "years_1_4"]
pooled_band <- tapply(person_time$deaths_w, person_time$band, sum)
pooled_within_infant <- pooled_band[infant_bands] / sum(pooled_band[infant_bands])
pooled_within_child <- pooled_band[child_bands] / sum(pooled_band[child_bands])

## ---- the burden, both denominators -------------------------------------------------------------------------
zero <- matrix(0, N_DRAWS, length(YEARS), dimnames = list(NULL, YEARS))
series_rows <- list()
ssa <- lapply(DENOMINATORS, function(...) list(total = zero, pn = zero, countries = 0L))
for (iso3 in countries) {
  blocks <- lapply(DENOMINATORS, function(f) f(iso3))
  if (is.null(blocks$IGME)) { message("  ", iso3, ": no UN IGME counts; skipped"); next }
  if (tier[[iso3]] == "regional (latest DHS)") {
    svkey <- latest$svkey[latest$iso3 == iso3]
    pt <- person_time[person_time$svkey == svkey, ]
    deaths_rb <- tapply(pt$deaths_w, list(pt$regkey, pt$band), sum); deaths_rb[is.na(deaths_rb)] <- 0
    deaths_rb <- deaths_rb[, AGE6B, drop = FALSE]
    regions <- rownames(deaths_rb)
    prev <- tryCatch(prevalence_matrix(iso3, regions, svkey), error = function(e) { message("  ", iso3, ": ", conditionMessage(e)); NULL })
    if (is.null(prev) || any(is.na(prev))) { message("  ", iso3, ": regional prevalence unavailable, using national"); tier[[iso3]] <- "national (no panel DHS)" }
  }
  if (tier[[iso3]] != "regional (latest DHS)") {
    regions <- "national"
    deaths_rb <- matrix(1, 1, length(AGE6B), dimnames = list("national", AGE6B))
    prev <- prevalence_matrix(iso3, regions)
    within_infant <- pooled_within_infant; within_child <- pooled_within_child
  } else {
    band_totals <- colSums(deaths_rb)
    within_infant <- band_totals[infant_bands] / sum(band_totals[infant_bands])
    within_child <- band_totals[child_bands] / sum(band_totals[child_bands])
    if (any(!is.finite(within_infant))) within_infant <- pooled_within_infant
    if (any(!is.finite(within_child))) within_child <- pooled_within_child
  }
  within <- setNames(c(1, within_infant[infant_bands], within_child[child_bands]), AGE6B)
  share_rb <- sweep(deaths_rb, 2, pmax(colSums(deaths_rb), 1e-9), "/")
  pfpr_weighted <- as.numeric(colSums(prev * rowSums(deaths_rb) / sum(deaths_rb)))
  # attributable fractions by band, region and year: the same for every denominator
  af <- lapply(AGE6B, function(band) { a <- af_draws(band, as.numeric(prev)); dim(a) <- c(N_DRAWS, length(regions), length(YEARS)); a })
  names(af) <- AGE6B
  for (den in names(DENOMINATORS)) {
    b <- blocks[[den]]
    if (is.null(b)) { message("  ", iso3, ": no ", den, " all-cause counts"); next }
    total <- zero; pn <- zero
    by_block <- matrix(0, length(BLOCKS), length(YEARS), dimnames = list(BLOCKS, YEARS))
    for (band in AGE6B) {
      national_band <- b[[band_block[[band]]]] * within[[band]]
      deaths <- sweep(af[[band]], c(2, 3), outer(share_rb[, band], national_band), "*")
      by_year <- apply(deaths, c(1, 3), sum)
      total <- total + by_year
      if (band != AGE6B[1]) pn <- pn + by_year
      by_block[band_block[[band]], ] <- by_block[band_block[[band]], ] + colMeans(by_year)
    }
    ssa[[den]]$total <- ssa[[den]]$total + total; ssa[[den]]$pn <- ssa[[den]]$pn + pn
    ssa[[den]]$countries <- ssa[[den]]$countries + 1L
    series_rows[[paste(iso3, den)]] <- data.frame(
      iso3 = iso3, tier = tier[[iso3]], denominator = den, year = YEARS, regions = length(regions),
      pfpr_death_weighted = pfpr_weighted,
      allcause_under5_deaths = b$neonatal + b$months_1_11 + b$years_1_4,
      allcause_neonatal = b$neonatal, allcause_months_1_11 = b$months_1_11, allcause_years_1_4 = b$years_1_4,
      malaria_deaths = colMeans(total), lo = apply(total, 2, stats::quantile, 0.025), hi = apply(total, 2, stats::quantile, 0.975),
      malaria_deaths_1_59m = colMeans(pn), lo_1_59m = apply(pn, 2, stats::quantile, 0.025), hi_1_59m = apply(pn, 2, stats::quantile, 0.975),
      malaria_neonatal = by_block["neonatal", ], malaria_months_1_11 = by_block["months_1_11", ],
      malaria_years_1_4 = by_block["years_1_4", ], stringsAsFactors = FALSE)
  }
  at <- as.character(COMPARE_YEAR)
  message(sprintf("  %s %-24s %2d regions  %d: IGME denominator %6.0fk, IHME denominator %s", iso3, tier[[iso3]],
                  length(regions), COMPARE_YEAR,
                  series_rows[[paste(iso3, "IGME")]]$malaria_deaths[YEARS == COMPARE_YEAR] / 1e3,
                  if (is.null(blocks$IHME)) "  n/a" else sprintf("%6.0fk", series_rows[[paste(iso3, "IHME")]]$malaria_deaths[YEARS == COMPARE_YEAR] / 1e3)))
}
series <- do.call(rbind, series_rows); rownames(series) <- NULL
series$country <- countrycode::countrycode(series$iso3, "iso3c", "country.name", warn = FALSE)

## ---- IHME malaria comparator, by block ---------------------------------------------------------------------
ihme_cmp <- ihme_malaria[, c("iso3", "year", "under5", "under5_lo", "under5_hi", "neonatal", "months_1_11", "years_1_4")]
names(ihme_cmp) <- c("iso3", "year", "ihme_u5", "ihme_lo", "ihme_hi", "ihme_neonatal", "ihme_months_1_11", "ihme_years_1_4")
series <- merge(series, ihme_cmp, by = c("iso3", "year"), all.x = TRUE)
series <- series[order(series$denominator, series$iso3, series$year), ]
write.csv(series, file.path(RESULTS_DIR, "person_time_burden_ssa_timeseries.csv"), row.names = FALSE)

## ---- country table at the comparison year ------------------------------------------------------------------
at <- series[series$year == COMPARE_YEAR, ]
pick <- function(den, col) { z <- at[at$denominator == den, ]; z[[col]][match(comparison$iso3, z$iso3)] }
base <- at[at$denominator == "IGME", ]
comparison <- data.frame(iso3 = base$iso3, country = base$country, tier = base$tier, regions = base$regions,
                         pfpr = round(base$pfpr_death_weighted, 1), stringsAsFactors = FALSE)
comparison$allcause_igme <- pick("IGME", "allcause_under5_deaths")
comparison$allcause_ihme <- pick("IHME", "allcause_under5_deaths")
comparison$ratio_allcause_igme_ihme <- comparison$allcause_igme / comparison$allcause_ihme
for (den in names(DENOMINATORS)) {
  key <- tolower(den)
  comparison[[paste0("model_", key)]] <- pick(den, "malaria_deaths")
  comparison[[paste0("model_", key, "_lo")]] <- pick(den, "lo")
  comparison[[paste0("model_", key, "_hi")]] <- pick(den, "hi")
  comparison[[paste0("model_", key, "_1_59m")]] <- pick(den, "malaria_deaths_1_59m")
  for (block in BLOCKS) comparison[[paste0("model_", key, "_", block)]] <- pick(den, paste0("malaria_", block))
}
for (col in c("ihme_u5", "ihme_lo", "ihme_hi", "ihme_neonatal", "ihme_months_1_11", "ihme_years_1_4"))
  comparison[[col]] <- pick("IGME", col)
comparison$ratio_model_igme_ihme <- comparison$model_igme / comparison$ihme_u5
comparison$ratio_model_ihme_ihme <- comparison$model_ihme / comparison$ihme_u5
prev_cols <- comparison_prev[, c("iso3", "best_model_deaths", "who_u5_proxy_2024", "who_all_age_2024")]
names(prev_cols) <- c("iso3", "region_pipeline_2024", "who_u5_proxy_2024", "who_all_age_2024")
comparison <- merge(comparison, prev_cols, by = "iso3", all.x = TRUE)
comparison <- comparison[order(-comparison$model_igme), ]
write.csv(comparison, file.path(RESULTS_DIR, "person_time_burden_ssa_2024_comparison.csv"), row.names = FALSE)

## ---- sub-Saharan totals ------------------------------------------------------------------------------------------
totals <- do.call(rbind, lapply(names(DENOMINATORS), function(den) {
  s <- series[series$denominator == den, ]
  data.frame(
    denominator = den, year = YEARS, countries = ssa[[den]]$countries,
    malaria_deaths = colMeans(ssa[[den]]$total), lo = apply(ssa[[den]]$total, 2, stats::quantile, 0.025),
    hi = apply(ssa[[den]]$total, 2, stats::quantile, 0.975),
    malaria_deaths_1_59m = colMeans(ssa[[den]]$pn), lo_1_59m = apply(ssa[[den]]$pn, 2, stats::quantile, 0.025),
    hi_1_59m = apply(ssa[[den]]$pn, 2, stats::quantile, 0.975),
    malaria_neonatal = tapply(s$malaria_neonatal, s$year, sum), malaria_months_1_11 = tapply(s$malaria_months_1_11, s$year, sum),
    malaria_years_1_4 = tapply(s$malaria_years_1_4, s$year, sum),
    allcause_under5 = tapply(s$allcause_under5_deaths, s$year, sum),
    # the IHME comparators and the model restricted to the countries IHME covers
    model_ihme_matched = tapply(ifelse(is.na(s$ihme_u5), NA, s$malaria_deaths), s$year, sum, na.rm = TRUE),
    ihme_matched_countries = tapply(s$ihme_u5, s$year, function(v) if (all(is.na(v))) NA else sum(v, na.rm = TRUE)),
    ihme_matched_n = tapply(s$ihme_u5, s$year, function(v) sum(!is.na(v))), stringsAsFactors = FALSE)
}))
rownames(totals) <- NULL
write.csv(totals, file.path(RESULTS_DIR, "person_time_burden_ssa_totals.csv"), row.names = FALSE)

## ---- the denominators compared: IGME against IHME all-cause deaths by block ----------------------------------------
igme_long <- igme[igme$iso3 %in% countries & igme$year %in% YEARS, c("iso3", "year", "block", "OBS_VALUE", "LOWER_BOUND", "UPPER_BOUND")]
names(igme_long) <- c("iso3", "year", "block", "igme", "igme_lo", "igme_hi")
ihme_long <- do.call(rbind, lapply(c(BLOCKS, "under5"), function(block) data.frame(
  iso3 = ihme_allcause$iso3, year = ihme_allcause$year, block = block, ihme = ihme_allcause[[block]],
  ihme_lo = if (block == "under5") ihme_allcause$under5_lo else NA, ihme_hi = if (block == "under5") ihme_allcause$under5_hi else NA,
  stringsAsFactors = FALSE)))
ihme_long <- ihme_long[ihme_long$iso3 %in% countries & ihme_long$year %in% YEARS, ]
allcause <- merge(igme_long, ihme_long, by = c("iso3", "year", "block"), all = TRUE)
allcause$ratio_igme_ihme <- allcause$igme / allcause$ihme
allcause$country <- countrycode::countrycode(allcause$iso3, "iso3c", "country.name", warn = FALSE)
allcause <- allcause[order(allcause$block, allcause$iso3, allcause$year), c("iso3", "country", "year", "block", "igme", "igme_lo", "igme_hi", "ihme", "ihme_lo", "ihme_hi", "ratio_igme_ihme")]
write.csv(allcause, file.path(RESULTS_DIR, "allcause_igme_vs_ihme_by_country.csv"), row.names = FALSE)
both <- allcause[complete.cases(allcause[, c("igme", "ihme")]), ]
matched_iso <- unique(both$iso3)
allcause_totals <- aggregate(cbind(igme, ihme) ~ year + block, data = both, FUN = sum)
allcause_totals$ratio_igme_ihme <- allcause_totals$igme / allcause_totals$ihme

## ---- console summary --------------------------------------------------------------------------------------------------
at_year <- function(den, col) totals[[col]][totals$denominator == den & totals$year == COMPARE_YEAR]
with_ihme <- comparison[is.finite(comparison$ihme_u5), ]
message(sprintf("\nAll-cause under-5 deaths, %d, %d matched countries: IGME %.0fk, IHME %.0fk (IGME/IHME %.2f)",
                COMPARE_YEAR, length(matched_iso),
                allcause_totals$igme[allcause_totals$year == COMPARE_YEAR & allcause_totals$block == "under5"] / 1e3,
                allcause_totals$ihme[allcause_totals$year == COMPARE_YEAR & allcause_totals$block == "under5"] / 1e3,
                allcause_totals$ratio_igme_ihme[allcause_totals$year == COMPARE_YEAR & allcause_totals$block == "under5"]))
for (block in BLOCKS) message(sprintf("  %-12s IGME %.0fk, IHME %.0fk, ratio %.2f", BLOCK_LABELS[[block]],
                                     allcause_totals$igme[allcause_totals$year == COMPARE_YEAR & allcause_totals$block == block] / 1e3,
                                     allcause_totals$ihme[allcause_totals$year == COMPARE_YEAR & allcause_totals$block == block] / 1e3,
                                     allcause_totals$ratio_igme_ihme[allcause_totals$year == COMPARE_YEAR & allcause_totals$block == block]))
u5 <- both[both$block == "under5" & both$year == COMPARE_YEAR, ]
message(sprintf("  Country ratios IGME/IHME under 5: median %.2f, IQR %.2f-%.2f; largest %s (%.2f), smallest %s (%.2f)",
                stats::median(u5$ratio_igme_ihme), stats::quantile(u5$ratio_igme_ihme, 0.25), stats::quantile(u5$ratio_igme_ihme, 0.75),
                u5$iso3[which.max(u5$ratio_igme_ihme)], max(u5$ratio_igme_ihme), u5$iso3[which.min(u5$ratio_igme_ihme)], min(u5$ratio_igme_ihme)))
for (den in names(DENOMINATORS)) {
  message(sprintf("\n%s denominator, %d countries. Sub-Saharan total %d: %.0fk (%.0fk-%.0fk), 1-59 months %.0fk; blocks neonatal %.0fk, 1-11 m %.0fk, 1-4 y %.0fk",
                  den, ssa[[den]]$countries, COMPARE_YEAR, at_year(den, "malaria_deaths") / 1e3, at_year(den, "lo") / 1e3, at_year(den, "hi") / 1e3,
                  at_year(den, "malaria_deaths_1_59m") / 1e3, at_year(den, "malaria_neonatal") / 1e3,
                  at_year(den, "malaria_months_1_11") / 1e3, at_year(den, "malaria_years_1_4") / 1e3))
  message(sprintf("  Same %d countries as IHME: model %.0fk against IHME %.0fk (ratio %.2f)", at_year(den, "ihme_matched_n"),
                  at_year(den, "model_ihme_matched") / 1e3, at_year(den, "ihme_matched_countries") / 1e3,
                  at_year(den, "model_ihme_matched") / at_year(den, "ihme_matched_countries")))
  change <- ssa[[den]]$total[, as.character(max(YEARS))] / ssa[[den]]$total[, "2005"] - 1
  message(sprintf("  2005 -> %d: %.0fk -> %.0fk, %+.0f%% (%+.0f%% to %+.0f%%)", max(YEARS),
                  totals$malaria_deaths[totals$denominator == den & totals$year == 2005] / 1e3,
                  totals$malaria_deaths[totals$denominator == den & totals$year == max(YEARS)] / 1e3,
                  100 * mean(change), 100 * stats::quantile(change, 0.025), 100 * stats::quantile(change, 0.975)))
  r <- with_ihme[[paste0("ratio_model_", tolower(den), "_ihme")]]
  message(sprintf("  Country ratio model/IHME: median %.2f, IQR %.2f-%.2f; %d of %d above IHME",
                  stats::median(r, na.rm = TRUE), stats::quantile(r, 0.25, na.rm = TRUE), stats::quantile(r, 0.75, na.rm = TRUE),
                  sum(r > 1, na.rm = TRUE), sum(is.finite(r))))
}
message(sprintf("IHME malaria, same countries, 2005 -> 2024: %+.0f%%",
                100 * (at_year("IGME", "ihme_matched_countries") / totals$ihme_matched_countries[totals$denominator == "IGME" & totals$year == 2005] - 1)))
ihme_blocks_2024 <- colSums(with_ihme[, c("ihme_neonatal", "ihme_months_1_11", "ihme_years_1_4")])
message(sprintf("IHME malaria %d by block, %d countries: neonatal %.0fk (%.0f%%), 1-11 m %.0fk (%.0f%%), 1-4 y %.0fk (%.0f%%)",
                COMPARE_YEAR, nrow(with_ihme), ihme_blocks_2024[1] / 1e3, 100 * ihme_blocks_2024[1] / sum(ihme_blocks_2024),
                ihme_blocks_2024[2] / 1e3, 100 * ihme_blocks_2024[2] / sum(ihme_blocks_2024),
                ihme_blocks_2024[3] / 1e3, 100 * ihme_blocks_2024[3] / sum(ihme_blocks_2024)))
message("\nCountry table, ", COMPARE_YEAR, " (thousands):")
print(transform(comparison,
                allcause = sprintf("%.0f / %s", allcause_igme / 1e3, ifelse(is.na(allcause_ihme), "-", sprintf("%.0f", allcause_ihme / 1e3))),
                ac_ratio = round(ratio_allcause_igme_ihme, 2),
                model_igme = sprintf("%.0f (%.0f-%.0f)", model_igme / 1e3, model_igme_lo / 1e3, model_igme_hi / 1e3),
                model_ihme = ifelse(is.na(model_ihme), "-", sprintf("%.0f", model_ihme / 1e3)),
                ihme = ifelse(is.na(ihme_u5), "-", sprintf("%.0f (%.0f-%.0f)", ihme_u5 / 1e3, ihme_lo / 1e3, ihme_hi / 1e3)),
                r_igme = round(ratio_model_igme_ihme, 2), r_ihme = round(ratio_model_ihme_ihme, 2),
                tier = substr(tier, 1, 8))[
                  , c("iso3", "tier", "regions", "pfpr", "allcause", "ac_ratio", "model_igme", "model_ihme", "ihme", "r_igme", "r_ihme")],
      row.names = FALSE, right = FALSE)

## ---- figure 33: totals against IHME, both denominators ----------------------------------------------------------------
theme_small <- ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(legend.position = "bottom")
# countries where IHME puts fewer than 100 under-5 malaria deaths sit at the
# bottom of a log axis and stretch it without adding information
scatter_df <- with_ihme[with_ihme$ihme_u5 >= 100 & with_ihme$model_igme >= 100, ]
scatter_df$tier_short <- ifelse(grepl("^regional", scatter_df$tier), "Admin-1 (latest DHS)", "National prevalence only")
tiny <- setdiff(with_ihme$iso3, scatter_df$iso3)
points_df <- rbind(
  data.frame(iso3 = scatter_df$iso3, x = scatter_df$ihme_u5, y = scatter_df$model_igme, denominator = DENOMINATOR_LABELS[["IGME"]], tier_short = scatter_df$tier_short),
  data.frame(iso3 = scatter_df$iso3, x = scatter_df$ihme_u5, y = scatter_df$model_ihme, denominator = DENOMINATOR_LABELS[["IHME"]], tier_short = scatter_df$tier_short))
points_df <- points_df[is.finite(points_df$y), ]
points_df$denominator <- factor(points_df$denominator, levels = DENOMINATOR_LABELS)
scatter <- ggplot2::ggplot(scatter_df, ggplot2::aes(ihme_u5 / 1e3, model_igme / 1e3)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = model_igme_lo / 1e3, ymax = model_igme_hi / 1e3), colour = "grey75", width = 0) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ihme_lo / 1e3, xmax = ihme_hi / 1e3), colour = "grey75", height = 0) +
  ggplot2::geom_segment(ggplot2::aes(xend = ihme_u5 / 1e3, y = model_igme / 1e3, yend = model_ihme / 1e3), colour = "#1D6F8B", linewidth = 0.4) +
  ggplot2::geom_point(data = points_df, ggplot2::aes(x / 1e3, y / 1e3, colour = denominator, shape = tier_short), size = 2.2) +
  ggplot2::geom_text(ggplot2::aes(label = iso3), size = 2.4, vjust = -0.9, colour = "grey25") +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = c("#1D6F8B", "#E08214"), name = NULL) +
  ggplot2::scale_shape_manual(values = c(16, 1), name = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), shape = ggplot2::guide_legend(nrow = 2)) +
  ggplot2::labs(x = sprintf("IHME/GBD under-5 malaria deaths, %d (thousands, log)", COMPARE_YEAR),
                y = "Model, all under-5 (thousands, log)",
                title = sprintf("Country estimates against IHME, %d", COMPARE_YEAR),
                subtitle = sprintf("Each country twice: attributable fractions applied to IGME and to IHME all-cause deaths.\nBars are 95%% intervals (model: dose-response only, IGME denominator). Not shown, under 100 deaths on either side: %s",
                                   paste(tiny, collapse = ", "))) +
  theme_small
trend_df <- rbind(
  data.frame(year = YEARS, deaths = totals$malaria_deaths[totals$denominator == "IGME"], lo = totals$lo[totals$denominator == "IGME"],
             hi = totals$hi[totals$denominator == "IGME"], series = "Model, IGME all-cause denominator"),
  data.frame(year = YEARS, deaths = totals$malaria_deaths[totals$denominator == "IHME"], lo = NA, hi = NA, series = "Model, IHME all-cause denominator"),
  data.frame(year = YEARS, deaths = totals$malaria_deaths_1_59m[totals$denominator == "IGME"], lo = NA, hi = NA, series = "Model, 1-59 months (IGME denominator)"),
  data.frame(year = YEARS, deaths = totals$ihme_matched_countries[totals$denominator == "IGME"], lo = NA, hi = NA, series = "IHME/GBD malaria, its countries"))
trend_df$series <- factor(trend_df$series, levels = unique(trend_df$series))
trend <- ggplot2::ggplot(trend_df, ggplot2::aes(year, deaths / 1e3, colour = series)) +
  ggplot2::geom_ribbon(data = trend_df[!is.na(trend_df$lo), ], ggplot2::aes(ymin = lo / 1e3, ymax = hi / 1e3), fill = "#1D6F8B", alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_vline(xintercept = LAST_DATA_YEAR + 0.5, linetype = "dotted", colour = "grey50") +
  ggplot2::scale_colour_manual(values = c("#1D6F8B", "#E08214", "#7FB3C8", "#B2182B"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
  ggplot2::labs(x = NULL, y = "Malaria-attributable under-5 deaths,\nsub-Saharan Africa (thousands)",
                title = sprintf("Sub-Saharan total: %d countries (IGME denominator), %d (IHME)",
                                ssa$IGME$countries, ssa$IHME$countries),
                subtitle = "Right of the dotted line MAP and IGME are 2024 carried forward; IHME runs to 2025") +
  theme_small
ggplot2::ggsave(file.path(RESULTS_DIR, "figure33_person_time_burden_ssa.png"), scatter | trend,
                width = 13, height = 6.4, dpi = 200, bg = "white")

## ---- figure 34: the denominators compared ------------------------------------------------------------------------------
ac_u5 <- both[both$block == "under5" & both$year == COMPARE_YEAR, ]
panel_a <- ggplot2::ggplot(ac_u5, ggplot2::aes(ihme / 1e3, igme / 1e3)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = igme_lo / 1e3, ymax = igme_hi / 1e3), colour = "grey75", width = 0) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ihme_lo / 1e3, xmax = ihme_hi / 1e3), colour = "grey75", height = 0) +
  ggplot2::geom_point(colour = "#1D6F8B", size = 2.2) +
  ggplot2::geom_text(ggplot2::aes(label = iso3), size = 2.4, vjust = -0.9, colour = "grey25") +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::labs(x = sprintf("IHME/GBD all-cause under-5 deaths, %d (thousands, log)", COMPARE_YEAR),
                y = "UN IGME 2025 all-cause under-5 deaths (thousands, log)",
                title = "All-cause under-5 deaths: the two denominators",
                subtitle = "Dashed line is equality; bars are each source's 95% interval") +
  theme_small
ratio_df <- both[both$year == COMPARE_YEAR & both$block %in% BLOCKS, ]
ratio_df$block <- factor(BLOCK_LABELS[ratio_df$block], levels = BLOCK_LABELS)
order_iso <- ac_u5$iso3[order(ac_u5$ratio_igme_ihme)]
ratio_df$iso3 <- factor(ratio_df$iso3, levels = order_iso)
u5_pts <- data.frame(iso3 = factor(ac_u5$iso3, levels = order_iso), ratio_igme_ihme = ac_u5$ratio_igme_ihme)
panel_b <- ggplot2::ggplot(ratio_df, ggplot2::aes(ratio_igme_ihme, iso3)) +
  ggplot2::geom_vline(xintercept = 1, colour = "grey50") +
  ggplot2::geom_point(data = u5_pts, shape = 124, size = 4, colour = "grey30") +
  ggplot2::geom_point(ggplot2::aes(colour = block), size = 1.9, alpha = 0.9) +
  ggplot2::scale_x_log10(breaks = c(0.5, 0.7, 1, 1.5, 2, 3)) +
  ggplot2::scale_colour_manual(values = c("#7570B3", "#1B9E77", "#D95F02"), name = NULL) +
  ggplot2::labs(x = "IGME / IHME all-cause deaths (log scale)", y = NULL,
                title = sprintf("Ratio by age block, %d", COMPARE_YEAR),
                subtitle = "Vertical tick is the under-5 total; countries ordered by it") +
  theme_small + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
tot_df <- allcause_totals[allcause_totals$block %in% BLOCKS & allcause_totals$year <= LAST_DATA_YEAR, ]
tot_long <- rbind(data.frame(year = tot_df$year, block = tot_df$block, deaths = tot_df$igme, source = "UN IGME 2025"),
                  data.frame(year = tot_df$year, block = tot_df$block, deaths = tot_df$ihme, source = "IHME/GBD"))
tot_long$block <- factor(BLOCK_LABELS[tot_long$block], levels = BLOCK_LABELS)
panel_c <- ggplot2::ggplot(tot_long, ggplot2::aes(year, deaths / 1e6, colour = block, linetype = source)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::scale_colour_manual(values = c("#7570B3", "#1B9E77", "#D95F02"), name = NULL) +
  ggplot2::scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(x = NULL, y = "All-cause deaths (millions)",
                title = sprintf("Sub-Saharan totals by age block, %d matched countries", length(matched_iso))) +
  theme_small
ggplot2::ggsave(file.path(RESULTS_DIR, "figure34_allcause_igme_vs_ihme.png"),
                (panel_a | panel_b) / panel_c + patchwork::plot_layout(heights = c(1.25, 0.75)),
                width = 12.5, height = 11, dpi = 200, bg = "white")

## ---- figure 35: model against IHME by age block ------------------------------------------------------------------------
block_bars <- rbind(
  do.call(rbind, lapply(names(DENOMINATORS), function(den) data.frame(
    estimate = paste0("Model, ", DENOMINATOR_LABELS[[den]]), block = BLOCKS,
    deaths = colSums(with_ihme[, paste0("model_", tolower(den), "_", BLOCKS)], na.rm = TRUE)))),
  data.frame(estimate = "IHME/GBD malaria", block = BLOCKS, deaths = as.numeric(ihme_blocks_2024)))
block_bars$block <- factor(BLOCK_LABELS[block_bars$block], levels = rev(BLOCK_LABELS))
block_bars$estimate <- factor(block_bars$estimate, levels = unique(block_bars$estimate))
block_bars$share <- ave(block_bars$deaths, block_bars$estimate, FUN = function(v) v / sum(v))
panel_bars <- ggplot2::ggplot(block_bars, ggplot2::aes(estimate, deaths / 1e3, fill = block)) +
  ggplot2::geom_col(width = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0fk (%.0f%%)", deaths / 1e3, 100 * share)),
                     position = ggplot2::position_stack(vjust = 0.5), size = 2.8, colour = "white") +
  ggplot2::scale_fill_manual(values = rev(c("#7570B3", "#1B9E77", "#D95F02")), name = NULL) +
  ggplot2::scale_x_discrete(labels = function(l) gsub(", ", ",\n", l)) +
  ggplot2::labs(x = NULL, y = "Under-5 malaria deaths (thousands)",
                title = sprintf("Age composition, %d, %d countries with IHME estimates", COMPARE_YEAR, nrow(with_ihme))) +
  theme_small
block_scatter <- do.call(rbind, lapply(BLOCKS, function(block) data.frame(
  iso3 = with_ihme$iso3, block = block, ihme = with_ihme[[paste0("ihme_", block)]],
  model = with_ihme[[paste0("model_igme_", block)]], stringsAsFactors = FALSE)))
block_scatter <- block_scatter[is.finite(block_scatter$ihme) & is.finite(block_scatter$model) &
                                 block_scatter$ihme >= 20 & block_scatter$model >= 20, ]
block_scatter$block <- factor(BLOCK_LABELS[block_scatter$block], levels = BLOCK_LABELS)
panel_blocks <- ggplot2::ggplot(block_scatter, ggplot2::aes(ihme / 1e3, model / 1e3)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = block), size = 2, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = iso3), size = 2.2, vjust = -0.9, colour = "grey25") +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = c("#7570B3", "#1B9E77", "#D95F02")) +
  ggplot2::facet_wrap(~ block, nrow = 1) +
  ggplot2::labs(x = sprintf("IHME/GBD malaria deaths in the block, %d (thousands, log)", COMPARE_YEAR),
                y = "Model, IGME denominator (thousands, log)",
                title = "Country estimates by age block",
                subtitle = "Dashed line is equality; countries with fewer than 20 deaths on either side omitted") +
  theme_small
ggplot2::ggsave(file.path(RESULTS_DIR, "figure35_burden_by_age_block_vs_ihme.png"),
                panel_bars / panel_blocks + patchwork::plot_layout(heights = c(1, 1)),
                width = 12.5, height = 10, dpi = 200, bg = "white")
message("\nWrote person_time_burden_ssa_timeseries.csv, person_time_burden_ssa_2024_comparison.csv, ",
        "person_time_burden_ssa_totals.csv, allcause_igme_vs_ihme_by_country.csv and figures 33-35")
