# =============================================================================
# 49_person_time_burden_ssa.R — malaria-attributable under-5 deaths for every
# sub-Saharan African country from the age-band person-time models, 2005-2025,
# against IHME/GBD.
#
# Two tiers. Countries with a DHS/MIS in the panel (36) follow script 48: the
# latest survey's admin-1 regions carry MAP prevalence by year, the survey's
# regional death shares by age band allocate UN IGME's national death counts
# (three age blocks, DHS shares within blocks), and each band's fitted
# dose-response gives the attributable fraction. Countries without a panel survey
# (8) use the national population-weighted MAP prevalence and the pooled DHS age
# shares, so their attributable fraction ignores within-country heterogeneity.
#
# The dose-response coefficients are drawn once per band (1,000 draws) and
# reused for every country, so country estimates are positively correlated and
# the sub-Saharan total's interval is the sum's, not a root-sum-of-squares.
# Intervals cover the dose-response only. MAP and IGME end in 2024; 2025 carries
# 2024 forward. The IHME comparison is made at 2024, the last year with actual
# inputs on both sides.
#
# Outputs (results/dhs_rebuild)
#   person_time_burden_ssa_timeseries.csv       country x year, with IHME where available
#   person_time_burden_ssa_2024_comparison.csv  country table against IHME, WHO, region-level pipeline
#   person_time_burden_ssa_totals.csv           sub-Saharan totals by year
#   figure33_person_time_burden_ssa.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "terra", "sf", "ggplot2", "MASS", "countrycode", "patchwork"))

YEARS <- 2005:2025
LAST_DATA_YEAR <- 2024L
COMPARE_YEAR <- 2024L
N_DRAWS <- 1000L
CACHE <- file.path(DATA_DIR, "map_country_series_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
set.seed(20260828)

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
registry <- registry[registry$svkey %in% map$svkey, ]
latest <- do.call(rbind, lapply(split(registry, registry$iso3), function(z) z[which.max(z$year), ]))
person_time <- read.csv(file.path(DERIVED_DIR, "person_time_region_window_segment.csv"), stringsAsFactors = FALSE)
person_time$band <- band_from_segment(person_time$seg_lo, c(1, 4, 12, 24, 36), AGE6B)
national_pfpr <- read.csv(file.path(DATA_DIR, "pfpr_by_country_year.csv"), stringsAsFactors = FALSE)
national_wb <- read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"), stringsAsFactors = FALSE)
comparison_prev <- read.csv(file.path(RESULTS_DIR, "latest_country_burden_comparison.csv"), stringsAsFactors = FALSE)
fits <- readRDS(file.path(DERIVED_DIR, "person_time_model_bundle.rds"))$bands6b_smooth_fits
model_data <- read.csv(file.path(RESULTS_DIR, "person_time_model_data.csv"), stringsAsFactors = FALSE)

# the country set: script 10's burden countries plus any panel DHS country it lacks
countries <- sort(union(comparison_prev$iso3, latest$iso3))
countries <- countries[countries %in% national_pfpr$iso3]
tier <- ifelse(countries %in% latest$iso3, "regional (latest DHS)", "national (no panel DHS)")
names(tier) <- countries
message(length(countries), " countries: ", sum(tier == "regional (latest DHS)"), " regional, ",
        sum(tier != "regional (latest DHS)"), " national")

## ---- UN IGME 2025 death counts by age block, all countries ------------------------------
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
igme_blocks <- function(iso3) {
  z <- igme[igme$iso3 == iso3 & igme$Indicator %in% IGME_BLOCKS, c("year", "Indicator", "OBS_VALUE")]
  if (!nrow(z)) return(NULL)
  w <- reshape(z, idvar = "year", timevar = "Indicator", direction = "wide")
  names(w) <- c("year", names(IGME_BLOCKS)[match(sub("OBS_VALUE.", "", names(w)[-1]), IGME_BLOCKS)])
  w <- w[w$year %in% YEARS, ]
  missing <- setdiff(YEARS, w$year)
  if (length(missing)) { last <- w[w$year == max(w$year), ]; w <- rbind(w, do.call(rbind, lapply(missing, function(y) { l <- last; l$year <- y; l }))) }
  w[order(w$year), ]
}

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

## ---- dose-response draws, once per band ---------------------------------------------------------
band_draws <- lapply(AGE6B, function(band) MASS::mvrnorm(N_DRAWS, coef(fits[[band]]), vcov(fits[[band]])))
names(band_draws) <- AGE6B
af_draws <- function(band, prevalence) {
  fit <- fits[[band]]
  band_data <- model_data[model_data$age6b == band, ]
  segment_levels <- levels(droplevels(factor(band_data$segment)))
  frame <- function(p) {
    out <- data.frame(pfpr10 = p / 10,
                      segment_f = factor(segment_levels[1], levels = segment_levels),
                      window_f = factor("1", levels = as.character(1:5)), year_c = 0, log_pm = 0,
                      country = factor(levels(factor(model_data$iso3))[1], levels = levels(factor(model_data$iso3))),
                      survey = factor(levels(factor(model_data$svkey))[1], levels = levels(factor(model_data$svkey))))
    G <- matrix(0, nrow(out), length(grep("^G", names(coef(fit)))))
    colnames(G) <- sub("^G", "", grep("^G", names(coef(fit)), value = TRUE))
    out$G <- G
    out
  }
  Xh <- predict(fit, frame(prevalence), type = "lpmatrix", discrete = FALSE)
  Xl <- predict(fit, frame(rep(AF_REFERENCE, length(prevalence))), type = "lpmatrix", discrete = FALSE)
  re <- grep("^s\\(country\\)|^s\\(survey\\)", colnames(Xh))
  Xh[, re] <- 0; Xl[, re] <- 0
  pmax(1 - exp(-(band_draws[[band]] %*% t(Xh - Xl))), 0)
}

## ---- pooled DHS age shares for countries without a survey ------------------------------------------
pooled_band <- tapply(person_time$deaths_w, person_time$band, sum)
infant_bands <- AGE6B[2:3]; child_bands <- AGE6B[4:6]
pooled_within_infant <- pooled_band[infant_bands] / sum(pooled_band[infant_bands])
pooled_within_child <- pooled_band[child_bands] / sum(pooled_band[child_bands])

## ---- the burden -----------------------------------------------------------------------------------------
series_rows <- list(); ssa_total <- matrix(0, N_DRAWS, length(YEARS), dimnames = list(NULL, YEARS)); ssa_pn <- ssa_total
for (iso3 in countries) {
  blocks <- igme_blocks(iso3)
  if (is.null(blocks)) { message("  ", iso3, ": no UN IGME counts; skipped"); next }
  if (tier[[iso3]] == "regional (latest DHS)") {
    svkey <- latest$svkey[latest$iso3 == iso3]
    pt <- person_time[person_time$svkey == svkey, ]
    deaths_rb <- tapply(pt$deaths_w, list(pt$regkey, pt$band), sum); deaths_rb[is.na(deaths_rb)] <- 0
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
  share_rb <- sweep(deaths_rb, 2, pmax(colSums(deaths_rb), 1e-9), "/")
  total <- matrix(0, N_DRAWS, length(YEARS), dimnames = list(NULL, YEARS)); pn <- total
  for (band in AGE6B) {
    national_band <- if (band == AGE6B[1]) blocks$neonatal else if (band %in% infant_bands) blocks$months_1_11 * within_infant[[band]] else blocks$years_1_4 * within_child[[band]]
    af <- af_draws(band, as.numeric(prev)); dim(af) <- c(N_DRAWS, length(regions), length(YEARS))
    deaths <- sweep(af, c(2, 3), outer(share_rb[, band], national_band), "*")
    by_year <- apply(deaths, c(1, 3), sum)
    total <- total + by_year
    if (band != AGE6B[1]) pn <- pn + by_year
  }
  ssa_total <- ssa_total + total; ssa_pn <- ssa_pn + pn
  allcause <- blocks$neonatal + blocks$months_1_11 + blocks$years_1_4
  series_rows[[iso3]] <- data.frame(
    iso3 = iso3, tier = tier[[iso3]], year = YEARS, regions = length(regions),
    pfpr_death_weighted = as.numeric(colSums(prev * rowSums(deaths_rb) / sum(deaths_rb))),
    allcause_under5_deaths = allcause,
    malaria_deaths = colMeans(total), lo = apply(total, 2, stats::quantile, 0.025), hi = apply(total, 2, stats::quantile, 0.975),
    malaria_deaths_1_59m = colMeans(pn), lo_1_59m = apply(pn, 2, stats::quantile, 0.025), hi_1_59m = apply(pn, 2, stats::quantile, 0.975),
    stringsAsFactors = FALSE)
  message(sprintf("  %s %-24s %2d regions  %d: %6.0fk (%s)", iso3, tier[[iso3]], length(regions), COMPARE_YEAR,
                  mean(total[, as.character(COMPARE_YEAR)]) / 1e3, "all under-5"))
}
series <- do.call(rbind, series_rows); rownames(series) <- NULL
series$country <- countrycode::countrycode(series$iso3, "iso3c", "country.name", warn = FALSE)

## ---- comparators ------------------------------------------------------------------------------------------
ihme <- read.csv(file.path(DATA_DIR, "ihme_malaria_u5_deaths_by_country_year.csv"), stringsAsFactors = FALSE, check.names = FALSE)
ihme <- ihme[ihme$Measure == "Deaths" & ihme$Age == "Under 5" & ihme$Unit == "Number", ]
ihme$iso3 <- countrycode::countrycode(ihme$Location, "country.name", "iso3c", warn = FALSE)
ihme <- ihme[!is.na(ihme$iso3), c("iso3", "Year", "Value", "Lower", "Upper")]
names(ihme) <- c("iso3", "year", "ihme_u5", "ihme_lo", "ihme_hi")
series <- merge(series, ihme, by = c("iso3", "year"), all.x = TRUE)
series <- series[order(series$iso3, series$year), ]
write.csv(series, file.path(RESULTS_DIR, "person_time_burden_ssa_timeseries.csv"), row.names = FALSE)

at <- series[series$year == COMPARE_YEAR, ]
comparison <- data.frame(
  iso3 = at$iso3, country = at$country, tier = at$tier, regions = at$regions,
  pfpr = round(at$pfpr_death_weighted, 1), allcause_under5 = at$allcause_under5_deaths,
  model_under5 = at$malaria_deaths, model_lo = at$lo, model_hi = at$hi,
  model_1_59m = at$malaria_deaths_1_59m,
  ihme_under5 = at$ihme_u5, ihme_lo = at$ihme_lo, ihme_hi = at$ihme_hi,
  stringsAsFactors = FALSE)
comparison$ratio_model_ihme <- comparison$model_under5 / comparison$ihme_under5
prev_cols <- comparison_prev[, c("iso3", "best_model_deaths", "who_u5_proxy_2024", "who_all_age_2024")]
names(prev_cols) <- c("iso3", "region_pipeline_2024", "who_u5_proxy_2024", "who_all_age_2024")
comparison <- merge(comparison, prev_cols, by = "iso3", all.x = TRUE)
comparison <- comparison[order(-comparison$model_under5), ]
write.csv(comparison, file.path(RESULTS_DIR, "person_time_burden_ssa_2024_comparison.csv"), row.names = FALSE)

totals <- data.frame(
  year = YEARS, countries = length(series_rows),
  malaria_deaths = colMeans(ssa_total), lo = apply(ssa_total, 2, stats::quantile, 0.025), hi = apply(ssa_total, 2, stats::quantile, 0.975),
  malaria_deaths_1_59m = colMeans(ssa_pn), lo_1_59m = apply(ssa_pn, 2, stats::quantile, 0.025), hi_1_59m = apply(ssa_pn, 2, stats::quantile, 0.975),
  allcause_under5 = tapply(series$allcause_under5_deaths, series$year, sum),
  ihme_matched_countries = tapply(series$ihme_u5, series$year, function(v) if (all(is.na(v))) NA else sum(v, na.rm = TRUE)),
  ihme_matched_n = tapply(series$ihme_u5, series$year, function(v) sum(!is.na(v))),
  stringsAsFactors = FALSE)
write.csv(totals, file.path(RESULTS_DIR, "person_time_burden_ssa_totals.csv"), row.names = FALSE)

with_ihme <- comparison[is.finite(comparison$ihme_under5), ]
message(sprintf("\nSub-Saharan total, %d: model %.0fk (%.0fk-%.0fk), 1-59 months %.0fk; IHME (%d matched countries) %.0fk; all-cause under-5 %.0fk",
                COMPARE_YEAR, totals$malaria_deaths[totals$year == COMPARE_YEAR] / 1e3,
                totals$lo[totals$year == COMPARE_YEAR] / 1e3, totals$hi[totals$year == COMPARE_YEAR] / 1e3,
                totals$malaria_deaths_1_59m[totals$year == COMPARE_YEAR] / 1e3,
                nrow(with_ihme), sum(with_ihme$ihme_under5) / 1e3, totals$allcause_under5[totals$year == COMPARE_YEAR] / 1e3))
change <- ssa_total[, as.character(max(YEARS))] / ssa_total[, "2005"] - 1
message(sprintf("2005 -> %d: %.0fk -> %.0fk, %+.0f%% (%+.0f%% to %+.0f%%); IHME matched countries 2005 -> 2024: %+.0f%%",
                max(YEARS), totals$malaria_deaths[1] / 1e3, totals$malaria_deaths[nrow(totals)] / 1e3,
                100 * mean(change), 100 * stats::quantile(change, 0.025), 100 * stats::quantile(change, 0.975),
                100 * (totals$ihme_matched_countries[totals$year == 2024] / totals$ihme_matched_countries[totals$year == 2005] - 1)))
message(sprintf("Ratio model/IHME across countries: median %.2f, IQR %.2f-%.2f; %d of %d countries above IHME",
                stats::median(with_ihme$ratio_model_ihme), stats::quantile(with_ihme$ratio_model_ihme, 0.25),
                stats::quantile(with_ihme$ratio_model_ihme, 0.75), sum(with_ihme$ratio_model_ihme > 1), nrow(with_ihme)))
message("\nCountry table, ", COMPARE_YEAR, " (thousands):")
print(transform(comparison, model = sprintf("%.0f (%.0f-%.0f)", model_under5 / 1e3, model_lo / 1e3, model_hi / 1e3),
                m1_59 = round(model_1_59m / 1e3), ihme = sprintf("%.0f (%.0f-%.0f)", ihme_under5 / 1e3, ihme_lo / 1e3, ihme_hi / 1e3),
                ratio = round(ratio_model_ihme, 2), region_pipeline = round(region_pipeline_2024 / 1e3),
                allcause = round(allcause_under5 / 1e3), tier = substr(tier, 1, 8))[
                  , c("iso3", "tier", "regions", "pfpr", "allcause", "model", "m1_59", "ihme", "ratio", "region_pipeline")],
      row.names = FALSE, right = FALSE)

## ---- figure ----------------------------------------------------------------------------------------------
# countries where IHME puts fewer than 100 under-5 malaria deaths sit at the
# bottom of a log axis and stretch it without adding information
scatter_df <- with_ihme[with_ihme$ihme_under5 >= 100 & with_ihme$model_under5 >= 100, ]
scatter_df$tier_short <- ifelse(grepl("^regional", scatter_df$tier), "Admin-1 (latest DHS)", "National prevalence only")
tiny <- setdiff(with_ihme$iso3, scatter_df$iso3)
scatter <- ggplot2::ggplot(scatter_df, ggplot2::aes(ihme_under5 / 1e3, model_under5 / 1e3)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = model_lo / 1e3, ymax = model_hi / 1e3), colour = "grey70", width = 0) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ihme_lo / 1e3, xmax = ihme_hi / 1e3), colour = "grey70", height = 0) +
  ggplot2::geom_point(ggplot2::aes(shape = tier_short), colour = "#1D6F8B", size = 2.2) +
  ggplot2::geom_text(ggplot2::aes(label = iso3), size = 2.4, vjust = -0.8, colour = "grey25") +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::scale_shape_manual(values = c(16, 1), name = NULL) +
  ggplot2::labs(x = sprintf("IHME/GBD under-5 malaria deaths, %d (thousands, log)", COMPARE_YEAR),
                y = "Model, all under-5 (thousands, log)",
                title = sprintf("Country estimates against IHME, %d", COMPARE_YEAR),
                subtitle = sprintf("Dashed line is equality; bars are 95%% intervals (model: dose-response only).\nNot shown, under 100 deaths on either side: %s",
                                   paste(tiny, collapse = ", "))) +
  ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(legend.position = "bottom")
trend_df <- rbind(
  data.frame(year = totals$year, deaths = totals$malaria_deaths, lo = totals$lo, hi = totals$hi, series = "Person-time model, all under-5"),
  data.frame(year = totals$year, deaths = totals$malaria_deaths_1_59m, lo = NA, hi = NA, series = "Person-time model, 1-59 months"),
  data.frame(year = totals$year, deaths = totals$ihme_matched_countries, lo = NA, hi = NA, series = "IHME/GBD, same countries"))
trend <- ggplot2::ggplot(trend_df, ggplot2::aes(year, deaths / 1e3, colour = series)) +
  ggplot2::geom_ribbon(data = trend_df[!is.na(trend_df$lo), ], ggplot2::aes(ymin = lo / 1e3, ymax = hi / 1e3), fill = "#1D6F8B", alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_vline(xintercept = LAST_DATA_YEAR + 0.5, linetype = "dotted", colour = "grey50") +
  ggplot2::scale_colour_manual(values = c("#1D6F8B", "#7FB3C8", "#B2182B"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(x = NULL, y = "Malaria-attributable under-5 deaths,\nsub-Saharan Africa (thousands)",
                title = sprintf("Sub-Saharan total, %d countries", length(series_rows)),
                subtitle = "Right of the dotted line MAP and IGME inputs are 2024 carried forward") +
  ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(RESULTS_DIR, "figure33_person_time_burden_ssa.png"), scatter | trend,
                width = 13, height = 6, dpi = 200, bg = "white")
message("\nWrote person_time_burden_ssa_timeseries.csv, person_time_burden_ssa_2024_comparison.csv, ",
        "person_time_burden_ssa_totals.csv and figure33_person_time_burden_ssa.png")
