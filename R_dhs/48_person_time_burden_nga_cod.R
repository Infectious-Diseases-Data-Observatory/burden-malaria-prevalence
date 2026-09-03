# =============================================================================
# 48_person_time_burden_nga_cod.R — malaria-attributable under-5 deaths in
# Nigeria and DR Congo from the age-band person-time models, aggregated over
# admin-1 regions, 2005-2025.
#
# For each country the latest DHS in the panel defines the regions (Nigeria
# 2024: six geopolitical zones, which are the DHS admin-1 units in this
# pipeline; DRC 2023: 26 provinces). For each region r, age band a and year y:
#
#   malaria deaths(r, a, y) = AF_a(PfPR_r(y)) x share(r, a) x D_a(y)
#
#   AF_a(p)     attributable fraction from band a's smooth dose-response (script
#               42), 1 - 1/HR against the 0% reference, floored at zero
#   PfPR_r(y)   population-weighted MAP PfPR2-10 for the region in year y
#   share(r, a) the region's share of the country's deaths in band a in the
#               latest survey (weighted, five windows), held fixed over time
#   D_a(y)      national all-cause deaths in band a. UN IGME publishes death
#               counts for three under-5 age blocks (neonatal, 1-11 months,
#               1-4 years); the survey's age-at-death shares split the 1-11
#               month block into 1-3 and 4-11 months and the 1-4 year block
#               into 12-23, 24-35 and 36-59 months. If the IGME file is absent
#               the script falls back to IGME rates x World Bank births with a
#               single DHS split of the post-neonatal total.
#
# The IGME counts come from the UN IGME 2025 "download all" archive on
# childmortality.org (UN-IGME-2025.zip, one CSV of every indicator, sex and
# wealth quintile, with source data). The site sits behind Cloudflare rate
# limiting, so the fetch below uses a browser user-agent and pauses between
# attempts; run it once and the extract is cached in data/derived_dhs.
#
# Uncertainty is the dose-response model's only: 1,000 draws of each band's
# coefficients from their fitted normal approximation. MAP surfaces and the IGME
# series end in 2024; 2025 carries 2024 forward and is labelled as such.
#
# Outputs (results/dhs_rebuild)
#   person_time_burden_nga_cod_timeseries.csv
#   person_time_burden_nga_cod_2025.csv          by band and by region
#   figure32_person_time_burden_nga_cod.png
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("mgcv", "terra", "sf", "ggplot2", "MASS"))

COUNTRIES <- c(NGA = "NG8BFL", COD = "CD81FL")
YEARS <- 2005:2025
LAST_DATA_YEAR <- 2024L
N_DRAWS <- 1000L
CACHE <- file.path(DATA_DIR, "map_country_series_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
set.seed(20260828)

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
person_time <- read.csv(file.path(DERIVED_DIR, "person_time_region_window_segment.csv"),
                        stringsAsFactors = FALSE)
person_time$band <- band_from_segment(person_time$seg_lo, c(1, 4, 12, 24, 36), AGE6B)
national <- read.csv(file.path(DATA_DIR, "wb_mortality_timeseries.csv"), stringsAsFactors = FALSE)

## ---- UN IGME 2025 death counts by age block ----------------------------------------------
IGME_ZIP_URL <- "https://childmortality.org/wp-content/uploads/UN-IGME-2025.zip"
IGME_DIR <- file.path(DATA_DIR, "igme")
IGME_EXTRACT <- file.path(DERIVED_DIR, "igme_2025_nga_cod.csv")
IGME_BLOCKS <- c(neonatal = "Neonatal deaths", months_1_11 = "Deaths age 1-11 months",
                 years_1_4 = "Child deaths age 1 to 4")
fetch_igme <- function() {
  dir.create(IGME_DIR, showWarnings = FALSE, recursive = TRUE)
  zip <- file.path(IGME_DIR, "UN-IGME-2025.zip")
  if (!file.exists(zip)) {
    required_packages("httr")
    for (attempt in 1:4) {
      r <- tryCatch(httr::GET(IGME_ZIP_URL, httr::user_agent("Mozilla/5.0"),
                              httr::write_disk(zip, overwrite = TRUE), httr::timeout(900)),
                    error = function(e) NULL)
      if (!is.null(r) && httr::status_code(r) == 200L) break
      unlink(zip); message("  UN IGME download attempt ", attempt, " refused (rate limit); pausing")
      Sys.sleep(75)
    }
    if (!file.exists(zip)) stop("Could not download ", IGME_ZIP_URL)
  }
  csv <- unzip(zip, exdir = IGME_DIR, overwrite = TRUE)
  x <- read.csv(csv[grepl("[.]csv$", csv)][1], stringsAsFactors = FALSE, check.names = TRUE)
  est <- x[x$REF_AREA %in% names(COUNTRIES) & x$Sex == "Total" & x$Wealth.Quintile == "Total" &
             x$Series.Name == "UN IGME estimate", ]
  est$year <- floor(est$REF_DATE)
  out <- est[, c("REF_AREA", "year", "Indicator", "OBS_VALUE", "LOWER_BOUND", "UPPER_BOUND")]
  names(out)[1] <- "iso3"
  write.csv(out, IGME_EXTRACT, row.names = FALSE)
  out
}
igme <- tryCatch({
  if (file.exists(IGME_EXTRACT)) read.csv(IGME_EXTRACT, stringsAsFactors = FALSE) else fetch_igme()
}, error = function(e) { message("UN IGME counts unavailable (", conditionMessage(e), "); using rates x births"); NULL })
igme_blocks <- function(iso3) {
  if (is.null(igme)) return(NULL)
  z <- igme[igme$iso3 == iso3 & igme$Indicator %in% IGME_BLOCKS, c("year", "Indicator", "OBS_VALUE")]
  if (!nrow(z)) return(NULL)
  w <- reshape(z, idvar = "year", timevar = "Indicator", direction = "wide")
  names(w) <- c("year", names(IGME_BLOCKS)[match(sub("OBS_VALUE.", "", names(w)[-1]), IGME_BLOCKS)])
  w[order(w$year), ]
}
fits <- readRDS(file.path(DERIVED_DIR, "person_time_model_bundle.rds"))$bands6b_smooth_fits
model_data <- read.csv(file.path(RESULTS_DIR, "person_time_model_data.csv"), stringsAsFactors = FALSE)

## ---- MAP prevalence by region and year --------------------------------------------------
raster_files <- list.files(MAP_RASTER_DIR, pattern = "pfpr2_10_[0-9]{4}[.]tif$", full.names = TRUE)
raster_years <- as.integer(regmatches(basename(raster_files), regexpr("[0-9]{4}", basename(raster_files))))
rasters <- setNames(lapply(raster_files, terra::rast), as.character(raster_years))
weights_full <- terra::resample(terra::crop(terra::rast(GPW_TIF), rasters[[1]]), rasters[[1]],
                                method = "bilinear")
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
    pfpr = 100 * numerator[[j]] / denominator, population_weight = denominator,
    stringsAsFactors = FALSE)))
  out <- out[nzchar(out$regkey) & is.finite(out$pfpr), ]
  saveRDS(out, cache)
  out
}

## ---- attributable fraction draws from a band's smooth fit --------------------------------
af_draws <- function(fit, prevalence, band_data) {
  segment_levels <- levels(droplevels(factor(band_data$segment)))
  frame <- function(p) {
    out <- data.frame(pfpr10 = p / 10,
                      segment_f = factor(segment_levels[1], levels = segment_levels),
                      window_f = factor("1", levels = as.character(1:5)),
                      year_c = 0, log_pm = 0,
                      country = factor(levels(factor(model_data$iso3))[1], levels = levels(factor(model_data$iso3))),
                      survey = factor(levels(factor(model_data$svkey))[1], levels = levels(factor(model_data$svkey))))
    G <- matrix(0, nrow(out), length(grep("^G", names(coef(fit)))))
    colnames(G) <- sub("^G", "", grep("^G", names(coef(fit)), value = TRUE))
    out$G <- G
    out
  }
  Xh <- predict(fit, frame(prevalence), type = "lpmatrix", discrete = FALSE)
  Xl <- predict(fit, frame(rep(PERSON_TIME_AF_REFERENCE, length(prevalence))), type = "lpmatrix", discrete = FALSE)
  re <- grep("^s\\(country\\)|^s\\(survey\\)", colnames(Xh))
  Xh[, re] <- 0; Xl[, re] <- 0
  dX <- Xh - Xl
  beta <- MASS::mvrnorm(N_DRAWS, coef(fit), vcov(fit))       # draws x coefficients
  log_hr <- beta %*% t(dX)                                    # draws x prevalence values
  pmax(1 - exp(-log_hr), 0)
}

## ---- the burden ----------------------------------------------------------------------------------
series_rows <- list(); band_rows <- list(); region_rows <- list()
for (iso3 in names(COUNTRIES)) {
  svkey <- COUNTRIES[[iso3]]
  message("\n", iso3, ": ", svkey)
  prevalence <- regional_prevalence(svkey)
  pt <- person_time[person_time$svkey == svkey, ]
  stopifnot(all(unique(pt$regkey) %in% unique(prevalence$regkey)))
  regions <- sort(unique(pt$regkey))
  message("  ", length(regions), " regions; ", length(unique(prevalence$year)), " years of MAP prevalence")

  # regional share of deaths within each band, and the survey's age-at-death shares
  # within the two IGME post-neonatal blocks
  deaths_rb <- tapply(pt$deaths_w, list(pt$regkey, pt$band), sum)
  deaths_rb[is.na(deaths_rb)] <- 0
  share_rb <- sweep(deaths_rb, 2, colSums(deaths_rb), "/")
  band_totals <- colSums(deaths_rb)
  infant_bands <- AGE6B[2:3]; child_bands <- AGE6B[4:6]
  within_infant <- band_totals[infant_bands] / sum(band_totals[infant_bands])
  within_child <- band_totals[child_bands] / sum(band_totals[child_bands])
  pn_share <- band_totals[AGE6B[-1]] / sum(band_totals[AGE6B[-1]])

  nat <- national[national$iso3 == iso3 & national$year %in% YEARS, ]
  carried <- setdiff(YEARS, nat$year)
  if (length(carried)) {
    last <- nat[nat$year == max(nat$year), ]
    nat <- rbind(nat, do.call(rbind, lapply(carried, function(y) { z <- last; z$year <- y; z })))
  }
  nat <- nat[order(nat$year), ]
  blocks <- igme_blocks(iso3)
  if (!is.null(blocks)) {
    blocks <- blocks[blocks$year %in% YEARS, ]
    carried_igme <- setdiff(YEARS, blocks$year)
    if (length(carried_igme)) {
      last <- blocks[blocks$year == max(blocks$year), ]
      blocks <- rbind(blocks, do.call(rbind, lapply(carried_igme, function(y) { z <- last; z$year <- y; z })))
    }
    blocks <- blocks[order(blocks$year), ]
    stopifnot(identical(blocks$year, nat$year))
    nat$neonatal_deaths <- blocks$neonatal
    nat$infant_1_11_deaths <- blocks$months_1_11
    nat$child_1_4_deaths <- blocks$years_1_4
    nat$postneonatal_deaths <- nat$infant_1_11_deaths + nat$child_1_4_deaths
    band_deaths <- function(band) {
      if (band == AGE6B[1]) nat$neonatal_deaths
      else if (band %in% infant_bands) nat$infant_1_11_deaths * within_infant[[band]]
      else nat$child_1_4_deaths * within_child[[band]]
    }
    message("  all-cause deaths from UN IGME 2025 counts (three age blocks), DHS shares within blocks")
  } else {
    nat$neonatal_deaths <- nat$nnmr / 1000 * nat$births
    nat$postneonatal_deaths <- (nat$u5mr - nat$nnmr) / 1000 * nat$births
    band_deaths <- function(band) {
      if (band == AGE6B[1]) nat$neonatal_deaths else nat$postneonatal_deaths * pn_share[[band]]
    }
    message("  all-cause deaths from IGME rates x World Bank births, DHS shares within 1-59 months")
  }

  # prevalence by region and year, 2025 carrying 2024
  prev_wide <- tapply(prevalence$pfpr, list(prevalence$regkey, prevalence$year), mean)
  prev_wide <- prev_wide[regions, , drop = FALSE]
  for (y in YEARS[YEARS > LAST_DATA_YEAR]) prev_wide <- cbind(prev_wide, setNames(prev_wide[, as.character(LAST_DATA_YEAR), drop = FALSE], y))
  colnames(prev_wide)[colnames(prev_wide) == "2024" & FALSE] <- NA
  colnames(prev_wide) <- c(sort(unique(prevalence$year)), YEARS[YEARS > LAST_DATA_YEAR])

  total_draws <- matrix(0, N_DRAWS, length(YEARS), dimnames = list(NULL, YEARS))
  postneonatal_draws <- total_draws
  for (band in AGE6B) {
    band_data <- model_data[model_data$age6b == band, ]
    # AF draws for every region-year prevalence at once
    p_vec <- as.numeric(prev_wide)                        # regions x years, column-major
    af <- af_draws(fits[[band]], p_vec, band_data)        # draws x (regions*years)
    dim(af) <- c(N_DRAWS, length(regions), length(YEARS))
    national_band <- band_deaths(band)
    # deaths(draw, region, year) = AF x share(region, band) x D_band(year)
    weights <- outer(share_rb[regions, band], national_band)   # regions x years
    deaths <- sweep(af, c(2, 3), weights, "*")
    band_year <- apply(deaths, c(1, 3), sum)                    # draws x years
    total_draws <- total_draws + band_year
    if (band != AGE6B[1]) postneonatal_draws <- postneonatal_draws + band_year
    band_rows[[paste(iso3, band)]] <- data.frame(
      iso3 = iso3, year = max(YEARS), age_band = band,
      allcause_deaths = national_band[nat$year == max(YEARS)],
      malaria_deaths = mean(band_year[, ncol(band_year)]),
      lo = unname(stats::quantile(band_year[, ncol(band_year)], 0.025)),
      hi = unname(stats::quantile(band_year[, ncol(band_year)], 0.975)),
      stringsAsFactors = FALSE)
    if (band == AGE6B[1]) region_acc <- array(0, c(N_DRAWS, length(regions)))
    region_acc <- region_acc + deaths[, , length(YEARS)]
  }
  for (j in seq_along(regions)) {
    region_rows[[paste(iso3, regions[j])]] <- data.frame(
      iso3 = iso3, year = max(YEARS), regkey = regions[j],
      pfpr = prev_wide[j, as.character(max(YEARS))],
      share_of_under5_deaths = sum(deaths_rb[regions[j], ]) / sum(deaths_rb),
      malaria_deaths = mean(region_acc[, j]),
      lo = unname(stats::quantile(region_acc[, j], 0.025)),
      hi = unname(stats::quantile(region_acc[, j], 0.975)),
      stringsAsFactors = FALSE)
  }
  pfpr_national <- colSums(prev_wide * as.numeric(rowSums(deaths_rb) / sum(deaths_rb)))
  series_rows[[iso3]] <- data.frame(
    iso3 = iso3, year = YEARS,
    pfpr_death_weighted = as.numeric(pfpr_national),
    allcause_under5_deaths = nat$neonatal_deaths + nat$postneonatal_deaths,
    allcause_postneonatal_deaths = nat$postneonatal_deaths,
    malaria_deaths = colMeans(total_draws),
    lo = apply(total_draws, 2, stats::quantile, 0.025),
    hi = apply(total_draws, 2, stats::quantile, 0.975),
    malaria_deaths_1_59m = colMeans(postneonatal_draws),
    lo_1_59m = apply(postneonatal_draws, 2, stats::quantile, 0.025),
    hi_1_59m = apply(postneonatal_draws, 2, stats::quantile, 0.975),
    carried_forward = YEARS > LAST_DATA_YEAR,
    stringsAsFactors = FALSE)
  change <- total_draws[, ncol(total_draws)] / total_draws[, 1] - 1
  message(sprintf("  %d: %.0fk malaria-attributable under-5 deaths (%.0fk-%.0fk), of which 1-59 months %.0fk; all-cause under-5 %.0fk; AF %.1f%%",
                  max(YEARS), mean(total_draws[, ncol(total_draws)]) / 1e3,
                  stats::quantile(total_draws[, ncol(total_draws)], 0.025) / 1e3,
                  stats::quantile(total_draws[, ncol(total_draws)], 0.975) / 1e3,
                  mean(postneonatal_draws[, ncol(postneonatal_draws)]) / 1e3,
                  (nat$neonatal_deaths + nat$postneonatal_deaths)[nat$year == max(YEARS)] / 1e3,
                  100 * mean(total_draws[, ncol(total_draws)]) / (nat$neonatal_deaths + nat$postneonatal_deaths)[nat$year == max(YEARS)]))
  message(sprintf("  %d -> %d: %.0fk -> %.0fk, change %+.0f%% (%+.0f%% to %+.0f%%); death-weighted PfPR %.1f%% -> %.1f%%",
                  min(YEARS), max(YEARS), mean(total_draws[, 1]) / 1e3, mean(total_draws[, ncol(total_draws)]) / 1e3,
                  100 * mean(change), 100 * stats::quantile(change, 0.025), 100 * stats::quantile(change, 0.975),
                  pfpr_national[1], pfpr_national[length(pfpr_national)]))
}
series <- do.call(rbind, series_rows); rownames(series) <- NULL
by_band <- do.call(rbind, band_rows); rownames(by_band) <- NULL
by_region <- do.call(rbind, region_rows); rownames(by_region) <- NULL
write.csv(series, file.path(RESULTS_DIR, "person_time_burden_nga_cod_timeseries.csv"), row.names = FALSE)
write.csv(rbind(cbind(level = "band", by_band[, c("iso3", "year", "age_band", "allcause_deaths", "malaria_deaths", "lo", "hi")]),
                cbind(level = "region", data.frame(iso3 = by_region$iso3, year = by_region$year,
                                                   age_band = by_region$regkey, allcause_deaths = NA,
                                                   malaria_deaths = by_region$malaria_deaths,
                                                   lo = by_region$lo, hi = by_region$hi))),
          file.path(RESULTS_DIR, "person_time_burden_nga_cod_2025.csv"), row.names = FALSE)

message("\n", max(YEARS), " by age band:")
print(transform(by_band, allcause_deaths = round(allcause_deaths), malaria_deaths = round(malaria_deaths),
                lo = round(lo), hi = round(hi), af_pct = round(100 * malaria_deaths / allcause_deaths, 1)),
      row.names = FALSE)
message("\n", max(YEARS), " by region:")
print(transform(by_region, pfpr = round(pfpr, 1), share_of_under5_deaths = round(share_of_under5_deaths, 3),
                malaria_deaths = round(malaria_deaths), lo = round(lo), hi = round(hi)), row.names = FALSE)

## ---- figure ----------------------------------------------------------------------------------------
labels <- c(NGA = "Nigeria (6 zones, DHS 2024)", COD = "DR Congo (26 provinces, DHS 2023)")
series$country <- labels[series$iso3]
plot_deaths <- ggplot2::ggplot(series, ggplot2::aes(year, malaria_deaths / 1e3)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo / 1e3, ymax = hi / 1e3), fill = "#1D6F8B", alpha = 0.18) +
  ggplot2::geom_line(colour = "#1D6F8B", linewidth = 1) +
  ggplot2::geom_line(ggplot2::aes(y = malaria_deaths_1_59m / 1e3), colour = "#1D6F8B", linetype = "dashed") +
  ggplot2::geom_vline(xintercept = LAST_DATA_YEAR + 0.5, linetype = "dotted", colour = "grey50") +
  ggplot2::facet_wrap(~country, scales = "free_y") +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(x = NULL, y = "Malaria-attributable under-5 deaths\n(thousands; dashed = ages 1-59 months)",
                title = "Malaria-attributable under-5 deaths from the age-band person-time models",
                subtitle = paste("Admin-1 prevalence x regional death shares x UN IGME 2025 all-cause deaths by age block",
                                 "(DHS shares within blocks);\nband = 95% interval from the dose-response fits only.",
                                 "Right of the dotted line MAP and IGME inputs are 2024 carried forward.")) +
  ggplot2::theme_minimal(base_size = 10)
plot_context <- ggplot2::ggplot(series, ggplot2::aes(year)) +
  ggplot2::geom_line(ggplot2::aes(y = allcause_under5_deaths / 1e3, colour = "All-cause under-5 deaths (thousands)"), linewidth = 0.9) +
  ggplot2::geom_line(ggplot2::aes(y = pfpr_death_weighted * 10, colour = "PfPR2-10, death-weighted mean (% x 10)"), linewidth = 0.9) +
  ggplot2::geom_vline(xintercept = LAST_DATA_YEAR + 0.5, linetype = "dotted", colour = "grey50") +
  ggplot2::facet_wrap(~country, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c("grey35", "#D95F0E"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::labs(x = NULL, y = "Inputs") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_blank())
required_packages("patchwork")
combined <- plot_deaths / plot_context + patchwork::plot_layout(heights = c(1.4, 1))
ggplot2::ggsave(file.path(RESULTS_DIR, "figure32_person_time_burden_nga_cod.png"), combined,
                width = 10, height = 6.8, dpi = 200, bg = "white")
message("\nWrote person_time_burden_nga_cod_timeseries.csv, person_time_burden_nga_cod_2025.csv and ",
        "figure32_person_time_burden_nga_cod.png")
