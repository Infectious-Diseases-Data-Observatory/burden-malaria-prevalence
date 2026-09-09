# =============================================================================
# 51_snow_polygon_prevalence.R
# Merge the Snow et al. (2017) polygon PfPR2-10 estimates (Stan re-fit, see
# ../Prevalence Model/R code) onto the DHS/MIS survey regions used throughout
# this pipeline.
#
# For every survey in the registry, each region polygon (data/dhs_boundaries/
# <SurveyId>.rds, DHSREGEN) is intersected with the 520 Snow administrative
# polygons. Each intersection piece is weighted by its 2020 GPW population
# (the same raster script 02 uses for MAP), and the region's PfPR2-10 in each
# period is the population-weighted mean of the overlapping Snow polygons.
# The weighting is applied to the thinned posterior draws, so the region-level
# median and 95% credible interval propagate the Snow model's uncertainty.
#
# Inputs
#   data/derived_dhs/survey_registry.csv, data/dhs_boundaries/<SurveyId>.rds
#   data/pop/gpw_v4_population_density_rev11_2020_2.5m.tif
#   Two sources of polygon estimates, chosen with SNOW_SOURCE:
#
#   SNOW_SOURCE = "draws" (default) - the five-year fits, with posterior draws:
#     <SNOW_OUTPUT_DIR>/polygons_520.gpkg, PfPR210_estimates_*.csv,
#     p_draws_thinned.rds. SNOW_OUTPUT_DIR defaults to
#     ../Prevalence Model/R code/output (1990-2015, 5 periods); point it at
#     output_1900_2015 for the paper's 16 periods. Region intervals are exact,
#     because the population weights are applied to the draws themselves.
#
#   SNOW_SOURCE = "annual" - the ANNUAL fit. Uses p_draws_annual.rds from
#     <SNOW_ANNUAL_DIR> when present, so region intervals are EXACT, and falls
#     back to the summary means when it is absent.
#
#   SNOW_SOURCE = "annual_csv" - the ANNUAL fit, forced to summary means only:
#     <SNOW_ANNUAL_DIR>/snow_pfpr_by_polygon_year.csv and polygons_520.gpkg
#     (SNOW_ANNUAL_DIR defaults to data/snow_prevalence_model). Time steps are
#     single years, and each survey is matched to its own year.
#
#     With means rather than draws, the region POINT ESTIMATE is still exact:
#     a population-weighted average is linear, so sum_j w_j E[p_j] = E[sum_j w_j p_j].
#     The region SD is NOT recoverable, because it needs the covariance between
#     polygons. Two SD summaries are reported:
#       snow_pfpr_sd_indep  = sqrt(sum w_j^2 sd_j^2)  - independence benchmark;
#                             a lower bound only if covariances are nonnegative
#       snow_pfpr_sd_upper  = sum w_j sd_j            - assumes perfect
#                             correlation (upper bound)
#     Use the draws path when correct region intervals matter.
#   data/derived_dhs/map_pfpr_by_survey_region.csv          (comparison only)
#   data/dhs_prevalence_by_region.csv, results/rdt_microscopy_conversion.csv
#
# Outputs (file names gain a suffix: "_1900_2015" for the full five-year fit,
# "_annual" for the annual draws source, "_annual_means" for the means-only one)
#   data/derived_dhs/snow_pfpr_by_survey_region_long.csv   one row per survey-region-period
#   data/derived_dhs/snow_pfpr_by_survey_region.csv        one row per survey-region, wide by
#                                                           period + the survey-year-matched period
#   results/dhs_rebuild/snow_extraction_status.csv          per-survey status
#   results/dhs_rebuild/snow_vs_map_agreement.csv           agreement with MAP, by period
#   results/dhs_rebuild/snow_vs_map_scatter.png
#   results/dhs_rebuild/snow_vs_measured_scatter.png        agreement with DHS-measured parasitaemia
#   results/dhs_rebuild/snow_vs_measured_agreement.csv, snow_vs_measured_by_survey.csv
# =============================================================================
source("R_dhs/00_config.R")
required_packages(c("sf", "terra", "ggplot2", "digest", "exactextractr"))
suppressPackageStartupMessages({ library(sf); library(ggplot2) })
sf::sf_use_s2(FALSE)

SNOW_SOURCE <- Sys.getenv("SNOW_SOURCE", "draws")
if (!SNOW_SOURCE %in% c("draws", "annual", "annual_csv"))
  stop("SNOW_SOURCE must be 'draws', 'annual' or 'annual_csv'")
SNOW_OUTPUT_DIR <- Sys.getenv(
  "SNOW_OUTPUT_DIR",
  normalizePath(file.path(REPO_ROOT, "..", "Prevalence Model", "R code", "output"),
                mustWork = FALSE)
)
SNOW_ANNUAL_DIR <- Sys.getenv("SNOW_ANNUAL_DIR",
                              file.path(DATA_DIR, "snow_prevalence_model"))
## The cache holds computed prevalence values, not just geometry, so a means-based
## run and a draws-based run must never share it.
cache_key <- switch(SNOW_SOURCE, "annual" = "annual_draws",
                    "annual_csv" = "annual_means", basename(SNOW_OUTPUT_DIR))
SNOW_REGION_DIR <- file.path(DATA_DIR, "snow_region_cache", cache_key)
dir.create(SNOW_REGION_DIR, recursive = TRUE, showWarnings = FALSE)
## outputs from the 1990-2015 fit (folder "output") keep plain names; any other
## source gets a suffix ("_1900_2015", "_annual")
## the two annual paths get distinct suffixes so their outputs can coexist and
## be compared rather than overwriting one another
SNOW_SUFFIX     <- switch(SNOW_SOURCE, "annual" = "_annual",
                          "annual_csv" = "_annual_means",
                          sub("^output", "", basename(SNOW_OUTPUT_DIR)))
snow_out <- function(dir, stem, ext) file.path(dir, paste0(stem, SNOW_SUFFIX, ext))
SNOW_LONG_CSV   <- snow_out(DERIVED_DIR, "snow_pfpr_by_survey_region_long", ".csv")
SNOW_WIDE_CSV   <- snow_out(DERIVED_DIR, "snow_pfpr_by_survey_region", ".csv")
SNOW_STATUS_CSV <- snow_out(RESULTS_DIR, "snow_extraction_status", ".csv")
MIN_COVERAGE    <- 0.5   # flag regions with < 50% of population inside the Snow polygons

## ---- Snow polygon estimates -------------------------------------------------
## Both sources are reduced to a common representation:
##   periods / period_start / period_end   time labels and their year ranges
##   have_draws                            TRUE if draws are available
##   draws$p[draw, polygon, period]        (draws source only)
##   M$mean, M$sd, M$median, M$q2.5, M$q97.5   polygon x period matrices, per cent
##                                         (annual_csv source only)
if (SNOW_SOURCE == "draws") {
  snow_poly_file <- file.path(SNOW_OUTPUT_DIR, "polygons_520.gpkg")
  snow_draw_file <- file.path(SNOW_OUTPUT_DIR, "p_draws_thinned.rds")
  for (f in c(snow_poly_file, snow_draw_file, GPW_TIF)) {
    if (!file.exists(f)) stop("Required input missing: ", f)
  }
  snow_poly <- sf::st_make_valid(sf::st_read(snow_poly_file, quiet = TRUE))
  draws     <- readRDS(snow_draw_file)       # $Admin_ID, $periods, $p[draw, polygon, period]
  stopifnot(identical(as.character(draws$Admin_ID), as.character(snow_poly$Admin_ID)))
  periods      <- draws$periods
  period_start <- as.integer(substr(periods, 1, 4))
  period_end   <- as.integer(substr(periods, 6, 9))
  have_draws   <- TRUE
  n_draw       <- dim(draws$p)[1]
  message(sprintf("Snow estimates: %d polygons x %d periods (%s to %s), %d posterior draws, from %s",
                  nrow(snow_poly), length(periods), periods[1], tail(periods, 1),
                  n_draw, SNOW_OUTPUT_DIR))
} else if (SNOW_SOURCE == "annual" &&
           file.exists(file.path(SNOW_ANNUAL_DIR, "p_draws_annual.rds"))) {
  ## annual fit WITH draws: region intervals are exact
  snow_poly_file   <- file.path(SNOW_ANNUAL_DIR, "polygons_520.gpkg")
  annual_draw_file <- file.path(SNOW_ANNUAL_DIR, "p_draws_annual.rds")
  for (f in c(snow_poly_file, annual_draw_file, GPW_TIF))
    if (!file.exists(f)) stop("Required input missing: ", f)
  snow_poly <- sf::st_make_valid(sf::st_read(snow_poly_file, quiet = TRUE))
  draws     <- readRDS(annual_draw_file)
  stopifnot(identical(as.character(draws$Admin_ID), as.character(snow_poly$Admin_ID)),
            identical(draws$scale, "proportion"))
  periods      <- as.character(draws$years)
  period_start <- as.integer(draws$period_start)
  period_end   <- as.integer(draws$period_end)
  have_draws   <- TRUE
  n_draw       <- dim(draws$p)[1]
  message(sprintf("Snow estimates: %d polygons x %d years (%s to %s), %d posterior draws (fit '%s'), from %s",
                  nrow(snow_poly), length(periods), periods[1], tail(periods, 1),
                  n_draw, draws$fit_tag %||% "unknown", annual_draw_file))
  message("Region intervals are exact: the population weights are applied to the draws.")
} else {
  if (SNOW_SOURCE == "annual")
    message("p_draws_annual.rds not found in ", SNOW_ANNUAL_DIR,
            "; falling back to summary means without region quantiles.")
  snow_poly_file <- file.path(SNOW_ANNUAL_DIR, "polygons_520.gpkg")
  annual_csv     <- file.path(SNOW_ANNUAL_DIR, "snow_pfpr_by_polygon_year.csv")
  for (f in c(snow_poly_file, annual_csv, GPW_TIF))
    if (!file.exists(f)) stop("Required input missing: ", f)
  snow_poly <- sf::st_make_valid(sf::st_read(snow_poly_file, quiet = TRUE))
  ann <- read.csv(annual_csv, stringsAsFactors = FALSE)
  need <- c("Admin_ID", "year", "PfPR_mean", "PfPR_sd", "PfPR_median",
            "PfPR_q2.5", "PfPR_q97.5")
  miss <- setdiff(need, names(ann))
  if (length(miss)) stop("annual CSV lacks column(s): ", paste(miss, collapse = ", "))
  yrs <- sort(unique(ann$year))
  stopifnot(!anyDuplicated(ann[c("Admin_ID", "year")]),
            !anyDuplicated(snow_poly$Admin_ID),
            setequal(ann$Admin_ID, snow_poly$Admin_ID),
            identical(as.integer(yrs), 2000:2015),
            all(is.finite(ann$PfPR_mean) & ann$PfPR_mean >= 0 & ann$PfPR_mean <= 100),
            all(is.finite(ann$PfPR_sd) & ann$PfPR_sd >= 0),
            all(ann$Country_ID == snow_poly$Country_ID[match(ann$Admin_ID, snow_poly$Admin_ID)]))
  ## polygon x year matrices, rows aligned to snow_poly$Admin_ID
  to_mat <- function(col) {
    m <- matrix(NA_real_, nrow(snow_poly), length(yrs),
                dimnames = list(as.character(snow_poly$Admin_ID), as.character(yrs)))
    ri <- match(ann$Admin_ID, snow_poly$Admin_ID)
    ci <- match(ann$year, yrs)
    ok <- !is.na(ri) & !is.na(ci)
    m[cbind(ri[ok], ci[ok])] <- ann[[col]][ok]
    m
  }
  M <- list(mean = to_mat("PfPR_mean"), sd = to_mat("PfPR_sd"),
            median = to_mat("PfPR_median"), q2.5 = to_mat("PfPR_q2.5"),
            q97.5 = to_mat("PfPR_q97.5"))
  if (anyNA(M$mean))
    stop(sprintf("annual CSV is missing %d of the %d polygon-year combinations",
                 sum(is.na(M$mean)), length(M$mean)))
  periods      <- as.character(yrs)
  period_start <- yrs
  period_end   <- yrs
  have_draws   <- FALSE
  n_draw       <- NA_integer_
  message(sprintf("Snow estimates: %d polygons x %d years (%d to %d) from summary means, %s",
                  nrow(snow_poly), length(yrs), min(yrs), max(yrs), annual_csv))
  message("Region posterior means are exact for the fixed weights. Independence SD is a benchmark, ",
          "not a guaranteed lower bound; no region quantiles are inferred from summary CSVs.")
}

## Row order of the draws array / summary matrices, shared by both sources
snow_admin_id <- as.character(snow_poly$Admin_ID)
if (have_draws) {
  stopifnot(identical(snow_admin_id, as.character(draws$Admin_ID)))
} else {
  stopifnot(identical(snow_admin_id, rownames(M$mean)))
}
## the five-year draws are stored without explicit year bounds, so they are
## still parsed from the period labels above; the annual export carries them
stopifnot(!anyNA(period_start), !anyNA(period_end))

## population count per raster cell (density is persons per km^2)
density  <- terra::rast(GPW_TIF)
density  <- terra::crop(density, terra::ext(-20, 55, -37, 40))
pop_rast <- density * terra::cellSize(density, unit = "km")
source_files <- c(snow_poly_file, GPW_TIF, "R_dhs/51_snow_polygon_prevalence.R",
                 "R_dhs/00_config.R", if (have_draws) {
                   if (SNOW_SOURCE == "draws") snow_draw_file else annual_draw_file
                 } else annual_csv)
source_hashes <- vapply(source_files, digest::digest, "", file = TRUE, algo = "md5")
source_signature <- digest::digest(list(source_hashes, SNOW_SOURCE, have_draws,
  as.character(packageVersion("sf")), as.character(packageVersion("terra")),
  as.character(packageVersion("exactextractr"))), algo = "md5")

## ---- per-survey extraction ---------------------------------------------------
## Map a survey year to a time index. With five-year periods this finds the
## containing block; with annual steps period_start == period_end, so it is an
## exact year match. No carry-forward outside the available exposure period.
period_for_year <- function(y) {
  idx <- which(y >= period_start & y <= period_end)
  if (length(idx)) idx[1] else NA_integer_
}

extract_survey_regions <- function(survey) {
  cache_file <- file.path(SNOW_REGION_DIR, paste0(survey$svkey, ".rds"))
  if (!file.exists(survey$boundary_file)) stop("boundary unavailable")
  cache_signature <- digest::digest(list(source_signature, survey,
    digest::digest(file = survey$boundary_file, algo = "md5")), algo = "md5")
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(attr(cached, "signature"), cache_signature)) return(cached)
  }
  boundary <- readRDS(survey$boundary_file)
  if (!"DHSREGEN" %in% names(boundary)) stop("boundary lacks DHSREGEN")
  boundary <- sf::st_make_valid(sf::st_transform(boundary, sf::st_crs(snow_poly)))
  stopifnot(!anyDuplicated(rkey(boundary$DHSREGEN)))
  boundary$row_id <- seq_len(nrow(boundary))

  ## population of the whole region (denominator for coverage)
  region_pop <- exactextractr::exact_extract(pop_rast, boundary, "sum", progress = FALSE)

  ## intersect with the Snow polygons and weight each piece by population
  pieces <- suppressWarnings(sf::st_intersection(
    boundary[, c("row_id", "DHSREGEN")],
    snow_poly[snow_poly$Country_ID == survey$iso3, "Admin_ID"]
  ))
  pieces <- pieces[sf::st_geometry_type(pieces) %in%
                     c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"), ]
  if (any(sf::st_geometry_type(pieces) == "GEOMETRYCOLLECTION")) {
    pieces <- sf::st_collection_extract(pieces, "POLYGON")
  }
  if (!nrow(pieces)) stop("no overlap with the Snow polygons")
  ## zero-area slivers receive zero population below and drop out of the weights
  pieces$pop <- exactextractr::exact_extract(pop_rast, pieces, "sum", progress = FALSE)
  pieces$pop[!is.finite(pieces$pop)] <- 0
  pieces <- sf::st_drop_geometry(pieces)
  pieces$k <- match(pieces$Admin_ID, snow_admin_id)
  stopifnot(!anyNA(pieces$k))
  saveRDS(list(signature = cache_signature, pieces = pieces, region_population = region_pop,
              region = as.character(boundary$DHSREGEN)),
          file.path(SNOW_REGION_DIR, paste0(survey$svkey, "_weights.rds")))

  rows <- vector("list", nrow(boundary))
  for (r in seq_len(nrow(boundary))) {
    pr <- pieces[pieces$row_id == r & pieces$pop > 0, , drop = FALSE]
    if (!nrow(pr)) next
    w <- pr$pop / sum(pr$pop)
    common <- data.frame(
      svkey = survey$svkey, SurveyId = survey$SurveyId, iso3 = survey$iso3, year = survey$year,
      region = as.character(boundary$DHSREGEN[r]), regkey = rkey(boundary$DHSREGEN[r]),
      period = periods, period_index = seq_along(periods),
      population_weight = sum(pr$pop), region_population = region_pop[r],
      coverage = sum(pr$pop) / region_pop[r], n_snow_polygons = nrow(pr),
      main_snow_polygon = pr$Admin_ID[which.max(pr$pop)],
      stringsAsFactors = FALSE
    )
    if (have_draws) {
      ## population-weight the draws themselves, so the region interval is exact
      p_reg <- matrix(0, n_draw, length(periods))
      for (j in seq_len(nrow(pr))) p_reg <- p_reg + w[j] * draws$p[, pr$k[j], ]
      q <- apply(p_reg, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      rows[[r]] <- cbind(common, data.frame(
        snow_pfpr_median = 100 * q[3, ], snow_pfpr_mean = 100 * colMeans(p_reg),
        snow_pfpr_q2.5 = 100 * q[1, ], snow_pfpr_q25 = 100 * q[2, ],
        snow_pfpr_q75 = 100 * q[4, ], snow_pfpr_q97.5 = 100 * q[5, ],
        snow_pfpr_sd_indep = NA_real_, snow_pfpr_sd_upper = NA_real_,
        interval_is_exact = TRUE))
    } else {
      ## Summary means aggregate exactly by linearity. Quantiles and medians
      ## cannot be aggregated from marginal summaries, so leave them missing.
      idx <- pr$k
      m_mean <- as.numeric(crossprod(w, M$mean[idx, , drop = FALSE]))
      sd_ind <- sqrt(as.numeric(crossprod(w^2, M$sd[idx, , drop = FALSE]^2)))
      sd_up  <- as.numeric(crossprod(w, M$sd[idx, , drop = FALSE]))
      rows[[r]] <- cbind(common, data.frame(
        snow_pfpr_median = NA_real_, snow_pfpr_mean = m_mean,
        snow_pfpr_q2.5 = NA_real_, snow_pfpr_q25 = NA_real_,
        snow_pfpr_q75 = NA_real_, snow_pfpr_q97.5 = NA_real_,
        snow_pfpr_sd_indep = sd_ind, snow_pfpr_sd_upper = sd_up,
        interval_is_exact = FALSE))
    }
  }
  out <- do.call(rbind, rows)
  out <- out[nzchar(out$regkey) & is.finite(out$snow_pfpr_mean), , drop = FALSE]
  attr(out, "signature") <- cache_signature
  saveRDS(out, cache_file)
  out
}

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
region_rows <- list(); status_rows <- vector("list", nrow(registry))
for (i in seq_len(nrow(registry))) {
  survey <- registry[i, , drop = FALSE]
  message("Extracting ", survey$svkey, " (", i, "/", nrow(registry), ")")
  result <- tryCatch(extract_survey_regions(survey),
                     error = function(e) list(failed = TRUE, reason = conditionMessage(e)))
  ok <- is.data.frame(result) && nrow(result) > 0
  if (ok) region_rows[[survey$svkey]] <- result
  status_rows[[i]] <- data.frame(
    svkey = survey$svkey, SurveyId = survey$SurveyId, iso3 = survey$iso3, year = survey$year,
    success = ok,
    regions = if (ok) length(unique(result$regkey)) else 0L,
    regions_low_coverage = if (ok) sum(result$coverage[result$period_index == 1] < MIN_COVERAGE) else NA_integer_,
    reason = if (ok) "" else result$reason %||% "unknown failure",
    stringsAsFactors = FALSE
  )
  if (i %% 20 == 0) message("  ", i, " of ", nrow(registry), " surveys processed")
}
status <- do.call(rbind, status_rows)
write.csv(status, SNOW_STATUS_CSV, row.names = FALSE)
if (!length(region_rows)) stop("No survey-region Snow estimates were produced.")
long <- do.call(rbind, region_rows); rownames(long) <- NULL
long$low_coverage <- long$coverage < MIN_COVERAGE
stopifnot(!anyNA(long$snow_pfpr_mean))
write.csv(long, SNOW_LONG_CSV, row.names = FALSE)

## ---- wide table: one row per survey-region -----------------------------------
long$matched <- long$period_index == vapply(long$year, period_for_year, integer(1))
long$matched[is.na(long$matched)] <- FALSE
matched <- long[long$matched, c("svkey", "regkey", "period", "snow_pfpr_mean",
                                "snow_pfpr_q2.5", "snow_pfpr_q97.5", "snow_pfpr_sd_indep",
                                "snow_pfpr_sd_upper", "interval_is_exact")]
names(matched)[3:6] <- c("survey_period", "snow_pfpr_survey_period",
                         "snow_pfpr_survey_period_q2.5", "snow_pfpr_survey_period_q97.5")
matched$survey_period_note <- ifelse(
  long$year[long$matched] > max(period_end),
  sprintf("survey after %d; latest Snow %s used", max(period_end),
          if (have_draws) "period" else "year"), "")
base <- unique(long[, c("svkey", "SurveyId", "iso3", "year", "region", "regkey",
                        "population_weight", "region_population", "coverage",
                        "low_coverage", "n_snow_polygons", "main_snow_polygon")])
## spread the per-period value across columns: the posterior mean, which is the
## quantity the annual source supplies and is exact for a weighted average
wide_med <- reshape(long[, c("svkey", "regkey", "period", "snow_pfpr_mean")],
                    idvar = c("svkey", "regkey"), timevar = "period", direction = "wide")
names(wide_med) <- sub("^snow_pfpr_mean\\.", "snow_pfpr_", names(wide_med))
wide <- merge(base, wide_med, by = c("svkey", "regkey"))
wide <- merge(wide, matched, by = c("svkey", "regkey"), all.x = TRUE)
wide <- wide[order(wide$iso3, wide$year, wide$regkey), ]
write.csv(wide, SNOW_WIDE_CSV, row.names = FALSE)

message(sprintf(
  "Snow merge complete: %d survey-regions from %d of %d surveys; %d regions (%.0f%%) have < %.0f%% of their population inside the Snow polygons.",
  nrow(wide), sum(status$success), nrow(status), sum(wide$low_coverage),
  100 * mean(wide$low_coverage), 100 * MIN_COVERAGE))
write.csv(data.frame(file = source_files, md5 = unname(source_hashes)),
          snow_out(RESULTS_DIR, "snow_source_provenance", ".csv"), row.names = FALSE)
if (identical(Sys.getenv("SNOW_SKIP_COMPARISONS"), "1")) {
  message("Extraction complete; comparison figures skipped as requested.")
  quit(save = "no", status = 0)
}

## ---- comparison 1: MAP PfPR2-10 for the same survey regions ------------------
if (file.exists(MAP_REGION_CSV)) {
  map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
  cmp <- merge(wide[wide$year <= max(period_end), ],
               map[, c("svkey", "regkey", "pfpr2_10")], by = c("svkey", "regkey"))
  cmp <- cmp[is.finite(cmp$pfpr2_10) & is.finite(cmp$snow_pfpr_survey_period) & !cmp$low_coverage, ]
  agree <- function(z, label) {
    if (nrow(z) < 10) return(NULL)
    data.frame(group = label, n_regions = nrow(z), n_surveys = length(unique(z$svkey)),
               pearson_r = cor(z$snow_pfpr_survey_period, z$pfpr2_10),
               spearman_rho = cor(z$snow_pfpr_survey_period, z$pfpr2_10, method = "spearman"),
               mean_difference_snow_minus_map = mean(z$snow_pfpr_survey_period - z$pfpr2_10),
               mean_absolute_difference = mean(abs(z$snow_pfpr_survey_period - z$pfpr2_10)),
               ci95_covers_map = mean(z$pfpr2_10 >= z$snow_pfpr_survey_period_q2.5 &
                                        z$pfpr2_10 <= z$snow_pfpr_survey_period_q97.5))
  }
  rows <- c(list(agree(cmp, "All survey regions")),
            lapply(sort(unique(cmp$survey_period)), function(p)
              agree(cmp[cmp$survey_period == p, ], paste("Survey period", p))))
  agree_map <- do.call(rbind, rows[!vapply(rows, is.null, TRUE)])
  write.csv(agree_map, snow_out(RESULTS_DIR, "snow_vs_map_agreement", ".csv"), row.names = FALSE)
  cat("\n=== Snow polygon PfPR2-10 versus MAP PfPR2-10 (survey-year matched) ===\n")
  print(within(agree_map, { pearson_r <- round(pearson_r, 3); spearman_rho <- round(spearman_rho, 3)
    mean_difference_snow_minus_map <- round(mean_difference_snow_minus_map, 1)
    mean_absolute_difference <- round(mean_absolute_difference, 1)
    ci95_covers_map <- round(ci95_covers_map, 2) }), row.names = FALSE)
  g <- ggplot(cmp, aes(pfpr2_10, snow_pfpr_survey_period)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(ymin = snow_pfpr_survey_period_q2.5, ymax = snow_pfpr_survey_period_q97.5),
                  width = 0, colour = "grey75", linewidth = 0.3) +
    geom_point(aes(colour = iso3), alpha = 0.7, size = 1.8, show.legend = FALSE) +
    facet_wrap(~ survey_period) +
    labs(x = "MAP PfPR2-10 in the survey year (%)", y = "Snow polygon PfPR2-10, survey period (%)",
         title = "DHS/MIS survey regions: Snow et al. polygon model versus MAP",
         subtitle = if (have_draws) "Population-weighted means and posterior 95% intervals; one colour per country" else
           "Population-weighted posterior means; regional intervals unavailable from summary CSVs") +
    coord_equal(xlim = c(0, 90), ylim = c(0, 90)) + theme_minimal(base_size = 11)
  ggsave(snow_out(RESULTS_DIR, "snow_vs_map_scatter", ".png"), g, width = 10, height = 7, dpi = 200, bg = "white")
}

## ---- comparison 2: DHS-measured parasitaemia (2009 onwards) ------------------
MEASURED_CSV   <- file.path(DATA_DIR, "dhs_prevalence_by_region.csv")
CONVERSION_CSV <- file.path(REPO_ROOT, "results", "rdt_microscopy_conversion.csv")
if (file.exists(MEASURED_CSV) && file.exists(CONVERSION_CSV)) {
  required_packages("malariaAtlas")
  conversion <- read.csv(CONVERSION_CSV, stringsAsFactors = FALSE)
  micro_per_rdt <- conversion$slope[conversion$model == "through_origin"]
  measured <- read.csv(MEASURED_CSV, stringsAsFactors = FALSE)
  measured$micro_equivalent <- ifelse(is.finite(measured$mic), measured$mic, micro_per_rdt * measured$rdt)
  ok <- is.finite(measured$micro_equivalent)
  measured$measured_pfpr <- NA_real_
  measured$measured_pfpr[ok] <- 100 * suppressMessages(as.numeric(malariaAtlas::convertPrevalence(
    pmin(pmax(measured$micro_equivalent[ok], 0), 100) / 100,
    rep(0.5, sum(ok)), rep(5, sum(ok)), rep(2, sum(ok)), rep(10, sum(ok)))))
  ## reconcile the recode's region labels with the boundary's, survey by survey
  measured$svkey <- survey_key(measured$survey)
  meas_rows <- lapply(split(measured, measured$svkey), function(m) {
    b <- wide[wide$svkey == m$svkey[1], ]
    if (!nrow(b)) return(NULL)
    cw <- match_region_keys(m$region, b$region)
    m$regkey_boundary <- cw$to[match(rkey(m$region), cw$from)]
    m[!is.na(m$regkey_boundary), c("svkey", "regkey_boundary", "measured_pfpr", "rdt", "mic")]
  })
  meas <- do.call(rbind, meas_rows)
  names(meas)[2] <- "regkey"
  cmp2 <- merge(wide[wide$year <= max(period_end) & !wide$low_coverage, ], meas, by = c("svkey", "regkey"))
  cmp2 <- cmp2[is.finite(cmp2$measured_pfpr) & is.finite(cmp2$snow_pfpr_survey_period), ]
  cat(sprintf("\nDHS-measured parasitaemia matched to %d survey regions in %d surveys (%d-%d)\n",
              nrow(cmp2), length(unique(cmp2$svkey)), min(cmp2$year), max(cmp2$year)))
  agree2 <- data.frame(
    group = "Surveys with biomarker, up to 2015", n_regions = nrow(cmp2),
    n_surveys = length(unique(cmp2$svkey)),
    pearson_r = cor(cmp2$snow_pfpr_survey_period, cmp2$measured_pfpr),
    spearman_rho = cor(cmp2$snow_pfpr_survey_period, cmp2$measured_pfpr, method = "spearman"),
    mean_difference_snow_minus_measured = mean(cmp2$snow_pfpr_survey_period - cmp2$measured_pfpr),
    mean_absolute_difference = mean(abs(cmp2$snow_pfpr_survey_period - cmp2$measured_pfpr)),
    ci95_covers_measured = mean(cmp2$measured_pfpr >= cmp2$snow_pfpr_survey_period_q2.5 &
                                  cmp2$measured_pfpr <= cmp2$snow_pfpr_survey_period_q97.5))
  ## per-survey summary: which surveys drive the disagreement
  by_survey <- do.call(rbind, lapply(split(cmp2, cmp2$SurveyId), function(z) data.frame(
    SurveyId = z$SurveyId[1], iso3 = z$iso3[1], year = z$year[1], n_regions = nrow(z),
    snow_mean = mean(z$snow_pfpr_survey_period), measured_mean = mean(z$measured_pfpr),
    mean_difference_snow_minus_measured = mean(z$snow_pfpr_survey_period - z$measured_pfpr),
    spearman_rho = if (nrow(z) >= 4) cor(z$snow_pfpr_survey_period, z$measured_pfpr, method = "spearman") else NA_real_)))
  by_survey <- by_survey[order(-abs(by_survey$mean_difference_snow_minus_measured)), ]
  write.csv(by_survey, snow_out(RESULTS_DIR, "snow_vs_measured_by_survey", ".csv"), row.names = FALSE)
  write.csv(agree2, snow_out(RESULTS_DIR, "snow_vs_measured_agreement", ".csv"), row.names = FALSE)
  cat("=== Snow polygon PfPR2-10 versus DHS-measured PfPR2-10 (microscopy-equivalent, age-standardised) ===\n")
  print(within(agree2, { pearson_r <- round(pearson_r, 3); spearman_rho <- round(spearman_rho, 3)
    mean_difference_snow_minus_measured <- round(mean_difference_snow_minus_measured, 1)
    mean_absolute_difference <- round(mean_absolute_difference, 1)
    ci95_covers_measured <- round(ci95_covers_measured, 2) }), row.names = FALSE)
  g2 <- ggplot(cmp2, aes(measured_pfpr, snow_pfpr_survey_period)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(ymin = snow_pfpr_survey_period_q2.5, ymax = snow_pfpr_survey_period_q97.5),
                  width = 0, colour = "grey75", linewidth = 0.3) +
    geom_point(aes(colour = iso3), alpha = 0.7, size = 2, show.legend = FALSE) +
    facet_wrap(~ survey_period) +
    labs(x = "DHS-measured PfPR2-10 (%), microscopy-equivalent, age-standardised",
         y = "Snow polygon PfPR2-10, survey period (%)",
         title = "DHS/MIS survey regions: Snow et al. polygon model versus measured parasitaemia",
         subtitle = if (have_draws) "Surveys with the malaria biomarker module, 2009-2015; posterior 95% intervals" else
           "Surveys with the malaria biomarker module, 2009-2015; posterior means only") +
    coord_equal(xlim = c(0, 90), ylim = c(0, 90)) + theme_minimal(base_size = 11)
  ggsave(snow_out(RESULTS_DIR, "snow_vs_measured_scatter", ".png"), g2, width = 9, height = 5.5, dpi = 200, bg = "white")
}
message("Outputs: ", SNOW_LONG_CSV, ", ", SNOW_WIDE_CSV, ", ", SNOW_STATUS_CSV)
