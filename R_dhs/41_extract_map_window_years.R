# =============================================================================
# 41_extract_map_window_years.R — annual MAP PfPR2-10 per survey region for the
# six calendar years a survey's five 12-month windows can touch.
#
# Script 02 extracts the survey year only and script 17 lags 0-4 for the main
# sample. The person-time design (script 40) needs every region of every survey,
# including those below 1% prevalence, for survey_year - 5 to survey_year + 1
# (fieldwork can run into the following calendar year), so every window can be
# paired with the years it actually spans. Extraction follows script 02 exactly: exact fractional
# coverage of each admin-1 polygon, weighted by the time-invariant GPW
# population surface. Years before 2000 have no MAP surface; the pairing step
# (script 42) carries the 2000 value back to 1999 for windows that straddle it
# and drops windows centred before 2000.
#
# Output
#   data/derived_dhs/map_pfpr_window_years.csv   svkey, regkey, year, pfpr2_10
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("terra", "sf"))

YEARS_BACK <- 5L
# Fieldwork often runs into the calendar year after the registry's survey year,
# so the first window can touch survey_year + 1.
YEARS_FORWARD <- 2L
CACHE <- file.path(DATA_DIR, "map_window_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
OUTPUT <- file.path(DERIVED_DIR, "map_pfpr_window_years.csv")

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
map <- read.csv(MAP_REGION_CSV, stringsAsFactors = FALSE)
registry <- registry[registry$svkey %in% map$svkey, , drop = FALSE]

raster_files <- list.files(MAP_RASTER_DIR, pattern = "pfpr2_10_[0-9]{4}[.]tif$",
                           full.names = TRUE)
raster_years <- as.integer(regmatches(basename(raster_files),
                                      regexpr("[0-9]{4}", basename(raster_files))))
rasters <- setNames(lapply(raster_files, terra::rast), as.character(raster_years))
reference <- rasters[[1]]
weights_full <- terra::resample(
  terra::crop(terra::rast(GPW_TIF), reference), reference, method = "bilinear")
message("MAP surfaces ", min(raster_years), "-", max(raster_years), "; ",
        nrow(registry), " surveys")

extract_survey <- function(survey) {
  cache <- file.path(CACHE, sprintf("%s_b%df%d.rds", survey$svkey, YEARS_BACK, YEARS_FORWARD))
  if (file.exists(cache)) return(readRDS(cache))
  boundary_file <- file.path(BOUNDARY_DIR, paste0(survey$SurveyId, ".rds"))
  if (!file.exists(boundary_file)) stop("boundary unavailable")
  boundary <- sf::st_make_valid(readRDS(boundary_file))
  if (!"DHSREGEN" %in% names(boundary)) stop("boundary lacks DHSREGEN")
  polygons <- terra::makeValid(terra::vect(boundary))
  years <- (survey$year - YEARS_BACK):(survey$year + YEARS_FORWARD)
  years <- years[as.character(years) %in% names(rasters)]
  if (!length(years)) stop("no MAP surface in the window years")
  # Crop to the survey's own extent first: masking and extracting over the
  # whole-Africa grid made each survey take over a minute.
  extent <- terra::ext(polygons) + 0.5
  stack <- terra::crop(terra::rast(rasters[as.character(years)]), extent)
  w <- terra::mask(terra::crop(weights_full, extent), stack[[terra::nlyr(stack)]])
  denominator <- terra::extract(w, polygons, fun = sum, na.rm = TRUE,
                                exact = TRUE, ID = FALSE)[[1]]
  numerator <- terra::extract(stack * w, polygons, fun = sum, na.rm = TRUE,
                              exact = TRUE, ID = FALSE)
  regkeys <- rkey(boundary$DHSREGEN)
  out <- do.call(rbind, lapply(seq_along(years), function(j) data.frame(
    svkey = survey$svkey, iso3 = survey$iso3, survey_year = survey$year,
    year = years[j], regkey = regkeys,
    pfpr2_10 = 100 * numerator[[j]] / denominator, stringsAsFactors = FALSE)))
  out <- out[nzchar(out$regkey) & is.finite(out$pfpr2_10), , drop = FALSE]
  saveRDS(out, cache)
  out
}

rows <- list()
for (i in seq_len(nrow(registry))) {
  survey <- registry[i, , drop = FALSE]
  result <- tryCatch(extract_survey(survey), error = function(e) {
    message("  skip ", survey$svkey, ": ", conditionMessage(e)); NULL })
  if (!is.null(result) && nrow(result)) {
    rows[[survey$svkey]] <- result
    message(sprintf("  %s: %d regions x %d years", survey$svkey,
                    length(unique(result$regkey)), length(unique(result$year))))
  }
}
long <- do.call(rbind, rows)
rownames(long) <- NULL

# The survey-year value must reproduce script 02's extraction for the same
# region; anything else means the polygons or weights differ.
check <- merge(long[long$year == long$survey_year, c("svkey", "regkey", "pfpr2_10")],
               map[, c("svkey", "regkey", "pfpr2_10")],
               by = c("svkey", "regkey"), suffixes = c("_window", "_panel"))
message(sprintf("Survey-year check against script 02: %d regions, max |difference| %.4f points",
                nrow(check), max(abs(check$pfpr2_10_window - check$pfpr2_10_panel))))
write.csv(long, OUTPUT, row.names = FALSE)
message("Wrote ", basename(OUTPUT), ": ", nrow(long), " region-years across ",
        length(unique(long$svkey)), " surveys")
