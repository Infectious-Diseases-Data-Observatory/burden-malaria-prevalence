# =============================================================================
# 02_extract_map_pfpr.R
# Extract survey-year-matched, population-weighted MAP PfPR2-10 for every region
# in eligible DHS/MIS surveys conducted from 2000 through 2025.
#
# Inputs:
#   data/derived_dhs/survey_registry.csv
#   data/dhs_boundaries/<SurveyId>.rds
#   cached MAP rasters and GPW population density
#
# Outputs:
#   data/derived_dhs/map_pfpr_by_survey_region.csv
#   results/dhs_rebuild/map_extraction_status.csv
#
# Missing 2025 MAP availability is recorded as a status rather than silently
# substituting 2024.
# =============================================================================

source("R_dhs/00_config.R")
required_packages(c("malariaAtlas", "terra", "sf"))

if (!file.exists(SURVEY_REGISTRY_CSV)) {
  stop("Survey registry missing. Run R_dhs/01_access_dhs_data.R first.")
}
if (!file.exists(GPW_TIF)) stop("Population raster missing: ", GPW_TIF)

registry <- read.csv(SURVEY_REGISTRY_CSV, stringsAsFactors = FALSE)
registry <- registry[
  registry$year >= DHS_START_YEAR & registry$year <= DHS_END_YEAR,
  ,
  drop = FALSE
]
population <- terra::rast(GPW_TIF)

get_map_raster <- function(year) {
  output <- file.path(MAP_RASTER_DIR, sprintf("pfpr2_10_%d.tif", year))
  if (file.exists(output)) return(terra::rast(output))

  message("MAP raster not cached for ", year, "; requesting it from malariaAtlas.")
  download_dir <- file.path(tempdir(), paste0("map_", year))
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(tryCatch(
    malariaAtlas::getRaster(
      dataset_id = MAP_DATASET_ID,
      year = year,
      extent = MAP_AFRICA_EXTENT,
      file_path = download_dir
    ),
    error = function(e) NULL
  ))
  tif <- list.files(
    download_dir,
    pattern = "\\.tiff?$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (!length(tif)) return(NULL)

  raster <- terra::rast(tif)
  mean_layer <- grep("_1$", names(raster))
  raster <- raster[[if (length(mean_layer)) mean_layer[1] else 1]]
  terra::writeRaster(raster, output, overwrite = TRUE)
  terra::rast(output)
}

extract_survey_regions <- function(survey) {
  cache_file <- file.path(MAP_REGION_DIR, paste0(survey$svkey, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  if (!file.exists(survey$boundary_file)) {
    stop("boundary unavailable")
  }

  boundary <- readRDS(survey$boundary_file)
  if (!"DHSREGEN" %in% names(boundary)) {
    stop("boundary lacks DHSREGEN")
  }
  pfpr <- get_map_raster(survey$year)
  if (is.null(pfpr)) stop("MAP raster unavailable for survey year")

  boundary <- sf::st_make_valid(boundary)
  polygons <- terra::makeValid(terra::vect(boundary))
  weights <- terra::resample(
    terra::crop(population, pfpr),
    pfpr,
    method = "bilinear"
  )
  weights <- terra::mask(weights, pfpr)

  numerator <- terra::extract(
    pfpr * weights,
    polygons,
    fun = sum,
    na.rm = TRUE,
    exact = TRUE,
    ID = FALSE
  )[[1]]
  denominator <- terra::extract(
    weights,
    polygons,
    fun = sum,
    na.rm = TRUE,
    exact = TRUE,
    ID = FALSE
  )[[1]]

  out <- data.frame(
    svkey = survey$svkey,
    SurveyId = survey$SurveyId,
    iso3 = survey$iso3,
    year = survey$year,
    region = as.character(boundary$DHSREGEN),
    regkey = rkey(boundary$DHSREGEN),
    pfpr2_10 = 100 * numerator / denominator,
    population_weight = denominator
  )
  out <- out[
    nzchar(out$regkey) & is.finite(out$pfpr2_10) &
      is.finite(out$population_weight),
    ,
    drop = FALSE
  ]
  saveRDS(out, cache_file)
  out
}

region_rows <- list()
status_rows <- vector("list", nrow(registry))
for (i in seq_len(nrow(registry))) {
  survey <- registry[i, , drop = FALSE]
  result <- tryCatch(
    extract_survey_regions(survey),
    error = function(e) structure(NULL, reason = conditionMessage(e))
  )
  ok <- !is.null(result) && nrow(result) > 0
  if (ok) region_rows[[survey$svkey]] <- result
  status_rows[[i]] <- data.frame(
    svkey = survey$svkey,
    SurveyId = survey$SurveyId,
    iso3 = survey$iso3,
    year = survey$year,
    success = ok,
    regions = if (ok) nrow(result) else 0L,
    reason = if (ok) "" else attr(result, "reason") %||% "unknown failure"
  )
}

status <- do.call(rbind, status_rows)
write.csv(status, MAP_STATUS_CSV, row.names = FALSE)

if (!length(region_rows)) stop("No survey-region PfPR estimates were produced.")
region_pfpr <- do.call(rbind, region_rows)
rownames(region_pfpr) <- NULL
write.csv(region_pfpr, MAP_REGION_CSV, row.names = FALSE)

message(
  "MAP extraction complete: ", nrow(region_pfpr), " survey-regions from ",
  sum(status$success), " of ", nrow(status), " surveys."
)
if (any(!status$success & status$year == 2025)) {
  message(
    "One or more 2025 surveys lack a 2025 MAP surface; see ",
    MAP_STATUS_CSV, ". No 2024 substitution was made."
  )
}
