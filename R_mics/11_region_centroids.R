#!/usr/bin/env Rscript
# Boundary centroids of the MICS analysis regions, for the Sahel subgroup refit.
# Same method as the cached DHS table data/derived_dhs/dhs_region_centroids.csv:
# spherical (s2) centroid of the made-valid polygon, one row per survey region.
# Writes data/derived_mics/mics_region_centroids.csv with the DHS table's columns.
suppressPackageStartupMessages({ library(sf); library(data.table) })
sf_use_s2(TRUE)
source("R_cbh/R/geography.R")   # cbh_rkey()
out <- "data/derived_mics"
registry <- fread(file.path(out, "survey_registry_mics.csv"))
pfpr <- fread(file.path(out, "map_pfpr_by_survey_region_mics.csv"))
stopifnot(!anyDuplicated(registry$SurveyId), all(file.exists(registry$boundary_file)))
rows <- lapply(seq_len(nrow(registry)), function(i) {
  r <- registry[i]; b <- st_make_valid(readRDS(r$boundary_file))
  xy <- st_coordinates(suppressWarnings(st_centroid(st_geometry(b))))
  data.table(SurveyId = r$SurveyId, iso3 = r$iso3, year = r$year, regkey = cbh_rkey(b$DHSREGEN),
             lon = xy[, "X"], lat = xy[, "Y"])
})
centroids <- rbindlist(rows)
stopifnot(!anyDuplicated(centroids[, .(SurveyId, regkey)]), all(is.finite(centroids$lon)),
          all(is.finite(centroids$lat)), all(nzchar(centroids$regkey)),
          all(centroids$lat > -35 & centroids$lat < 25 & centroids$lon > -20 & centroids$lon < 52))
# Every region with MAP prevalence must have a centroid.
miss <- fsetdiff(unique(pfpr[, .(SurveyId, regkey)]), centroids[, .(SurveyId, regkey)])
if (nrow(miss)) stop("Regions without a centroid: ", paste(head(paste(miss$SurveyId, miss$regkey)), collapse = ", "))
fwrite(centroids, file.path(out, "mics_region_centroids.csv"))
cat("MICS region centroids:", nrow(centroids), "regions in", uniqueN(centroids$SurveyId), "surveys\n")
