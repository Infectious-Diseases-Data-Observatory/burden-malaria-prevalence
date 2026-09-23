#!/usr/bin/env Rscript
# Per-survey analysis-region polygons, annual regional MAP PfPR2-10, and the MICS
# registry, rules and boundary tables the child age-band builder needs.
# PfPR extraction replicates R_dhs/41_extract_map_window_years.R exactly: exact fractional
# polygon coverage, GPW v4 2020 density weights resampled to the MAP grid (no cell area),
# masked to MAP cells, for registry year - 5 to year + 2.
suppressPackageStartupMessages({ library(sf); library(terra); library(data.table) })
sf_use_s2(FALSE)
source("R_cbh/R/geography.R")   # cbh_rkey(), identical to R_dhs rkey()
root <- normalizePath(".")
source("R_mics/config/paths.R")   # mics_admin1_shp
snow <- mics_admin1_shp
out <- "data/derived_mics"; bdir <- file.path(out, "boundaries"); dir.create(bdir, showWarnings = FALSE)
setup <- fread(file.path(out, "survey_setup.csv"))
rmap <- fread(file.path(out, "region_map.csv"))
inv <- fread("results/mics_inventory/survey_inventory.csv")
setup[, year := inv$interview_year[match(folder, inv$survey)]]
setup[, iso3 := substr(folder, 1, 3)]
stopifnot(all(is.finite(setup$year)))

source_polygons <- function(b) {
  if (b %in% c("CAF_analysis", "GNB_analysis")) {
    x <- readRDS(sprintf("%s/%s_analysis_regions.rds", bdir, sub("_analysis", "", b))); x$name <- x$analysis_region; return(x) }
  if (startsWith(b, "SNOW:")) { a <- st_read(snow, quiet = TRUE); a <- a[a$Country_ID == sub("SNOW:", "", b), ]; a$name <- a$AFR_Admin; return(a) }
  x <- readRDS(sprintf("data/dhs_boundaries/%s.rds", b)); x$name <- x$DHSREGEN; x
}
# 1. Analysis-region polygons per survey: union of the listed source polygons.
for (i in seq_len(nrow(setup))) {
  s <- setup[i]; src <- st_make_valid(source_polygons(s$boundary))
  m <- unique(rmap[svkey == s$svkey, .(analysis_region, polygons)])
  stopifnot(!anyDuplicated(m$analysis_region))
  polys <- lapply(seq_len(nrow(m)), function(j) {
    want <- strsplit(m$polygons[j], ";", fixed = TRUE)[[1]]
    hit <- src[src$name %in% want, ]
    if (nrow(hit) != length(want)) stop(s$svkey, ": polygon(s) not found for ", m$analysis_region[j])
    st_sf(DHSREGEN = m$analysis_region[j], geometry = st_union(st_geometry(hit)))
  })
  b <- st_make_valid(do.call(rbind, polys)); st_crs(b) <- 4326
  stopifnot(!anyDuplicated(cbh_rkey(b$DHSREGEN)))
  saveRDS(b, file.path(bdir, paste0(s$svkey, ".rds")))
}
# 2. Annual PfPR per analysis region.
rf <- list.files("data/map_annual", pattern = "pfpr2_10_[0-9]{4}[.]tif$", full.names = TRUE)
ry <- as.integer(regmatches(basename(rf), regexpr("[0-9]{4}", basename(rf))))
rasters <- setNames(lapply(rf, rast), as.character(ry))
w_full <- resample(crop(rast("data/pop/gpw_v4_population_density_rev11_2020_2.5m.tif"), rasters[[1]]),
                   rasters[[1]], method = "bilinear")
annual <- list(); byreg <- list()
for (i in seq_len(nrow(setup))) {
  s <- setup[i]; b <- readRDS(file.path(bdir, paste0(s$svkey, ".rds")))
  v <- makeValid(vect(b)); years <- (s$year - 5L):(s$year + 2L); years <- years[as.character(years) %in% names(rasters)]
  e <- ext(v) + 0.5; stk <- crop(rast(rasters[as.character(years)]), e)
  w <- mask(crop(w_full, e), stk[[nlyr(stk)]])
  den <- extract(w, v, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)[[1]]
  num <- extract(stk * w, v, fun = sum, na.rm = TRUE, exact = TRUE, ID = FALSE)
  rk <- cbh_rkey(b$DHSREGEN)
  a <- rbindlist(lapply(seq_along(years), function(j) data.table(svkey = s$svkey, iso3 = s$iso3,
    survey_year = s$year, year = years[j], regkey = rk, pfpr2_10 = 100 * num[[j]] / den)))
  a <- a[nzchar(regkey) & is.finite(pfpr2_10)]
  annual[[s$svkey]] <- a
  sy <- a[year == s$year]
  byreg[[s$svkey]] <- data.table(svkey = s$svkey, SurveyId = s$svkey, iso3 = s$iso3, year = s$year,
    region = b$DHSREGEN, regkey = rk, pfpr2_10 = sy$pfpr2_10[match(rk, sy$regkey)], population_weight = den)[regkey %in% a$regkey]
  message(sprintf("%-14s %2d regions, %d years, PfPR %s", s$svkey, nrow(b), length(years),
    if (nrow(a)) sprintf("%.1f-%.1f", min(a$pfpr2_10), max(a$pfpr2_10)) else "none (no MAP surface)"))
}
fwrite(rbindlist(annual), file.path(out, "map_pfpr_window_years_mics.csv"))
fwrite(rbindlist(byreg), file.path(out, "map_pfpr_by_survey_region_mics.csv"))
# 3. Registry and builder rules. region_var is v024 in the DHS-shaped recode the adapter writes.
reg <- setup[, .(iso3, CountryName = iso3, year, SurveyType = "MICS", SurveyId = svkey,
  FileName = folder, no_extension = svkey, svkey,
  local_recode = file.path(root, out, "recodes", paste0(svkey, ".rds")),
  boundary_file = file.path(root, bdir, paste0(svkey, ".rds")),
  recode_available = TRUE, boundary_available = TRUE)]
fwrite(reg, file.path(out, "survey_registry_mics.csv"))
fwrite(setup[, .(svkey, region_var = "v024", cmc_offset_months = 0L, calendar = "gregorian",
  group_donor = "", group_fine = "", group_coarse = "", strata_var = "v022")], file.path(out, "survey_rules_mics.csv"))
cat("surveys:", nrow(setup), "| with PfPR:", uniqueN(rbindlist(annual)$svkey), "\n")
