#!/usr/bin/env Rscript
# Independent checks using only geographic aggregates and exposure draws.
source("R_cbh/load_pipeline.R")
library(sf)
library(terra)
out <- "results/cbh/age_band_snow_2000_2015_v1"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
provenance <- cbh_read_csv("results/dhs_rebuild/snow_source_provenance_annual_means.csv")
stopifnot(all(vapply(provenance$file,cbh_file_hash,"")==provenance$md5))
ann <- cbh_read_csv("data/snow_prevalence_model/snow_pfpr_by_polygon_year.csv")
poly <- st_read("data/snow_prevalence_model/polygons_520.gpkg",quiet=TRUE)
draws <- readRDS("data/snow_prevalence_model/p_draws_annual.rds")
cbh_unique(ann,c("Admin_ID","year"),"Annual polygon estimates")
stopifnot(nrow(ann)==520*16,identical(sort(unique(ann$year)),2000:2015),
  setequal(poly$Admin_ID,ann$Admin_ID),identical(draws$scale,"proportion"))
means <- 100*apply(draws$p,c(2,3),mean)
ix <- cbind(match(ann$Admin_ID,draws$Admin_ID),match(ann$year,draws$years))
csv_draw_error <- max(abs(ann$PfPR_mean-means[ix]))
stopifnot(csv_draw_error<1e-8)
panel <- cbh_read_csv("data/derived_dhs/snow_pfpr_by_survey_region_long_annual_means.csv")
cbh_unique(panel,c("svkey","regkey","period"),"Snow panel")
stopifnot(all(panel$period %in% 2000:2015),all(is.na(panel$snow_pfpr_q2.5)),
  all(is.na(panel$snow_pfpr_median)),all(!panel$interval_is_exact))
checks <- list()
for (sv in unique(panel$svkey)) {
  z <- panel[panel$svkey==sv,]
  w <- readRDS(file.path("data/snow_region_cache/annual_means",paste0(sv,"_weights.rds")))
  cached <- readRDS(file.path("data/snow_region_cache/annual_means",paste0(sv,".rds")))
  stopifnot(identical(w$signature,attr(cached,"signature")))
  p <- w$pieces
  stopifnot(all(poly$Country_ID[match(p$Admin_ID,poly$Admin_ID)]==unique(z$iso3)))
  for (r in unique(p$row_id)) {
    q <- p[p$row_id==r & p$pop>0,]; if(!nrow(q)) next
    weights <- q$pop/sum(q$pop)
    # Aggregate the draws, independently of the CSV-matrix path in script 51.
    region_draws <- matrix(0,nrow(draws$p),length(draws$years))
    for(j in seq_len(nrow(q))) region_draws <- region_draws +
      weights[j]*draws$p[,match(q$Admin_ID[j],draws$Admin_ID),]
    actual <- z[z$region==w$region[r],]
    stopifnot(nrow(actual)==16)
    difference <- actual$snow_pfpr_mean-100*colMeans(region_draws)[match(as.integer(actual$period),draws$years)]
    stopifnot(max(abs(difference))<1e-8,abs(sum(weights)-1)<1e-12)
    checks[[length(checks)+1L]] <- data.frame(survey=sv,region=w$region[r],
      polygons=nrow(q),coverage=unique(actual$coverage),maximum_mean_error=max(abs(difference)))
  }
}
checks <- do.call(rbind,checks)
cbh_atomic_csv(checks,file.path(out,"aggregation_checks.csv"))
# Independently verify fractional-cell weighting against an analytic example.
r <- rast(nrows=1,ncols=2,xmin=0,xmax=2,ymin=0,ymax=1,crs="EPSG:3857")
values(r) <- c(100,300)
g <- st_sfc(st_polygon(list(rbind(c(.5,0),c(1.25,0),c(1.25,1),c(.5,1),c(.5,0)))),crs=3857)
ee <- exactextractr::exact_extract(r,g,"sum",progress=FALSE)
te <- terra::extract(r,vect(g),fun=sum,exact=TRUE,ID=FALSE)[[1]]
stopifnot(abs(ee-125)<1e-6,abs(te-125)<1e-6,abs(ee-te)<1e-6)
# Also compare extractors on actual longitude/latitude geometry. Tiny numerical
# differences are expected because terra uses geodesic cell fractions.
sf_use_s2(FALSE)
registry <- cbh_read_csv("data/derived_dhs/survey_registry.csv")
survey <- registry[registry$svkey=="GM61FL",]
boundary <- st_make_valid(st_transform(readRDS(survey$boundary_file),4326))
pieces <- suppressWarnings(st_intersection(boundary[1,"DHSREGEN"],
  st_make_valid(poly[poly$Country_ID=="GMB","Admin_ID"])))
pop_path <- "data/pop/gpw_v4_population_density_rev11_2020_2.5m.tif"
pop <- crop(rast(pop_path),vect(boundary)); pop <- pop*cellSize(pop,unit="km")
ee <- exactextractr::exact_extract(pop,pieces,"sum",progress=FALSE)
te <- terra::extract(pop,vect(pieces),fun=sum,exact=TRUE,ID=FALSE)[[1]]
extractor <- data.frame(polygon=pieces$Admin_ID,exactextractr=ee,terra=te,
  relative_difference=(ee-te)/pmax(te,1))
stopifnot(max(abs(extractor$relative_difference))<.001)
cbh_atomic_csv(extractor,file.path(out,"extractor_comparison.csv"))
wide <- cbh_read_csv("data/derived_dhs/snow_pfpr_by_survey_region_annual_means.csv")
stopifnot(all(is.na(wide$snow_pfpr_survey_period[wide$year>2015])),
  all(panel$coverage>=0),all(panel$coverage<=1.01))
old_path <- "data/derived_dhs/snow_pfpr_by_survey_region_long_annual.csv"
if(file.exists(old_path)) {
  old <- cbh_read_csv(old_path)
  j <- match(cbh_key(panel,c("svkey","regkey","period")),cbh_key(old,c("svkey","regkey","period")))
  cmp <- panel[c("svkey","iso3","regkey","period","coverage","snow_pfpr_mean")]
  cmp$previous_mean <- old$snow_pfpr_mean[j]
  cmp$difference <- cmp$snow_pfpr_mean-cmp$previous_mean
  cbh_atomic_csv(cmp,file.path(out,"previous_extraction_comparison.csv"))
}
writeLines(c("Annual polygon CSV / posterior draw maximum mean error (percentage points):",
  format(csv_draw_error,digits=12),paste("Region weight checks:",nrow(checks)),
  paste("Maximum region mean error:",max(checks$maximum_mean_error)),
  "Fractional raster-cell population sum: 125, verified analytically and against terra.",
  "Country identifiers, complete annual keys, source/cache signatures and no post-2015 carry-forward: passed.",
  "Population coverage is retained; rows below 50% are excluded in dataset preparation."),
  file.path(out,"extraction_validation.txt"))
message("Extraction audit passed: ",nrow(checks)," survey-regions")
inputs <- c("R_cbh/snow/00_audit_extraction.R",provenance$file,
  "data/snow_prevalence_model/p_draws_annual.rds",
  "data/derived_dhs/snow_pfpr_by_survey_region_long_annual_means.csv",
  survey$boundary_file,if(file.exists(old_path)) old_path)
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),
  file.path(out,"extraction_audit_provenance.csv"))
