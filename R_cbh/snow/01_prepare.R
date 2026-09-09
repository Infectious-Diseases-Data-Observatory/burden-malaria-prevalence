#!/usr/bin/env Rscript
# Reuse validated child-band shards, replacing only the exposure and time range.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% "--retain-later-surveys"))
retain_later <- "--retain-later-surveys" %in% args
id <- if (retain_later) "age_band_snow_2000_2015_retrospective_v1" else "age_band_snow_2000_2015_v1"
out <- file.path("results/cbh", id)
private <- file.path("data/derived_cbh/models", id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(private, recursive = TRUE, showWarnings = FALSE)
spec <- cbh_trial_spec()
snow_path <- "data/derived_dhs/snow_pfpr_by_survey_region_long_annual_means.csv"
snow <- cbh_read_csv(snow_path)
snow$exposure_year <- as.integer(snow$period)
cbh_unique(snow, c("svkey", "regkey", "exposure_year"), "Snow region-year panel")
stopifnot(all(snow$exposure_year %in% 2000:2015), all(is.finite(snow$snow_pfpr_mean)),
  all(snow$snow_pfpr_mean >= 0 & snow$snow_pfpr_mean <= 100))
meta <- readRDS("data/derived_cbh/manifest.rds")
stopifnot(isTRUE(meta$complete), !any(meta$manifest$status == "failed"))
m <- meta$manifest
cbh_atomic_csv(m[!m$status %in% c("built", "cached"), c("survey", "country", "survey_year", "status")],
  file.path(out, "unavailable_shards.csv"))
m <- m[m$status %in% c("built", "cached"), ]
incidence <- cbh_read_csv(spec$incidence_panel)
required <- c("death", "age_band", "calendar_year", "band_years", "survey", "country", "region", spec$covariates)
retained <- unique(c(required, "pfpr_pct", "child_id", "regkey", "entry_year", "band_entry_cmc",
  "band_end_cmc", "interview_cmc", "survey_year", "survey_weight", "hiv_incidence_status"))
rows <- selection <- region_audit <- list()
cutoff_end <- (2016L - 1900L) * 12L + 1L # Jan 2016; band-end boundary is exclusive.
for (i in seq_len(nrow(m))) {
  obj <- readRDS(file.path("data/derived_cbh", m$file[i]))
  stopifnot(identical(obj$signature, m$signature[i]))
  d <- cbh_attach_incidence(obj$data, incidence)
  d$map_pfpr_pct <- d$pfpr_pct
  period_ok <- d$entry_year >= 2000 & d$entry_year <= 2015 & d$band_end_cmc <= cutoff_end
  if (!retain_later) period_ok <- period_ok & d$survey_year <= 2015 & d$interview_cmc < cutoff_end
  j <- match(cbh_key(d, c("survey", "regkey", "entry_year")),
    cbh_key(snow, c("svkey", "regkey", "exposure_year")))
  d$snow_coverage <- snow$coverage[j]
  d$snow_pfpr_pct <- snow$snow_pfpr_mean[j]
  has_snow <- !is.na(j) & is.finite(d$snow_pfpr_pct)
  stopifnot(all(d$country[has_snow] == snow$iso3[j[has_snow]]),
    all(d$entry_year[has_snow] == snow$exposure_year[j[has_snow]]))
  # Carry forward script 51's explicit 50% population-coverage threshold.
  good_coverage <- has_snow & is.finite(d$snow_coverage) & d$snow_coverage >= .5 & d$snow_coverage <= 1.01
  complete <- complete.cases(d[required])
  for (v in required) if (is.numeric(d[[v]])) complete <- complete & is.finite(d[[v]])
  keep <- period_ok & good_coverage & complete
  map_available <- is.finite(d$map_pfpr_pct) & d$map_pfpr_pct >= 0 & d$map_pfpr_pct <= 100
  selection[[i]] <- data.frame(survey=m$survey[i], country=m$country[i], survey_year=m$survey_year[i],
    eligible_rows=nrow(d), within_2000_2015=sum(period_ok),
    missing_snow_rows=sum(period_ok & !has_snow),
    inadequate_snow_coverage_rows=sum(period_ok & has_snow & !good_coverage),
    incomplete_confounder_rows=sum(period_ok & good_coverage & !complete),
    snow_rows=sum(keep), snow_deaths=sum(d$death[keep]),
    matched_map_rows=sum(keep & map_available), matched_map_deaths=sum(d$death[keep & map_available]))
  audit <- as.data.table(d[period_ok, c("survey", "country", "regkey", "region", "entry_year",
    "snow_coverage", "snow_pfpr_pct", "map_pfpr_pct")])
  region_audit[[i]] <- audit[, .(rows=.N, snow_coverage=first(snow_coverage),
    snow_pfpr_pct=first(snow_pfpr_pct), map_pfpr_pct=first(map_pfpr_pct)),
    by=.(survey,country,regkey,region,entry_year)]
  d$pfpr_pct <- d$snow_pfpr_pct
  rows[[i]] <- d[keep, c(retained, "map_pfpr_pct", "snow_pfpr_pct", "snow_coverage"), drop=FALSE]
  rm(obj,d); if (i %% 20L == 0L) message("Prepared ",i,"/",nrow(m)," surveys")
}
d <- as.data.frame(rbindlist(rows)); rm(rows); gc(FALSE)
cbh_unique(d, c("child_id", "age_band"), "Snow child-band data")
stopifnot(nrow(d)>0, all(d$entry_year %in% 2000:2015), all(d$band_end_cmc<=cutoff_end),
  all(d$band_end_cmc<=d$interview_cmc), all(d$band_entry_cmc>=d$interview_cmc-60),
  all(d$pfpr_pct==d$snow_pfpr_pct))
if (!retain_later) stopifnot(all(d$survey_year<=2015), all(cbh_cmc_year(d$interview_cmc)<=2015))
prepared <- cbh_trial_prepare(d,meta$config$age_bands$age_band,spec)
# Preserve the established scaling; both exposures receive identical covariates.
original <- cbh_read_csv(file.path("results/cbh",spec$id,"scaling.csv"))
for (i in seq_len(nrow(original))) {
  v <- original$variable[i]
  prepared$data[[paste0("z_",v)]] <- (prepared$data[[v]]-original$mean[i])/original$sd[i]
}
prepared$scaling <- original
selection <- as.data.frame(rbindlist(selection))
audit <- as.data.frame(rbindlist(region_audit))
inputs <- c(snow_path,spec$incidence_panel,"data/derived_cbh/manifest.rds",
  file.path("results/cbh",spec$id,"scaling.csv"),"R_cbh/snow/01_prepare.R","R_cbh/analysis/model.R")
hashes <- vapply(inputs,cbh_file_hash,"")
prepared$signature <- cbh_hash(list(hashes,m$signature,retain_later))
prepared$spec <- spec; prepared$retain_later_surveys <- retain_later
prepared$selection <- selection
prepared$input_manifest <- meta
cbh_atomic_rds(prepared,file.path(private,"dataset.rds"))
cbh_atomic_csv(selection,file.path(out,"selection.csv"))
cbh_atomic_csv(audit,file.path(out,"region_year_join_audit.csv"))
cbh_atomic_csv(original,file.path(out,"scaling.csv"))
cbh_atomic_csv(data.frame(file=inputs,md5=unname(hashes)),file.path(out,"preparation_provenance.csv"))
counts <- as.data.table(prepared$data)[,.(rows=.N,deaths=sum(death),countries=uniqueN(country),surveys=uniqueN(survey)),by=age_band]
cbh_atomic_csv(counts,file.path(out,"sample_by_age_band.csv"))
cbh_atomic_csv(as.data.table(prepared$data)[,.(rows=.N,deaths=sum(death),surveys=uniqueN(survey)),by=country],file.path(out,"sample_by_country.csv"))
message("Snow dataset: ",nrow(d)," records, ",sum(d$death)," deaths. Without MAP: ",sum(!is.finite(d$map_pfpr_pct)))
