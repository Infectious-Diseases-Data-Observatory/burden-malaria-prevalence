#!/usr/bin/env Rscript
# Verify fitted exposures and nuisance inputs against the exact prepared data.
# Only aggregate pass/fail records leave the ignored private model directory.
source("R_cbh/load_pipeline.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--retain-later-surveys"))
id <- if("--retain-later-surveys" %in% args) "age_band_snow_2000_2015_retrospective_v1" else "age_band_snow_2000_2015_v1"
private <- file.path("data/derived_cbh/models",id); out <- file.path("results/cbh",id)
prepared <- readRDS(file.path(private,"dataset.rds"))
signature <- prepared$signature
d <- prepared$data; rm(prepared); gc(FALSE)
jobs <- cbh_read_csv(file.path(out,"selected_fit_manifest.csv"))
results <- list()
for(i in seq_len(nrow(jobs))) {
  job <- jobs[i,]
  use <- d$age_band==job$age_band & if(job$series=="snow") TRUE else is.finite(d$map_pfpr_pct)
  model_path <- file.path(private,job$model_file)
  stopifnot(identical(cbh_file_hash(model_path),job$model_file_md5))
  saved <- readRDS(model_path)
  f <- saved$fit
  stopifnot(identical(saved$input_signature,signature),nrow(f$model)==sum(use))
  expected <- if(job$series=="map_matched") d$map_pfpr_pct[use] else d$snow_pfpr_pct[use]
  stopifnot(identical(f$model$pfpr_pct,expected),
    max(abs(exp(f$model[["offset(log(band_years))"]])-d$band_years[use]))<1e-12)
  nuisance <- setdiff(names(f$model),c("pfpr_pct","offset(log(band_years))"))
  for(v in nuisance) stopifnot(isTRUE(all.equal(as.vector(f$model[[v]]),as.vector(d[[v]][use]),check.attributes=FALSE)))
  results[[i]] <- data.frame(fit_id=job$fit_id,rows=sum(use),exposure_matches_prepared=TRUE,
    outcome_and_all_nuisance_variables_match=TRUE,offset_matches=TRUE,
    prepared_signature=signature,fit_file_md5=cbh_file_hash(model_path))
  rm(saved,f);gc(FALSE)
}
cbh_atomic_csv(do.call(rbind,results),file.path(out,"saved_fit_input_checks.csv"))
message("All ",nrow(jobs)," saved fits match their intended prepared inputs exactly")
