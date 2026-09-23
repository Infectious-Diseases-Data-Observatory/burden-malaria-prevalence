#!/usr/bin/env Rscript
# Assemble the audited regional overlay; no source extraction or HIV refit.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/primary/specification.R")
source("R_cbh/covariates/settings.R")
source("R_cbh/primary/settings.R")
library(data.table)
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"));cs <- cbh_covariate_settings();cfg <- cbh_config()
stopifnot(Sys.getenv("CBH_PRIMARY_VERSION","regional") %in% c("regional","regional18","regional_imputed","regional_mics"))
spec <- cbh_primary_regional_spec(settings);stopifnot(length(spec$covariates)==if(settings$nutrition)17L else 18L,"urban_pct" %in% spec$covariates)
for(p in c(settings$out,settings$private))dir.create(p,recursive=TRUE,showWarnings=FALSE)
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds"))
if(settings$imputed) {
  # Imputed-covariate sensitivity: overlay and national panel from the imputation
  # stage (covariates/15_impute_national.R, 16_impute_regional.R) and the extended
  # HIV panel. Their provenance files must match the files on disk.
  cs$private <- dirname(settings$imputed_overlay)
  cs$out <- settings$imputation_results
  provenance <- rbind(cbh_read_csv(file.path(cs$out,"regional_imputation_provenance.csv")),
    cbh_read_csv(file.path(cs$out,"national_imputation_provenance.csv")))
  stopifnot(meta$complete,identical(unname(vapply(provenance$file,cbh_file_hash,"")),provenance$md5),
    file.exists(settings$hiv_panel))
} else if(settings$nutrition) {
  cs$private <- "data/derived_cbh/regional_adjustment/planned17_audit"
  cs$out <- "results/cbh/planned17_covariate_missingness"
  provenance <- cbh_read_csv(file.path(cs$out,"provenance.csv"))
  stopifnot(meta$complete,identical(unname(vapply(provenance$file,cbh_file_hash,"")),provenance$md5))
} else {
  overlay_meta <- readRDS(file.path(cs$private,"manifest.rds"))
  stopifnot(meta$complete,overlay_meta$complete,identical(overlay_meta$specification$covariates,spec$covariates))
}
overlay_path <- if(settings$imputed) settings$imputed_overlay else file.path(cs$private,"regional_covariates_wide.csv")
hiv_panel_path <- if(is.null(settings$hiv_panel)) cbh_trial_spec()$incidence_panel else settings$hiv_panel
source_files <- c(file.path(cfg$output_dir,"manifest.rds"),
  if(settings$imputed) file.path(cs$out,c("regional_imputation_provenance.csv","national_imputation_provenance.csv"))
  else if(settings$nutrition) file.path(cs$out,"provenance.csv") else file.path(cs$private,"manifest.rds"),
  overlay_path,hiv_panel_path,
  if(settings$imputed) settings$imputed_national else file.path(cs$out,"complete_case_summary.csv"),
  "R_cbh/primary/00_prepare_regional.R","R_cbh/primary/specification.R","R_cbh/covariates/regional.R","R_cbh/analysis/model.R",
  if(settings$mics) c(file.path(settings$mics_output_dir,"manifest.rds"),settings$mics_overlay))
signature <- cbh_hash(list(files=source_files,hashes=vapply(source_files,cbh_file_hash,""),spec=spec))
signature_path <- file.path(settings$private,"preparation_signature.rds")
if(file.exists(settings$data) && file.exists(signature_path)) {
  cache <- readRDS(signature_path)
  if(identical(cache$signature,signature) && identical(cache$data_md5,cbh_file_hash(settings$data))) {
    message("Verified current regional prepared dataset; retaining it.");quit(save="no",status=0)
  }
}
wide <- cbh_read_csv(overlay_path)
cbh_unique(wide,c("survey","regkey"),"Regional overlay")
if(!settings$nutrition && !settings$imputed)stopifnot(identical(cbh_file_hash(overlay_path),overlay_meta$source_hashes$md5[match(overlay_path,overlay_meta$source_hashes$file)]))
hiv <- cbh_read_csv(hiv_panel_path)
if(settings$imputed) {
  national <- cbh_read_csv(settings$imputed_national); cbh_unique(national,c("iso3","year"),"Imputed national panel")
  flag_cols <- grep("_model_imputed$",names(wide),value=TRUE)
  imputed_counts <- list()
}
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
required <- c("death","age_band","pfpr_pct","calendar_year","band_years","survey","country","region",spec$covariates)
keep_cols <- unique(c(required,"child_id","entry_year","survey_year","regkey"))
pieces <- selections <- list()
for(i in seq_len(nrow(m))) {
  object <- readRDS(file.path(cfg$output_dir,m$file[i]));stopifnot(identical(object$signature,m$signature[i]))
  d <- cbh_attach_incidence(object$data,hiv)
  j <- match(paste(d$survey,d$regkey),paste(wide$survey,wide$regkey))
  for(v in spec$regional)d[[v]] <- wide[[v]][j]
  if(settings$imputed) {
    # National series with the 2001 WGI interpolation and health-expenditure model fills.
    k <- match(paste(d$country,d$entry_year),paste(national$iso3,national$year))
    d$political_stability <- national$political_stability[k]
    d$log_health_expenditure_pc <- national$log_health_expenditure_pc[k]
    ready <- d$model_ready
    imputed_counts[[i]] <- data.frame(survey=m$survey[i],country=m$country[i],records=sum(ready),
      regional_any_model_imputed=sum(rowSums(as.matrix(wide[j[ready],flag_cols,drop=FALSE]))>0),
      political_stability_interpolated=sum(national$political_stability_source[k[ready]]!="wgi_observed"),
      health_expenditure_imputed=sum(national$log_health_expenditure_pc_source[k[ready]]!="ghed_observed"),
      hiv_no_adolescent_series=sum(d$hiv_incidence_status[ready]=="imputed_no_adolescent_series"))
  }
  complete <- d$model_ready & complete.cases(d[required])
  for(v in required)if(is.numeric(d[[v]]))complete <- complete & is.finite(d[[v]])
  selections[[i]] <- data.frame(survey=m$survey[i],country=m$country[i],records=sum(complete),deaths=sum(d$death[complete]))
  pieces[[i]] <- d[complete,keep_cols,drop=FALSE]
  if(i%%20L==0L)message("Assembled ",i,"/",nrow(m)," survey shards")
}
if(settings$mics) {
  # The DHS part must reproduce the complete-case v3 sample exactly before MICS is added.
  dhs <- rbindlist(pieces)
  stopifnot(nrow(dhs)==settings$expected_dhs_records,sum(dhs$death)==settings$expected_dhs_deaths,
    uniqueN(dhs$region)==settings$expected_dhs_regions)
  message("DHS part reproduces v3: ",nrow(dhs)," records, ",sum(dhs$death)," deaths, ",uniqueN(dhs$region)," regions")
  rm(dhs)
  mmeta <- readRDS(file.path(settings$mics_output_dir,"manifest.rds")); stopifnot(mmeta$complete)
  mm <- mmeta$manifest[mmeta$manifest$status %in% c("built","cached"),]
  mwide <- cbh_read_csv(settings$mics_overlay); cbh_unique(mwide,c("survey","regkey"),"MICS regional overlay")
  for(i in seq_len(nrow(mm))) {
    object <- readRDS(file.path(settings$mics_output_dir,mm$file[i]));stopifnot(identical(object$signature,mm$signature[i]))
    d <- cbh_attach_incidence(object$data,hiv)
    j <- match(paste(d$survey,d$regkey),paste(mwide$survey,mwide$regkey))
    # Unlike the DHS join, a model-ready MICS record without an overlay row is an error, not a silent drop.
    if(any(d$model_ready & is.na(j))) stop("MICS overlay has no row for ",
      paste(unique(paste(d$survey,d$regkey)[d$model_ready & is.na(j)]),collapse=", "))
    for(v in spec$regional)d[[v]] <- mwide[[v]][j]
    complete <- d$model_ready & complete.cases(d[required])
    for(v in required)if(is.numeric(d[[v]]))complete <- complete & is.finite(d[[v]])
    selections[[length(selections)+1]] <- data.frame(survey=mm$survey[i],country=mm$country[i],records=sum(complete),deaths=sum(d$death[complete]))
    pieces[[length(pieces)+1]] <- d[complete,keep_cols,drop=FALSE]
  }
  message("Added ",nrow(mm)," MICS survey shards")
}
d <- as.data.frame(rbindlist(pieces));rm(pieces,object);gc(FALSE)
ages <- cfg$age_bands$age_band
d$age_band <- factor(d$age_band,levels=ages)
for(v in c("survey","country","region"))d[[v]] <- factor(d[[v]])
stopifnot(!anyNA(d[required]),all(d$death %in% 0:1),all(d$band_years>0),
  uniqueN(d,by=c("child_id","age_band"))==nrow(d),all(table(d$age_band)>0))
counts <- data.frame(records=nrow(d),children=uniqueN(d$child_id),deaths=sum(d$death),
  regions=nlevels(d$region),surveys=nlevels(d$survey),countries=nlevels(d$country))
if(settings$imputed) {
  # With every gap imputed, the sample must equal the MAP-eligible total: no record
  # may be lost to covariate availability.
  stopifnot(counts$records==settings$expected_records,counts$deaths==settings$expected_deaths,
    counts$regions==settings$expected_regions,counts$surveys==nrow(m),counts$countries==length(unique(m$country)))
  imputed_counts <- cbh_bind(imputed_counts)
  cbh_atomic_csv(imputed_counts,file.path(settings$out,"imputation_record_counts_by_survey.csv"))
  cbh_atomic_csv(data.frame(records=counts$records,
    regional_any_model_imputed=sum(imputed_counts$regional_any_model_imputed),
    political_stability_interpolated=sum(imputed_counts$political_stability_interpolated),
    health_expenditure_imputed=sum(imputed_counts$health_expenditure_imputed),
    hiv_no_adolescent_series=sum(imputed_counts$hiv_no_adolescent_series)),
    file.path(settings$out,"imputation_record_counts.csv"))
} else if(settings$mics) {
  cbh_atomic_csv(counts,file.path(settings$out,"prepared_sample_dhs_mics.csv"))
} else {
  expected <- cbh_read_csv(file.path(cs$out,"complete_case_summary.csv"))
  if(!settings$nutrition)expected <- expected[expected$baseline=="eligible_MAP",]
  for(v in names(counts))stopifnot(counts[[v]]==expected[[paste0("retained_",v)]])
}
scaling <- data.frame(variable=spec$covariates,mean=vapply(d[spec$covariates],mean,0),
  sd=vapply(d[spec$covariates],sd,0),row.names=NULL)
stopifnot(all(is.finite(scaling$sd) & scaling$sd>0))
for(i in seq_len(nrow(scaling))) {
  v <- scaling$variable[i];d[[paste0("z_",v)]] <- (d[[v]]-scaling$mean[i])/scaling$sd[i]
}
stopifnot("z_urban_pct" %in% names(d),all(is.finite(d$z_urban_pct)),
  !any(c("hib3_pct","pcv3_pct","rotavirus_pct","exclusive_breastfeeding_pct") %in% names(d)))
if(settings$nutrition)stopifnot(!any(c("male_pct","multiple_birth_pct","mean_birth_order","mean_maternal_age_birth") %in% names(d)),
  all(c("z_mean_maternal_age_first_birth","z_wasting_pct","z_stunting_pct") %in% names(d)))
prepared <- list(signature=signature,data=d,scaling=scaling,specification=spec,counts=counts,
  selection=cbh_bind(selections),source_hashes=data.frame(file=source_files,md5=vapply(source_files,cbh_file_hash,"")))
cbh_atomic_rds(prepared,settings$data)
cbh_atomic_rds(list(signature=signature,data_md5=cbh_file_hash(settings$data)),signature_path)
cbh_atomic_csv(counts,file.path(settings$out,"prepared_sample.csv"))
cbh_atomic_csv(scaling,file.path(settings$out,"covariate_scaling.csv"))
cbh_atomic_csv(prepared$selection,file.path(settings$out,"prepared_selection_by_survey.csv"))
cbh_atomic_csv(prepared$source_hashes,file.path(settings$out,"preparation_provenance.csv"))
print(counts)
message("Regional modelling dataset prepared and checked against the availability audit.")
