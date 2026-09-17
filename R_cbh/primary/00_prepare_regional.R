#!/usr/bin/env Rscript
# Assemble the audited regional overlay; no source extraction or HIV refit.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/covariates/regional.R")
source("R_cbh/covariates/settings.R")
source("R_cbh/primary/settings.R")
library(data.table)
settings <- cbh_primary_settings("regional");cs <- cbh_covariate_settings();cfg <- cbh_config()
spec <- cbh_regional_spec();stopifnot(length(spec$covariates)==18,"urban_pct" %in% spec$covariates)
for(p in c(settings$out,settings$private))dir.create(p,recursive=TRUE,showWarnings=FALSE)
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds"))
overlay_meta <- readRDS(file.path(cs$private,"manifest.rds"))
stopifnot(meta$complete,overlay_meta$complete,identical(overlay_meta$specification$covariates,spec$covariates))
overlay_path <- file.path(cs$private,"regional_covariates_wide.csv")
source_files <- c(file.path(cfg$output_dir,"manifest.rds"),file.path(cs$private,"manifest.rds"),
  overlay_path,cbh_trial_spec()$incidence_panel,file.path(cs$out,"complete_case_summary.csv"),
  "R_cbh/primary/00_prepare_regional.R","R_cbh/covariates/regional.R","R_cbh/analysis/model.R")
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
stopifnot(identical(cbh_file_hash(overlay_path),overlay_meta$source_hashes$md5[match(overlay_path,overlay_meta$source_hashes$file)]))
hiv <- cbh_read_csv(cbh_trial_spec()$incidence_panel)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
required <- c("death","age_band","pfpr_pct","calendar_year","band_years","survey","country","region",spec$covariates)
keep_cols <- unique(c(required,"child_id","entry_year","survey_year","regkey"))
pieces <- selections <- list()
for(i in seq_len(nrow(m))) {
  object <- readRDS(file.path(cfg$output_dir,m$file[i]));stopifnot(identical(object$signature,m$signature[i]))
  d <- cbh_attach_incidence(object$data,hiv)
  j <- match(paste(d$survey,d$regkey),paste(wide$survey,wide$regkey))
  for(v in spec$regional)d[[v]] <- wide[[v]][j]
  complete <- d$model_ready & complete.cases(d[required])
  for(v in required)if(is.numeric(d[[v]]))complete <- complete & is.finite(d[[v]])
  selections[[i]] <- data.frame(survey=m$survey[i],country=m$country[i],records=sum(complete),deaths=sum(d$death[complete]))
  pieces[[i]] <- d[complete,keep_cols,drop=FALSE]
  if(i%%20L==0L)message("Assembled ",i,"/",nrow(m)," survey shards")
}
d <- as.data.frame(rbindlist(pieces));rm(pieces,object);gc(FALSE)
ages <- cfg$age_bands$age_band
d$age_band <- factor(d$age_band,levels=ages)
for(v in c("survey","country","region"))d[[v]] <- factor(d[[v]])
stopifnot(!anyNA(d[required]),all(d$death %in% 0:1),all(d$band_years>0),
  uniqueN(d,by=c("child_id","age_band"))==nrow(d),all(table(d$age_band)>0))
counts <- data.frame(records=nrow(d),children=uniqueN(d$child_id),deaths=sum(d$death),
  regions=nlevels(d$region),surveys=nlevels(d$survey),countries=nlevels(d$country))
expected <- cbh_read_csv(file.path(cs$out,"complete_case_summary.csv"))
expected <- expected[expected$baseline=="eligible_MAP",]
for(v in names(counts))stopifnot(counts[[v]]==expected[[paste0("retained_",v)]])
scaling <- data.frame(variable=spec$covariates,mean=vapply(d[spec$covariates],mean,0),
  sd=vapply(d[spec$covariates],sd,0),row.names=NULL)
stopifnot(all(is.finite(scaling$sd) & scaling$sd>0))
for(i in seq_len(nrow(scaling))) {
  v <- scaling$variable[i];d[[paste0("z_",v)]] <- (d[[v]]-scaling$mean[i])/scaling$sd[i]
}
stopifnot("z_urban_pct" %in% names(d),all(is.finite(d$z_urban_pct)),
  !any(c("hib3_pct","pcv3_pct","rotavirus_pct","exclusive_breastfeeding_pct") %in% names(d)))
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
