#!/usr/bin/env Rscript
# Temporal accounting on the fixed previous-primary sample; no mortality fits.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/covariates/regional.R")
source("R_cbh/covariates/settings.R")
source("R_cbh/covariates/unicef.R")
library(data.table)
cfg <- cbh_config()
spec <- cbh_regional_spec(expanded=TRUE)
settings <- cbh_covariate_settings(expanded=TRUE)
previous <- cbh_covariate_settings(TRUE)
wide <- list(before=fread(file.path(previous$private,"regional_covariates_wide.csv")),
             after=fread(file.path(settings$private,"regional_covariates_wide.csv")))
panel <- cbh_read_csv(settings$panel)
hiv <- cbh_read_csv(cbh_trial_spec()$incidence_panel)
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds"))
stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
required <- c("death","age_band","pfpr_pct","calendar_year","band_years",
              "survey","country","region",cbh_trial_spec()$covariates)
pieces <- list()
for(i in seq_len(nrow(m))) {
  object <- readRDS(file.path(cfg$output_dir,m$file[i]))
  stopifnot(identical(object$signature,m$signature[i]))
  d <- cbh_attach_incidence(object$data,hiv)
  keep <- d$model_ready & complete.cases(d[required])
  for(v in required) if(is.numeric(d[[v]])) keep <- keep & is.finite(d[[v]])
  d <- d[keep,,drop=FALSE]
  if(!nrow(d)) next
  for(stage in names(wide)) {
    z <- d
    w <- wide[[stage]]
    j <- match(paste(z$survey,z$regkey),paste(w$survey,w$regkey))
    for(v in spec$regional) z[[v]] <- w[[v]][j]
    if(stage=="after") z <- cbh_unicef_fill(z,panel,c("hib3_pct","pcv3_pct","rotavirus_pct"))
    bad <- vapply(spec$covariates,function(v)!is.finite(z[[v]]),logical(nrow(z)))
    bad <- cbind(bad,any_required=rowSums(bad)>0)
    for(v in colnames(bad)) {
      a <- data.table(year=z$entry_year,missing=bad[,v])
      a <- a[,.(records=.N,missing_records=sum(missing)),by=year]
      a[,`:=`(stage=stage,variable=v)]
      pieces[[length(pieces)+1L]] <- a
    }
  }
  if(i%%20L==0L)message("Audited ",i,"/",nrow(m)," surveys")
}
annual <- rbindlist(pieces)[,.(records=sum(records),missing_records=sum(missing_records)),
                            by=.(stage,variable,year)]
annual[,missing_pct:=100*missing_records/records]
setorder(annual,stage,variable,year)
# Reconcile every covariate and the joint exclusion count with the full audits.
for(stage_name in names(wide)) {
  folder <- if(stage_name=="after")settings$out else previous$out
  reference <- fread(file.path(folder,"missingness_summary.csv"))[baseline=="previous_primary"]
  totals <- annual[stage==stage_name & variable!="any_required",
                   .(records=sum(records),missing_records=sum(missing_records)),by=variable]
  j <- match(totals$variable,reference$variable)
  stopifnot(!anyNA(j),all(totals$records==reference$records[j]),
            all(totals$missing_records==reference$missing_records[j]))
  joint <- fread(file.path(folder,"complete_case_summary.csv"))[baseline=="previous_primary"]
  a <- annual[stage==stage_name & variable=="any_required"]
  stopifnot(sum(a$records)==joint$records,
            sum(a$missing_records)==joint$records-joint$retained_records)
}
annual[,period:=as.character(cut(year,c(2000,2005,2010,2015,2020,2025),right=FALSE,
  labels=c("2000–2004","2005–2009","2010–2014","2015–2019","2020–2023")))]
stopifnot(!anyNA(annual$period),max(annual$year)==2023)
period <- annual[,.(records=sum(records),missing_records=sum(missing_records)),by=.(stage,variable,period)]
period[,missing_pct:=100*missing_records/records]
cbh_atomic_csv(annual,file.path(settings$out,"missingness_by_entry_year.csv"))
cbh_atomic_csv(period,file.path(settings$out,"missingness_by_period.csv"))
print(dcast(period[stage=="after" & variable %in% c(names(cbh_unicef_indicators()),"any_required")],
            period+records~variable,value.var="missing_pct"))
