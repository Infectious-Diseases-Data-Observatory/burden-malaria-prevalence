#!/usr/bin/env Rscript
# Aggregate in-sample outcome checks; no row-level predictions are exported.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
library(data.table)
settings <- cbh_primary_settings();out <- settings$out
manifest <- cbh_read_csv(file.path(out,"fit_manifest.csv"))
provenance <- cbh_read_csv(file.path(out,"fit_input_provenance.csv"))
stopifnot(nrow(manifest)==7,all(manifest$series=="map_full"),
  identical(cbh_file_hash(settings$data),provenance$md5[provenance$file==settings$data]))
checks <- groups <- coverage <- list()
for(i in seq_len(nrow(manifest))) {
  m <- manifest[i,]
  stopifnot(identical(cbh_file_hash(m$model_file),m$md5))
  saved <- readRDS(m$model_file);f <- saved$fit
  coverage[[i]] <- unique(f$model[c("survey","country","region")])
  p <- f$fitted.values;y <- f$model$death
  stopifnot(saved$gamma==2,all(is.finite(p)),all(p>0 & p<1),length(p)==nrow(f$model))
  checks[[i]] <- data.frame(age_band=m$age_band,records=length(y),observed_deaths=sum(y),
    fitted_deaths=sum(p),mean_observed_probability=mean(y),mean_fitted_probability=mean(p),
    brier_score=mean((y-p)^2),bernoulli_log_loss=-mean(y*log(p)+(1-y)*log1p(-p)),
    smoothing_corrected_covariance_available=!is.null(f$Vc),model_md5_verified=TRUE)
  # Aggregate residuals flag patterns; they are not tests of survey influence or held-out performance.
  for(v in c("country","survey","calendar_year","pfpr_pct","fitted_probability")) {
    g <- switch(v,pfpr_pct=floor(f$model$pfpr_pct/10)*10,
      calendar_year=floor(f$model$calendar_year),
      fitted_probability=cut(p,unique(quantile(p,seq(0,1,.1))),include.lowest=TRUE),f$model[[v]])
    dt <- data.table(group=as.character(g),observed=y,predicted=p)
    a <- dt[,.(records=.N,observed_deaths=sum(observed),fitted_deaths=sum(predicted)),by=group]
    a[,`:=`(age_band=m$age_band,grouping=v,observed_probability=observed_deaths/records,
      fitted_probability=fitted_deaths/records,residual_deaths=observed_deaths-fitted_deaths)]
    # Small cells are suppressed, including their observed/fitted summaries.
    a <- a[records>=100]
    groups[[paste(i,v)]] <- as.data.frame(a)
  }
  message("Verified stored model and fitted outcomes for ",m$age_band," months")
  rm(saved,f,p,y,dt);gc(FALSE)
}
coverage <- unique(rbindlist(coverage))
coverage <- coverage[,.(regions=uniqueN(region)),by=.(survey,country)]
sample <- cbh_read_csv(file.path(out,"primary_sample.csv"))
stopifnot(nrow(coverage)==sample$surveys,uniqueN(coverage$country)==sample$countries)
cbh_atomic_csv(as.data.frame(coverage),file.path(out,"survey_coverage.csv"))
cbh_atomic_csv(do.call(rbind,checks),file.path(out,"fitted_outcome_checks.csv"))
cbh_atomic_csv(do.call(rbind,groups),file.path(out,"grouped_outcome_checks.csv"))
cbh_atomic_csv(data.frame(file=c(manifest$model_file,"R_cbh/primary/04_diagnostics.R"),
  md5=vapply(c(manifest$model_file,"R_cbh/primary/04_diagnostics.R"),cbh_file_hash,"")),
  file.path(out,"diagnostic_provenance.csv"))
