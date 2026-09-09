#!/usr/bin/env Rscript
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
library(mgcv)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% c("--force","--retain-later-surveys")))
id <- if ("--retain-later-surveys" %in% args) "age_band_snow_2000_2015_retrospective_v1" else "age_band_snow_2000_2015_v1"
private <- file.path("data/derived_cbh/models",id); out <- file.path("results/cbh",id)
prepared <- readRDS(file.path(private,"dataset.rds"))
spec <- prepared$spec
form <- cbh_sensitivity_formula(spec,TRUE)
input_signature <- prepared$signature
vars <- unique(c(all.vars(form),"age_band","map_pfpr_pct","snow_pfpr_pct"))
d <- prepared$data[vars]
rm(prepared); gc(FALSE)
age <- levels(d$age_band)
matched <- is.finite(d$map_pfpr_pct)
# Fit a matched Snow comparator only if MAP is missing among eligible Snow rows.
series <- c("snow", "map_matched", if (any(!matched)) "snow_matched")
jobs <- expand.grid(age_band=age,series=series,stringsAsFactors=FALSE)
jobs$fit_id <- paste0(jobs$series,"_age_",match(jobs$age_band,age))
cbh_atomic_csv(jobs,file.path(out,"fit_manifest.csv"))
code <- c("R_cbh/snow/02_fit.R","R_cbh/sensitivity/model.R","R_cbh/analysis/model.R")
signature <- cbh_hash(list(input_signature,vapply(code,cbh_file_hash,""),deparse(form),
  R.version.string,as.character(packageVersion("mgcv"))))
diags <- curves <- counts <- list()
for (k in seq_len(nrow(jobs))) {
  job <- jobs[k,]; fit_id <- job$fit_id
  use <- d$age_band==job$age_band & if (job$series=="snow") TRUE else matched
  dd <- droplevels(d[use,,drop=FALSE])
  dd$pfpr_pct <- if (job$series=="map_matched") dd$map_pfpr_pct else dd$snow_pfpr_pct
  dd <- dd[unique(c(all.vars(form),"age_band"))]
  stopifnot(nrow(dd)>0,!anyNA(dd),nlevels(dd$age_band)==1L)
  path <- file.path(private,paste0(fit_id,".rds"))
  sig <- cbh_hash(list(signature,job))
  saved <- if (file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
  message("Starting ",fit_id,": ",nrow(dd)," records, ",sum(dd$death)," deaths")
  if (is.null(saved) || !identical(saved$signature,sig)) {
    warnings <- character(); set.seed(spec$seed)
    timing <- system.time(fit <- withCallingHandlers(bam(form,data=dd,
      family=binomial(link="cloglog"),method="fREML",discrete=TRUE,nthreads=spec$nthreads,
      gc.level=1,na.action=na.fail,control=gam.control(trace=FALSE,maxit=100)),
      warning=function(w){warnings <<- unique(c(warnings,conditionMessage(w)));invokeRestart("muffleWarning")}))
    saved <- list(fit=fit,elapsed_seconds=unname(timing["elapsed"]),warnings=warnings,
      signature=sig,input_signature=input_signature,session_info=sessionInfo())
    cbh_atomic_rds(saved,path); rm(fit)
  }
  diag <- cbh_sensitivity_diagnostics(saved,dd,fit_id)
  diag$series <- job$series; diag$age_band <- job$age_band
  diags[[fit_id]] <- diag
  cbh_atomic_csv(do.call(rbind,diags),file.path(out,"fit_diagnostics.csv"))
  stopifnot(diag$converged,diag$finite_coefficients,diag$finite_covariance,
    diag$rank==diag$coefficients,length(saved$fit$fitted.values)==nrow(dd))
  cc <- cbh_sensitivity_curves(saved$fit,dd,fit_id); cc$series <- job$series
  curves[[fit_id]] <- cc
  cbh_atomic_csv(do.call(rbind,curves),file.path(out,"pfpr_curves.csv"))
  st <- as.data.frame(summary(saved$fit)$s.table); st$term <- rownames(st)
  cbh_atomic_csv(st,file.path(out,paste0(fit_id,"_smooth_summary.csv")))
  message("Completed ",fit_id," in ",round(diag$elapsed_seconds,1)," fit seconds")
  rm(saved,dd); gc(FALSE)
}
stopifnot(length(diags)==nrow(jobs))
writeLines(trimws(deparse(form),which="right"),file.path(out,"model_formula.txt"))
cbh_atomic_csv(data.frame(file=code,md5=vapply(code,cbh_file_hash,"")),file.path(out,"fitting_provenance.csv"))
message("All Snow and matched MAP fits complete")
