#!/usr/bin/env Rscript
# A bounded numerical sensitivity for fits with a nonpositive smoothing Hessian.
# Preserve original fits; select a tighter restart only on numerical criteria.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/model.R")
library(mgcv)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--retain-later-surveys"))
id <- if("--retain-later-surveys" %in% args) "age_band_snow_2000_2015_retrospective_v1" else "age_band_snow_2000_2015_v1"
private <- file.path("data/derived_cbh/models",id); out <- file.path("results/cbh",id)
diag <- cbh_read_csv(file.path(out,"fit_diagnostics.csv"))
manifest <- cbh_read_csv(file.path(out,"fit_manifest.csv"))
stopifnot(nrow(diag)==nrow(manifest))
curves <- cbh_read_csv(file.path(out,"pfpr_curves.csv"))
manifest$model_file <- paste0(manifest$fit_id,".rds")
manifest$selected_version <- "original"
selected_diag <- diag
jobs <- diag[is.finite(diag$min_smoothing_hessian_eigenvalue) & diag$min_smoothing_hessian_eigenvalue<=0,]
checks <- list()
for(i in seq_len(nrow(jobs))) {
  job <- jobs[i,];path <- file.path(private,paste0(job$fit_id,".rds"))
  original <- readRDS(path); f <- original$fit
  dd <- f$model; dd$band_years <- exp(dd[["offset(log(band_years))"]])
  dd$age_band <- factor(rep(job$age_band,nrow(dd)))
  signature <- cbh_hash(list(cbh_file_hash(path),cbh_file_hash("R_cbh/snow/05_check_smoothing.R")))
  cache <- file.path(private,"stability",paste0(job$fit_id,"_strict.rds"))
  saved <- if(file.exists(cache)) readRDS(cache) else NULL
  if(is.null(saved) || !identical(saved$signature,signature)) {
    warnings <- character();set.seed(20260907L)
    timing <- system.time(strict <- withCallingHandlers(bam(formula(f),data=dd,
      family=binomial(link="cloglog"),method="fREML",discrete=TRUE,nthreads=2,
      coef=coef(f),gc.level=1,na.action=na.fail,
      control=gam.control(epsilon=1e-9,mgcv.tol=1e-9,efs.tol=.001,maxit=200)),
      warning=function(w){warnings <<- unique(c(warnings,conditionMessage(w)));invokeRestart("muffleWarning")}))
    saved <- list(fit=strict,elapsed_seconds=unname(timing["elapsed"]),warnings=warnings,
      signature=signature,input_signature=original$input_signature,session_info=sessionInfo())
    cbh_atomic_rds(saved,cache);rm(strict)
  }
  detail <- cbh_sensitivity_diagnostics(saved,dd,paste0(job$fit_id,"_strict"))
  stopifnot(detail$converged,detail$finite_covariance,detail$finite_coefficients)
  a <- cbh_sensitivity_curves(f,dd,"original")
  b <- cbh_sensitivity_curves(saved$fit,dd,"strict")
  stopifnot(identical(a$pfpr_pct,b$pfpr_pct))
  use <- a$within_central_support
  difference <- max(abs(a$log_hazard_ratio[use]-b$log_hazard_ratio[use]))
  choose_strict <- detail$converged && detail$rank==detail$coefficients &&
    is.finite(detail$min_smoothing_hessian_eigenvalue) && detail$min_smoothing_hessian_eigenvalue>0 &&
    detail$max_smoothing_gradient<job$max_smoothing_gradient && saved$fit$gcv.ubre<=f$gcv.ubre
  checks[[i]] <- data.frame(fit_id=job$fit_id,original_min_hessian=job$min_smoothing_hessian_eigenvalue,
    strict_min_hessian=detail$min_smoothing_hessian_eigenvalue,strict_max_gradient=detail$max_smoothing_gradient,
    strict_converged=detail$converged,max_central_log_hr_change=difference,
    original_hr40to20=exp(-a$log_hazard_ratio[a$pfpr_pct==40]),
    strict_hr40to20=exp(-b$log_hazard_ratio[b$pfpr_pct==40]),elapsed_seconds=detail$elapsed_seconds,
    original_fREML=unname(f$gcv.ubre),strict_fREML=unname(saved$fit$gcv.ubre),strict_selected=choose_strict,
    original_fit_md5=cbh_file_hash(path),strict_fit_md5=cbh_file_hash(cache))
  cbh_atomic_csv(do.call(rbind,checks),file.path(out,"smoothing_stability.csv"))
  cbh_atomic_csv(rbind(a,b),file.path(out,paste0(job$fit_id,"_stability_curves.csv")))
  if(choose_strict) {
    j <- match(job$fit_id,manifest$fit_id)
    manifest$model_file[j] <- file.path("stability",paste0(job$fit_id,"_strict.rds"))
    manifest$selected_version[j] <- "strict_restart"
    b$fit_id <- job$fit_id; b$series <- job$series
    curves <- rbind(curves[curves$fit_id!=job$fit_id,],b)
    detail$fit_id <- job$fit_id; detail$series <- job$series; detail$age_band <- job$age_band
    selected_diag[selected_diag$fit_id==job$fit_id,] <- detail[names(selected_diag)]
  }
  print(checks[[i]][1:8],row.names=FALSE)
  rm(saved,original,f,dd);gc(FALSE)
}
manifest$model_file_md5 <- vapply(file.path(private,manifest$model_file),cbh_file_hash,"")
cbh_atomic_csv(manifest,file.path(out,"selected_fit_manifest.csv"))
cbh_atomic_csv(curves,file.path(out,"pfpr_curves_selected.csv"))
cbh_atomic_csv(selected_diag,file.path(out,"fit_diagnostics_selected.csv"))
inputs <- c("R_cbh/snow/05_check_smoothing.R",file.path(out,c("fit_manifest.csv","fit_diagnostics.csv","pfpr_curves.csv")))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"selection_provenance.csv"))
