#!/usr/bin/env Rscript
# Two separate Snow-only changes: gamma=1.4, or a PfPR shrinkage cubic basis.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/model.R")
library(mgcv)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--force"))
id <- "age_band_snow_2000_2015_v1"
base_private <- file.path("data/derived_cbh/models",id)
base_out <- file.path("results/cbh",id)
private <- file.path(base_private,"penalty_sensitivity")
out <- file.path(base_out,"penalty_sensitivity")
dir.create(private,recursive=TRUE,showWarnings=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
manifest <- cbh_read_csv(file.path(base_out,"selected_fit_manifest.csv"))
manifest <- manifest[manifest$series=="snow",]
stopifnot(nrow(manifest)==7)
specs <- data.frame(series=c("reference","gamma14","cs"),basis=c("cr","cr","cs"),gamma=c(1,1.4,1))
cbh_atomic_csv(specs,file.path(out,"specifications.csv"))
code <- c("R_cbh/snow/06_fit_penalties.R","R_cbh/sensitivity/model.R")
code_hash <- vapply(code,cbh_file_hash,"")
change_pfpr_basis <- function(x) {
  if(is.call(x)) {
    if(identical(x[[1]],as.name("s")) && identical(x[[2]],as.name("pfpr_pct"))) x$bs <- "cs"
    else for(j in seq_along(x)[-1L]) x[[j]] <- change_pfpr_basis(x[[j]])
  }
  x
}
run_fit <- function(form,d,gamma,knots,initial=NULL,strict=FALSE) {
  warnings <- character();set.seed(20260907L)
  control <- if(strict) gam.control(epsilon=1e-9,mgcv.tol=1e-9,efs.tol=.001,maxit=200) else
    gam.control(maxit=100)
  timing <- system.time(f <- withCallingHandlers(bam(form,data=d,knots=knots,
    family=binomial(link="cloglog"),method="fREML",discrete=TRUE,gamma=gamma,select=FALSE,
    nthreads=2,gc.level=1,na.action=na.fail,coef=initial,control=control),
    warning=function(w){warnings <<- unique(c(warnings,conditionMessage(w)));invokeRestart("muffleWarning")}))
  list(fit=f,elapsed_seconds=unname(timing["elapsed"]),warnings=warnings,session_info=sessionInfo())
}
curves <- diagnostics <- summaries <- selections <- restarts <- list()
for(i in seq_len(nrow(manifest))) {
  base_path <- file.path(base_private,manifest$model_file[i])
  stopifnot(identical(cbh_file_hash(base_path),manifest$model_file_md5[i]))
  baseline <- readRDS(base_path)
  d <- baseline$fit$model
  d$band_years <- exp(d[["offset(log(band_years))"]])
  d$age_band <- factor(rep(manifest$age_band[i],nrow(d)))
  form <- formula(baseline$fit)
  knots <- setNames(lapply(baseline$fit$smooth[1:2],function(s)s$xp),c("pfpr_pct","calendar_year"))
  stopifnot(length(knots$pfpr_pct)==5,length(knots$calendar_year)==6,
    all(floor(d$calendar_year) %in% 2000:2015),!anyNA(d))
  for(k in seq_len(nrow(specs))) {
    job <- specs[k,]; fit_id <- paste0(job$series,"_age_",i)
    ff <- form
    if(job$basis=="cs") ff[[3]] <- change_pfpr_basis(ff[[3]])
    signature <- cbh_hash(list(manifest$model_file_md5[i],job,knots,deparse(ff),code_hash,
      R.version.string,as.character(packageVersion("mgcv"))))
    selected_version <- "original"
    path <- if(job$series=="reference") base_path else file.path(private,paste0(fit_id,".rds"))
    message("Starting ",fit_id,": ",nrow(d)," rows, ",sum(d$death)," deaths")
    saved <- if(job$series=="reference") baseline else if(file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
    if(job$series!="reference" && (is.null(saved) || !identical(saved$signature,signature))) {
      saved <- run_fit(ff,d,job$gamma,knots)
      saved$signature <- signature; saved$input_signature <- baseline$input_signature
      saved$specification <- job; saved$knots <- knots
      cbh_atomic_rds(saved,path)
    }
    diag <- cbh_sensitivity_diagnostics(saved,d,fit_id)
    if(job$series!="reference" && (!diag$converged ||
      (is.finite(diag$min_smoothing_hessian_eigenvalue) && diag$min_smoothing_hessian_eigenvalue<=0))) {
      strict_path <- file.path(private,paste0(fit_id,"_strict.rds"))
      strict_sig <- cbh_hash(list(signature,cbh_file_hash(path),"strict-restart"))
      strict <- if(file.exists(strict_path) && !"--force" %in% args) readRDS(strict_path) else NULL
      if(is.null(strict) || !identical(strict$signature,strict_sig)) {
        strict <- run_fit(ff,d,job$gamma,knots,coef(saved$fit),TRUE)
        strict$signature <- strict_sig; strict$input_signature <- baseline$input_signature
        strict$specification <- job;strict$knots <- knots
        cbh_atomic_rds(strict,strict_path)
      }
      sd <- cbh_sensitivity_diagnostics(strict,d,fit_id)
      take <- sd$converged && sd$finite_covariance && sd$rank==sd$coefficients &&
        is.finite(sd$min_smoothing_hessian_eigenvalue) && sd$min_smoothing_hessian_eigenvalue>0 &&
        sd$max_smoothing_gradient<diag$max_smoothing_gradient && strict$fit$gcv.ubre<=saved$fit$gcv.ubre
      restarts[[fit_id]] <- data.frame(fit_id=fit_id,original_min_hessian=diag$min_smoothing_hessian_eigenvalue,
        strict_min_hessian=sd$min_smoothing_hessian_eigenvalue,original_fREML=unname(saved$fit$gcv.ubre),
        strict_fREML=unname(strict$fit$gcv.ubre),strict_selected=take)
      cbh_atomic_csv(do.call(rbind,restarts),file.path(out,"restarts.csv"))
      if(take) {saved <- strict;diag <- sd;path <- strict_path;selected_version <- "strict_restart"}
      rm(strict)
    }
    stopifnot(diag$converged,diag$finite_coefficients,diag$finite_covariance,diag$rank==diag$coefficients,
      nrow(saved$fit$model)==nrow(d))
    # The original Snow model frame is the input. Verify every stored predictor
    # and outcome after fitting; no MAP-availability restriction is introduced.
    for(v in names(baseline$fit$model)) stopifnot(isTRUE(all.equal(
      saved$fit$model[[v]],baseline$fit$model[[v]],check.attributes=FALSE)))
    stopifnot(inherits(saved$fit$smooth[[1]],if(job$basis=="cs") "cs.smooth" else "cr.smooth"),
      inherits(saved$fit$smooth[[2]],"cr.smooth"),
      identical(saved$fit$smooth[[1]]$xp,knots$pfpr_pct),
      identical(saved$fit$smooth[[2]]$xp,knots$calendar_year))
    cc <- cbh_sensitivity_curves(saved$fit,d,fit_id);cc$series <- job$series
    curves[[fit_id]] <- cc
    diag$series <- job$series;diag$age_band <- manifest$age_band[i]
    diag$pfpr_basis <- job$basis;diag$gamma <- job$gamma;diag$input_verified <- TRUE
    diagnostics[[fit_id]] <- diag
    st <- as.data.frame(summary(saved$fit)$s.table);st$term <- rownames(st)
    st$series <- job$series;st$age_band <- manifest$age_band[i];st$fit_id <- fit_id
    summaries[[fit_id]] <- st
    selections[[fit_id]] <- data.frame(fit_id=fit_id,series=job$series,age_band=manifest$age_band[i],
      selected_version=selected_version,model_file=path,md5=cbh_file_hash(path),
      input_signature=baseline$input_signature)
    cbh_atomic_csv(do.call(rbind,curves),file.path(out,"pfpr_curves.csv"))
    cbh_atomic_csv(do.call(rbind,diagnostics),file.path(out,"fit_diagnostics.csv"))
    cbh_atomic_csv(do.call(rbind,summaries),file.path(out,"smooth_summaries.csv"))
    cbh_atomic_csv(do.call(rbind,selections),file.path(out,"fit_manifest.csv"))
    writeLines(trimws(deparse(ff),which="right"),file.path(out,paste0(job$series,"_formula.txt")))
    message("Completed ",fit_id,"; fit time ",round(diag$elapsed_seconds,1)," seconds")
    if(job$series!="reference") rm(saved)
    gc(FALSE)
  }
  rm(baseline,d);gc(FALSE)
}
stopifnot(length(diagnostics)==21L)
cbh_atomic_csv(data.frame(file=code,md5=unname(code_hash)),file.path(out,"code_provenance.csv"))
message("Snow penalty sensitivities complete")
