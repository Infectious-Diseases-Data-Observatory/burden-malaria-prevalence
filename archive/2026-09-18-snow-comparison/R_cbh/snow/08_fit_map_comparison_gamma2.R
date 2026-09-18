#!/usr/bin/env Rscript
# Compare exposure sources at gamma=2, with matched and full fitting samples.
source("R_cbh/load_pipeline.R")
source("R_cbh/sensitivity/model.R")
library(mgcv)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--force"))
id <- "map_snow_gamma2_v1"
private <- file.path("data/derived_cbh/models",id)
out <- file.path("results/cbh",id)
dir.create(private,recursive=TRUE,showWarnings=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ages <- cbh_config()$age_bands$age_band
snow_base <- "data/derived_cbh/models/age_band_snow_2000_2015_v1"
old <- cbh_read_csv("results/cbh/age_band_snow_2000_2015_v1/selected_fit_manifest.csv")
old <- old[old$series %in% c("map_matched","snow_matched"),]
jobs <- data.frame(series=old$series,age_band=old$age_band,
  source_path=file.path(snow_base,old$model_file),source_md5=old$model_file_md5,reuse=FALSE)
prev <- cbh_read_csv("results/cbh/age_band_separate_v1/country_burden_2024/provenance.csv")
paths <- file.path("data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/sensitivity_single_imputation",paste0("age_",seq_along(ages),".rds"))
stopifnot(all(paths %in% prev$file))
jobs <- rbind(jobs,data.frame(series="map_full",age_band=ages,source_path=paths,
  source_md5=prev$md5[match(paths,prev$file)],reuse=FALSE))
sn <- cbh_read_csv("results/cbh/age_band_snow_2000_2015_v1/penalty_sensitivity/gamma2/fit_manifest.csv")
sn <- sn[sn$series=="gamma2",]
jobs <- rbind(jobs,data.frame(series="snow_full",age_band=sn$age_band,source_path=sn$model_file,
  source_md5=sn$md5,reuse=TRUE))
jobs$fit_id <- paste0(jobs$series,"_age_",match(jobs$age_band,ages))
jobs <- jobs[order(match(jobs$age_band,ages),match(jobs$series,c("map_matched","snow_matched","map_full","snow_full"))),]
cbh_atomic_csv(jobs,file.path(out,"jobs.csv"))
code <- c("R_cbh/snow/08_fit_map_comparison_gamma2.R","R_cbh/sensitivity/model.R")
code_hash <- vapply(code,cbh_file_hash,"")
run <- function(form,d,knots,initial=NULL,strict=FALSE) {
  warnings <- character();set.seed(20260907L)
  control <- if(strict) gam.control(epsilon=1e-9,mgcv.tol=1e-9,efs.tol=.001,maxit=200) else gam.control(maxit=100)
  timing <- system.time(f <- withCallingHandlers(bam(form,data=d,knots=knots,
    family=binomial(link="cloglog"),method="fREML",discrete=TRUE,gamma=2,select=FALSE,
    nthreads=2,gc.level=1,na.action=na.fail,coef=initial,control=control),
    warning=function(w){warnings <<- unique(c(warnings,conditionMessage(w)));invokeRestart("muffleWarning")}))
  list(fit=f,elapsed_seconds=unname(timing["elapsed"]),warnings=warnings,session_info=sessionInfo())
}
diags <- curves <- summaries <- selections <- restarts <- components <- list()
for(i in seq_len(nrow(jobs))) {
  job <- jobs[i,];fit_id <- job$fit_id
  stopifnot(identical(cbh_file_hash(job$source_path),job$source_md5))
  base <- readRDS(job$source_path)
  d <- base$fit$model
  d$band_years <- exp(d[["offset(log(band_years))"]])
  d$age_band <- factor(rep(job$age_band,nrow(d)))
  stopifnot(!anyNA(d),max(abs(d$band_years-c(1,5,6,12,12,12,12)[match(job$age_band,ages)]/12))<1e-12)
  if(job$series!="map_full") stopifnot(all(floor(d$calendar_year) %in% 2000:2015))
  form <- formula(base$fit)
  knots <- setNames(lapply(base$fit$smooth[1:2],function(s)s$xp),c("pfpr_pct","calendar_year"))
  stopifnot(length(knots$pfpr_pct)==5,length(knots$calendar_year)==6)
  nuisance <- setdiff(names(base$fit$model),"pfpr_pct")
  nuisance_hash <- cbh_hash(lapply(base$fit$model[nuisance],as.vector))
  if(job$series=="map_matched") paired_nuisance_hash <- nuisance_hash
  if(job$series=="snow_matched") stopifnot(identical(nuisance_hash,paired_nuisance_hash))
  sig <- cbh_hash(list(job,knots,deparse(form),code_hash,R.version.string,as.character(packageVersion("mgcv"))))
  path <- if(job$reuse) job$source_path else file.path(private,paste0(fit_id,".rds"))
  message("Starting ",fit_id,": ",nrow(d)," rows, ",sum(d$death)," deaths")
  saved <- if(job$reuse) base else if(file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
  if(!job$reuse && (is.null(saved) || !identical(saved$signature,sig))) {
    saved <- run(form,d,knots)
    saved$signature <- sig;saved$input_signature <- base$input_signature
    saved$gamma <- 2;saved$knots <- knots
    cbh_atomic_rds(saved,path)
  }
  diag <- cbh_sensitivity_diagnostics(saved,d,fit_id)
  version <- if(job$reuse) "reused_gamma2" else "original"
  if(!job$reuse && (!diag$converged || !is.finite(diag$min_smoothing_hessian_eigenvalue) || diag$min_smoothing_hessian_eigenvalue<=0)) {
    spath <- file.path(private,paste0(fit_id,"_strict.rds"))
    ssig <- cbh_hash(list(sig,cbh_file_hash(path),"strict"))
    strict <- if(file.exists(spath) && !"--force" %in% args) readRDS(spath) else NULL
    if(is.null(strict) || !identical(strict$signature,ssig)) {
      strict <- run(form,d,knots,coef(saved$fit),TRUE)
      strict$signature <- ssig;strict$input_signature <- base$input_signature;strict$gamma <- 2;strict$knots <- knots
      cbh_atomic_rds(strict,spath)
    }
    sd <- cbh_sensitivity_diagnostics(strict,d,fit_id)
    take <- sd$converged && sd$finite_covariance && sd$rank==sd$coefficients &&
      is.finite(sd$min_smoothing_hessian_eigenvalue) && sd$min_smoothing_hessian_eigenvalue>0 &&
      sd$max_smoothing_gradient<diag$max_smoothing_gradient && strict$fit$gcv.ubre<=saved$fit$gcv.ubre
    restarts[[fit_id]] <- data.frame(fit_id,original_min_hessian=diag$min_smoothing_hessian_eigenvalue,
      strict_min_hessian=sd$min_smoothing_hessian_eigenvalue,original_fREML=unname(saved$fit$gcv.ubre),
      strict_fREML=unname(strict$fit$gcv.ubre),strict_selected=take)
    cbh_atomic_csv(do.call(rbind,restarts),file.path(out,"restarts.csv"))
    if(take) {saved <- strict;diag <- sd;path <- spath;version <- "strict_restart"}
    rm(strict)
  }
  f <- saved$fit
  stopifnot(diag$converged,diag$finite_coefficients,diag$finite_covariance,diag$rank==diag$coefficients,
    inherits(f$smooth[[1]],"cr.smooth"),inherits(f$smooth[[2]],"cr.smooth"),
    identical(f$smooth[[1]]$xp,knots$pfpr_pct),identical(f$smooth[[2]]$xp,knots$calendar_year))
  for(v in names(base$fit$model)) stopifnot(isTRUE(all.equal(f$model[[v]],base$fit$model[[v]],check.attributes=FALSE)))
  cc <- cbh_sensitivity_curves(f,d,fit_id);cc$series <- job$series;curves[[fit_id]] <- cc
  diag$series <- job$series;diag$age_band <- job$age_band;diag$gamma <- 2;diag$input_verified <- TRUE
  diag$nuisance_hash <- nuisance_hash;diags[[fit_id]] <- diag
  st <- as.data.frame(summary(f)$s.table);st$term <- rownames(st);st$fit_id <- fit_id
  st$series <- job$series;st$age_band <- job$age_band;summaries[[fit_id]] <- st
  # Compact PfPR basis, checked against full predictions, for later burden evaluation.
  s <- f$smooth[[1]];ix <- s$first.para:s$last.para
  nd <- d[rep(1L,5L),];nd$pfpr_pct <- c(0,20,40,60,max(d$pfpr_pct));ref <- nd;ref$pfpr_pct <- 20
  L <- predict(f,nd,type="lpmatrix",discrete=FALSE)-predict(f,ref,type="lpmatrix",discrete=FALSE)
  C <- PredictMat(s,nd)-PredictMat(s,ref)
  stopifnot(max(abs(drop(L%*%coef(f))-drop(C%*%coef(f)[ix])))<1e-8,
    max(abs(rowSums((L%*%f$Vp)*L)-rowSums((C%*%f$Vp[ix,ix])*C)))<1e-8)
  components[[fit_id]] <- list(fit_id=fit_id,series=job$series,age_band=job$age_band,
    smooth=s,coef=coef(f)[ix],covariance=f$Vp[ix,ix],model_md5=cbh_file_hash(path),
    support=quantile(d$pfpr_pct,c(0,.025,.975,1),names=FALSE))
  selections[[fit_id]] <- data.frame(fit_id,series=job$series,age_band=job$age_band,model_file=path,
    md5=cbh_file_hash(path),selected_version=version,source_path=job$source_path,source_md5=job$source_md5)
  cbh_atomic_csv(do.call(rbind,diags),file.path(out,"fit_diagnostics.csv"))
  cbh_atomic_csv(do.call(rbind,curves),file.path(out,"pfpr_curves.csv"))
  cbh_atomic_csv(do.call(rbind,summaries),file.path(out,"smooth_summaries.csv"))
  cbh_atomic_csv(do.call(rbind,selections),file.path(out,"fit_manifest.csv"))
  cbh_atomic_rds(components,file.path(private,"pfpr_components.rds"))
  writeLines(deparse(form),file.path(out,"model_formula.txt"))
  message("Completed ",fit_id,"; ",round(diag$elapsed_seconds,1)," fit seconds; Hessian minimum ",signif(diag$min_smoothing_hessian_eigenvalue,3))
  rm(base,saved,f,d);gc(FALSE)
}
stopifnot(length(diags)==28L)
cbh_atomic_csv(data.frame(file=code,md5=unname(code_hash)),file.path(out,"code_provenance.csv"))
message("MAP/Snow gamma=2 fits complete")
