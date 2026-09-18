#!/usr/bin/env Rscript
# Refit only the primary MAP models from the saved prepared dataset.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
library(mgcv)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--force"))
settings <- cbh_primary_settings()
if(settings$regional) source("R_cbh/primary/specification.R")
id <- settings$id; private <- settings$private; out <- settings$out
for(path in c(private,out)) dir.create(path,recursive=TRUE,showWarnings=FALSE)
ages <- cbh_config()$age_bands$age_band
jobs <- cbh_read_csv(file.path(settings$reference,"fit_manifest.csv"))
jobs <- jobs[jobs$series=="map_full",]
jobs <- jobs[match(ages,jobs$age_band),]
jobs$source_path <- jobs$model_file; jobs$source_md5 <- jobs$md5; jobs$reuse <- FALSE
stopifnot(nrow(jobs)==7L,!anyNA(jobs))
# The committed knot snapshot preserves the promoted MAP bases without fitting
# any reference or supplementary model. On first migration verify the old fits.
if(!file.exists(settings$knots)) {
  kk <- lapply(seq_len(nrow(jobs)),function(i) {
    stopifnot(identical(cbh_file_hash(jobs$source_path[i]),jobs$source_md5[i]))
    old <- readRDS(jobs$source_path[i])
    do.call(rbind,lapply(names(old$knots),function(v)
      data.frame(age_band=jobs$age_band[i],variable=v,index=seq_along(old$knots[[v]]),value=old$knots[[v]])))
  })
  cbh_atomic_csv(do.call(rbind,kk),settings$knots)
}
knot_table <- cbh_read_csv(settings$knots)
message("Loading existing prepared sample; no data setup or HIV refitting")
prepared <- readRDS(settings$data)
input_signature <- prepared$signature
input <- prepared$data
expected <- c(settings$expected_records,settings$expected_deaths)
stopifnot(nrow(input)==expected[1],sum(input$death)==expected[2],
  !anyNA(input$child_id),data.table::uniqueN(input,by=c("child_id","age_band"))==nrow(input))
cbh_atomic_csv(data.frame(distinct_children=data.table::uniqueN(input$child_id),records=nrow(input),
  deaths=sum(input$death),surveys=data.table::uniqueN(input$survey),countries=data.table::uniqueN(input$country)),
  file.path(out,"primary_sample.csv"))
cbh_atomic_csv(prepared$scaling,file.path(out,"covariate_scaling.csv"))
form <- if(settings$regional) cbh_primary_regional_formula(settings) else cbh_sensitivity_formula(cbh_trial_spec(),single_age=TRUE)
if(settings$regional) stopifnot(nrow(prepared$scaling)==length(cbh_primary_regional_spec(settings)$covariates),
  identical(prepared$scaling$variable,cbh_primary_regional_spec(settings)$covariates),
  "z_urban_pct" %in% all.vars(form))
input <- input[unique(c(all.vars(form),"age_band"))]
rm(prepared);gc(FALSE)
code <- c("R_cbh/primary/01_fit.R","R_cbh/primary/settings.R",settings$knots,
  "R_cbh/analysis/model.R","R_cbh/sensitivity/model.R")
if(settings$regional) code <- c(code,"R_cbh/covariates/regional.R","R_cbh/primary/specification.R","R_cbh/primary/00_prepare_regional.R")
code_hash <- vapply(code,cbh_file_hash,"")
data_hash <- cbh_file_hash(settings$data)
cbh_atomic_csv(data.frame(file=c(settings$data,code),md5=c(data_hash,unname(code_hash))),file.path(out,"fit_input_provenance.csv"))
cbh_atomic_csv(jobs,file.path(out,"jobs.csv"))
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
  d <- droplevels(input[input$age_band==job$age_band,,drop=FALSE])
  stopifnot(!anyNA(d),all(d$death %in% 0:1),all(d$calendar_year>=2000 & d$calendar_year<2025),
    max(abs(d$band_years-c(1,5,6,12,12,12,12)[match(job$age_band,ages)]/12))<1e-12)
  kt <- knot_table[knot_table$age_band==job$age_band,]
  knots <- setNames(lapply(c("pfpr_pct","calendar_year"),function(v) {
    z <- kt[kt$variable==v,];z$value[order(z$index)]
  }),c("pfpr_pct","calendar_year"))
  stopifnot(length(knots$pfpr_pct)==5L,length(knots$calendar_year)==6L)
  nuisance_hash <- cbh_hash(lapply(d[setdiff(names(d),c("pfpr_pct","age_band"))],as.vector))
  sig <- cbh_hash(list(job,knots,deparse(form),data_hash,code_hash,R.version.string,as.character(packageVersion("mgcv"))))
  path <- if(job$reuse) job$source_path else file.path(private,paste0(fit_id,".rds"))
  message("Starting ",fit_id,": ",nrow(d)," rows, ",sum(d$death)," deaths")
  saved <- if(file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
  if(!job$reuse && (is.null(saved) || !identical(saved$signature,sig))) {
    started <- as.character(Sys.time())
    saved <- run(form,d,knots)
    saved$started_at <- started;saved$finished_at <- as.character(Sys.time())
    saved$signature <- sig;saved$input_signature <- input_signature
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
      strict$signature <- ssig;strict$input_signature <- input_signature;strict$gamma <- 2;strict$knots <- knots
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
  stopifnot(diag$min_smoothing_hessian_eigenvalue>0,all(f$prior.weights==1),
    identical(f$family$link,"cloglog"),length(f$fitted.values)==nrow(d))
  for(v in names(f$model)) {
    want <- if(v=="offset(log(band_years))") log(d$band_years) else d[[v]]
    stopifnot(isTRUE(all.equal(f$model[[v]],want,check.attributes=FALSE)))
  }
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
    md5=cbh_file_hash(path),selected_version=version,source_path=job$source_path,source_md5=job$source_md5,prepared_data_md5=data_hash)
  cbh_atomic_csv(do.call(rbind,diags),file.path(out,"fit_diagnostics.csv"))
  cbh_atomic_csv(do.call(rbind,curves),file.path(out,"pfpr_curves.csv"))
  cbh_atomic_csv(do.call(rbind,summaries),file.path(out,"smooth_summaries.csv"))
  cbh_atomic_csv(do.call(rbind,selections),file.path(out,"fit_manifest.csv"))
  cbh_atomic_rds(components,file.path(private,"pfpr_components.rds"))
  writeLines(deparse(form),file.path(out,"model_formula.txt"))
  message("Completed ",fit_id,"; ",round(diag$elapsed_seconds,1)," fit seconds; Hessian minimum ",signif(diag$min_smoothing_hessian_eigenvalue,3))
  rm(saved,f,d);gc(FALSE)
}
stopifnot(length(diags)==7L)
cbh_atomic_csv(data.frame(file=code,md5=unname(code_hash)),file.path(out,"code_provenance.csv"))
stopifnot(identical(data_hash,cbh_file_hash(settings$data)))
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(out,"session_info.txt"))
message("Primary MAP gamma=2 fits complete")
