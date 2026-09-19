#!/usr/bin/env Rscript
# Multiple imputation of the remaining survey-region covariate gaps (whole-survey
# gaps in published indicators and the pre-2002 wealth index) by chained equations
# at the survey-region level. Input is the current overlay after the UNICEF and
# within-survey region-mean fallbacks; observed values are never changed.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
suppressMessages(library(mice))
out_private <- "data/derived_cbh/regional_adjustment/imputed_v4"
out <- "results/cbh/covariate_imputation_v4"
for(p in c(out_private,out)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
cfg <- cbh_config()
spec <- cbh_primary_regional_spec(cbh_primary_settings("regional"))
regional <- spec$regional
overlay_path <- "data/derived_cbh/regional_adjustment/planned17_audit/regional_covariates_wide.csv"
wide <- cbh_read_csv(overlay_path); cbh_unique(wide,c("survey","regkey"),"Regional overlay")
reg <- cbh_read_csv(cfg$registry)
wide$country <- reg$iso3[match(wide$survey,reg$svkey)]
wide$survey_year <- reg$year[match(wide$survey,reg$svkey)]
stopifnot(!anyNA(wide$country),!anyNA(wide$survey_year))
# Auxiliary predictors: regional MAP PfPR at the survey year and the national series
# (already imputed for 2001 and 2024 by 15_impute_national.R).
map <- cbh_read_csv(cfg$annual_map)
wide$pfpr_survey_year <- map$pfpr2_10[match(paste(wide$survey,wide$regkey,wide$survey_year),paste(map$svkey,map$regkey,map$year))]
nat <- cbh_read_csv(file.path(out_private,"national_covariates_imputed.csv"))
k <- match(paste(wide$country,wide$survey_year),paste(nat$iso3,nat$year))
stopifnot(!anyNA(k))
wide$log_gdp_pc <- nat$log_gdp_pc[k]; wide$log_health_expenditure_pc <- nat$log_health_expenditure_pc[k]
wide$political_stability <- nat$political_stability[k]

## ---- imputation model -------------------------------------------------------------------
vars <- c(regional,"pfpr_survey_year")
X <- wide[c(vars,"survey_year","log_gdp_pc","log_health_expenditure_pc","political_stability","country")]
X$country <- factor(X$country)
before <- X[regional]
meth <- make.method(X); meth[vars] <- "pmm"; meth[setdiff(names(X),vars)] <- ""
pred <- make.predictorMatrix(X)
m <- 10L; maxit <- 20L; seed <- 20260918L
imp <- mice(X,m=m,maxit=maxit,method=meth,predictorMatrix=pred,seed=seed,printFlag=FALSE)
stopifnot(!any(is.na(complete(imp,1)[vars])))
completed <- lapply(seq_len(m),function(i) complete(imp,i)[regional])
point <- Reduce(`+`,completed)/m
# Percentages stay in [0,100] and wealth in [1,5] because pmm donates observed values;
# the mean of donors keeps the same bounds.
for(v in regional) stopifnot(all(point[[v]]>=min(before[[v]],na.rm=TRUE)-1e-9),all(point[[v]]<=max(before[[v]],na.rm=TRUE)+1e-9))

## ---- outputs -----------------------------------------------------------------------------
res <- wide[c("survey","country","survey_year","regkey")]
for(v in regional) {
  res[[v]] <- point[[v]]
  res[[paste0(v,"_before_model_imputation")]] <- before[[v]]
  res[[paste0(v,"_model_imputed")]] <- !is.finite(before[[v]])
}
stopifnot(all(vapply(regional,function(v) all(abs(res[[v]][!res[[paste0(v,"_model_imputed")]]]-before[[v]][!is.na(before[[v]])])<1e-12),logical(1))))
cbh_atomic_csv(res,file.path(out_private,"regional_covariates_wide.csv"))
cbh_atomic_rds(list(imputations=completed,keys=wide[c("survey","regkey")],imputed=as.data.frame(lapply(before,function(z)!is.finite(z))),
  m=m,maxit=maxit,seed=seed,method="pmm",predictors=names(X),loggedEvents=imp$loggedEvents),
  file.path(out_private,"regional_imputations.rds"))
by_var <- data.frame(variable=regional,survey_regions=nrow(res),
  imputed_regions=vapply(regional,function(v) sum(res[[paste0(v,"_model_imputed")]]),0L),
  imputed_surveys=vapply(regional,function(v) length(unique(res$survey[res[[paste0(v,"_model_imputed")]]])),0L),
  observed_mean=vapply(regional,function(v) mean(before[[v]],na.rm=TRUE),0),
  imputed_mean=vapply(regional,function(v) if(any(res[[paste0(v,"_model_imputed")]])) mean(point[[v]][res[[paste0(v,"_model_imputed")]]]) else NA_real_,0),
  between_imputation_sd=vapply(regional,function(v) { ii <- res[[paste0(v,"_model_imputed")]]; if(!any(ii)) return(NA_real_)
    mean(apply(sapply(completed,function(d) d[[v]][ii]),1,sd)) },0),row.names=NULL)
cbh_atomic_csv(by_var,file.path(out,"regional_imputation_by_variable.csv"))
by_survey <- do.call(rbind,lapply(regional,function(v) { ii <- res[[paste0(v,"_model_imputed")]]
  if(!any(ii)) return(NULL)
  t <- as.data.frame(table(survey=res$survey[ii]),stringsAsFactors=FALSE); names(t)[2] <- "regions_imputed"
  t$country <- res$country[match(t$survey,res$survey)]; t$variable <- v; t[c("survey","country","variable","regions_imputed")] }))
cbh_atomic_csv(by_survey,file.path(out,"regional_imputation_by_survey.csv"))
ragg::agg_png(file.path(out,"mice_convergence.png"),width=1800,height=2400,res=200)
print(plot(imp,layout=c(4,7))); dev.off()
ragg::agg_png(file.path(out,"mice_density.png"),width=2200,height=1600,res=200)
plotted <- regional[vapply(regional,function(v) sum(!is.finite(before[[v]]))>=2L,logical(1))]
print(densityplot(imp,as.formula(paste("~",paste(plotted,collapse="+"))))); dev.off()
writeLines(c("# Regional covariate imputation (imputed_v4)","",
  sprintf("Chained-equation multiple imputation (`mice`, predictive mean matching, m = %d, maxit = %d, seed %d) at the survey-region level on %d survey-regions from %d surveys. Targets: the 13 regional covariates and regional MAP PfPR at the survey year; predictors: all of these plus survey year, the national series at survey year (log GDP, log health expenditure, political stability; 2001/2024 gaps already filled by `15_impute_national.R`) and a country factor. Observed values are unchanged. The point overlay is the mean of the 10 imputations; all 10 are kept for propagation.",m,maxit,seed,nrow(res),length(unique(res$survey))),"",
  "| Variable | Survey-regions imputed | Surveys | Observed mean | Imputed mean | Between-imputation SD |","|---|---:|---:|---:|---:|---:|",
  sprintf("| %s | %d | %d | %.2f | %s | %s |",by_var$variable,by_var$imputed_regions,by_var$imputed_surveys,by_var$observed_mean,
    ifelse(is.na(by_var$imputed_mean),"—",sprintf("%.2f",by_var$imputed_mean)),ifelse(is.na(by_var$between_imputation_sd),"—",sprintf("%.2f",by_var$between_imputation_sd))),"",
  "Whole-survey gaps borrow from other surveys of the same country (country factor) and from the covariate relationships; predictive mean matching donates observed values, so imputations stay within observed ranges. Including PfPR as an auxiliary variable follows the usual multiple-imputation advice to include the analysis model's variables; the imputed covariates are adjustment variables, not the exposure or outcome. Imputation uncertainty is not reflected in the point-imputation fit; it is propagated in the separate multiple-imputation check.","",
  "[Convergence traces](mice_convergence.png) · [Observed versus imputed densities](mice_density.png) · [By variable](regional_imputation_by_variable.csv) · [By survey](regional_imputation_by_survey.csv)"),
  file.path(out,"REGIONAL_IMPUTATION.md"))
paths <- c(overlay_path,cfg$registry,cfg$annual_map,file.path(out_private,"national_covariates_imputed.csv"),
  "R_cbh/covariates/16_impute_regional.R","R_cbh/primary/specification.R")
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"regional_imputation_provenance.csv"))
print(by_var[,c("variable","imputed_regions","imputed_surveys","observed_mean","imputed_mean")]); message("Regional covariate imputation complete")
