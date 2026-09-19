#!/usr/bin/env Rscript
# Child HIV incidence imputation model.
#   default      : the frozen panel used by the current primary (refit only with --force)
#   --extended   : adds countries with no UNAIDS adolescent series (latent adolescent
#                  path) and writes *_extended outputs for the imputed-covariate version
#   --cv         : five-fold country-held-out validation of child series
#   --validate-latent : (with --extended) treat two observed West and Central Africa
#                  countries' adolescent series as latent and withhold their child
#                  series, to check imputation for countries without adolescent data
source("R_cbh/load_pipeline.R")
source("R_cbh/hiv/model.R")
args<-commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% c("--cv","--force","--extended","--validate-latent")))
extended<-"--extended" %in% args
private<-"data/derived_cbh/hiv_incidence"
out<-"results/cbh/hiv_incidence"
dir.create(private,recursive=TRUE,showWarnings=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_file<-"data/HIV_Epidemiology_Children_Adolescents_2025.xlsx"
files<-c(source_file,"R_cbh/hiv/incidence.stan","R_cbh/hiv/model.R","R_cbh/hiv/01_fit_incidence.R")
# Countries in the DHS sample with no UNAIDS incidence series. Liberia's UNICEF region
# comes from its other workbook rows; Sao Tome and Principe has no rows and is assigned
# to West and Central Africa (UNICEF regional classification).
extra<-if(extended) c(LBR="West and Central Africa",STP="West and Central Africa") else NULL
suffix<-if(extended) "_extended" else ""
signature<-cbh_hash(list(vapply(files,cbh_file_hash,character(1)),as.character(packageVersion("rstan")),extra))
d<-hiv_incidence_data(source_file,extra)
sd<-hiv_stan_data(d)
cbh_atomic_csv(d,file.path(out,paste0("source_country_year",suffix,".csv")))
message("Joint incidence model: ",sd$N," adolescent entries (",sd$NA_miss," latent), ",sd$M," child entries, ",sd$C," countries.")
model<-rstan::stan_model(file="R_cbh/hiv/incidence.stan",auto_write=FALSE)
run_fit<-function(tag,held_out=character(),latent=character(),chains=4L,iter=2000L,seed=20260908L) {
  path<-file.path(private,paste0(tag,".rds"))
  ss<-hiv_stan_data(d,held_out,latent)
  cache<-if(file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
  if(is.null(cache) || !identical(cache$signature,signature)) {
    if(!extended && tag=="fit" && !is.null(cache))
      stop("The base HIV panel is frozen for the current primary version; its cached fit predates the ",
           "current model code. Use --extended for the imputed-covariate panel, or --force to refit the base panel.")
    f<-hiv_sample(model,ss,seed,chains,iter)
    cache<-list(fit=f,signature=signature,held_out=held_out,latent_adolescent=latent,stan_data=ss,session_info=sessionInfo())
    cbh_atomic_rds(cache,path)
  }
  diag<-hiv_fit_diagnostics(cache$fit)
  cbh_atomic_csv(diag,file.path(out,paste0(tag,"_diagnostics.csv")))
  print(diag)
  if(diag$max_rhat>1.05 || diag$min_ess<100 || diag$divergences>0 || diag$treedepth_hits>0)
    stop("Numerical diagnostics need resolution for ",tag,"; fitted object retained.")
  cache
}
main<-run_fit(paste0("fit",suffix))
draws<-hiv_child_draws(main$fit,sd,d,n_draws=200)
panel<-d[c("iso3","country_name","region","year","child_rate","child_source_value","adolescent_source_value")]
panel$hiv_incidence_per1000<-exp(apply(draws,2,median))
panel$lower_95<-exp(apply(draws,2,quantile,.025))
panel$upper_95<-exp(apply(draws,2,quantile,.975))
panel$hiv_incidence_status<-ifelse(is.finite(d$child_rate),"reported_numeric",
  ifelse(d$child_censored,"censored_posterior",
    ifelse(d$adolescent_missing,"imputed_no_adolescent_series","imputed_child_series")))
panel$imputation_signature<-signature
stopifnot(max(abs(panel$hiv_incidence_per1000[is.finite(d$child_rate)]-d$child_rate[is.finite(d$child_rate)]))<1e-8)
cbh_atomic_csv(panel,file.path(private,paste0("child_incidence_country_year",suffix,".csv")))
cbh_atomic_csv(panel,file.path(out,paste0("child_incidence_country_year",suffix,".csv")))
cbh_atomic_rds(list(log_incidence=draws,keys=d[c("iso3","year")],signature=signature),
  file.path(private,paste0("child_incidence_draws",suffix,".rds")))
sm<-rstan::summary(main$fit,pars=c("beta_between","beta_within","country_c_sd","sigma_c","rho_c"))$summary
cbh_atomic_csv(data.frame(parameter=rownames(sm),sm,row.names=NULL),file.path(out,paste0("model_parameters",suffix,".csv")))
if(extended) {
  # Compare with the frozen base panel: shared country-years should agree within MCMC noise.
  base_path<-file.path(private,"child_incidence_country_year.csv")
  if(file.exists(base_path)) {
    base<-cbh_read_csv(base_path)
    j<-match(paste(panel$iso3,panel$year),paste(base$iso3,base$year))
    shared<-!is.na(j)
    cmp<-data.frame(iso3=panel$iso3[shared],year=panel$year[shared],status=panel$hiv_incidence_status[shared],
      base=base$hiv_incidence_per1000[j[shared]],extended=panel$hiv_incidence_per1000[shared])
    cmp$log_ratio<-log(cmp$extended/cmp$base)
    cbh_atomic_csv(cmp,file.path(out,"extended_vs_base_panel.csv"))
    imp<-cmp[cmp$status!="reported_numeric",]
    message(sprintf("Extended vs base panel: %d imputed/censored country-years; median |log ratio| %.3f, max %.3f",
      nrow(imp),median(abs(imp$log_ratio)),max(abs(imp$log_ratio))))
  }
  new<-panel[panel$hiv_incidence_status=="imputed_no_adolescent_series",]
  cbh_atomic_csv(new,file.path(out,"imputed_no_adolescent_series.csv"))
}
if("--validate-latent" %in% args) {
  stopifnot(extended)
  # Two observed West and Central Africa countries with complete child series.
  wca<-unique(d$iso3[d$region=="West and Central Africa" & !d$adolescent_missing])
  full<-vapply(wca,function(i) all(is.finite(d$child_rate[d$iso3==i])),logical(1))
  latent<-sort(wca[full])[1:2]
  message("Latent-adolescent validation countries: ",paste(latent,collapse=", "))
  f<-run_fit("fit_extended_latent_validation",held_out=latent,latent=latent,chains=2L,iter=2000L,seed=20260918L)
  pred<-hiv_child_draws(f$fit,f$stan_data,d,n_draws=400L,seed=20260918L)
  ii<-which(d$iso3 %in% latent)
  v<-data.frame(d[ii,c("iso3","region","year","child_rate")],
    predicted=exp(apply(pred[,ii,drop=FALSE],2,median)),
    lower_95=exp(apply(pred[,ii,drop=FALSE],2,quantile,.025)),
    upper_95=exp(apply(pred[,ii,drop=FALSE],2,quantile,.975)))
  v$covered<-v$child_rate>=v$lower_95 & v$child_rate<=v$upper_95
  v$log_error<-log(v$predicted/v$child_rate)
  cbh_atomic_csv(v,file.path(out,"latent_adolescent_validation.csv"))
  message(sprintf("Latent-adolescent validation: coverage %.0f%%, median |log error| %.2f",
    100*mean(v$covered),median(abs(v$log_error))))
}
if("--cv" %in% args) {
  # Five folds stratified by region; whole child series withheld, including censoring.
  countries<-unique(d[is.finite(d$child_rate) | d$child_censored,c("iso3","region")])
  set.seed(20260908)
  countries$fold<-0L
  for(r in unique(countries$region)) {
    ii<-which(countries$region==r)
    countries$fold[ii]<-sample(rep(1:5,length.out=length(ii)))
  }
  cbh_atomic_csv(countries,file.path(out,paste0("validation_folds",suffix,".csv")))
  cv<-list()
  for(k in 1:5) {
    message("Country-held-out validation fold ",k,"/5")
    held<-countries$iso3[countries$fold==k]
    f<-run_fit(paste0("cv",k,suffix),held,chains=2L,iter=2000L,seed=20260908L+k)
    pred<-hiv_child_draws(f$fit,f$stan_data,d,n_draws=400L,seed=20260908L+k)
    ii<-which(d$iso3 %in% held & (is.finite(d$child_rate)|d$child_censored))
    cv[[k]]<-data.frame(d[ii,c("iso3","region","year","child_rate","child_censored")],fold=k,
      predicted=exp(apply(pred[,ii,drop=FALSE],2,median)),
      lower_95=exp(apply(pred[,ii,drop=FALSE],2,quantile,.025)),
      upper_95=exp(apply(pred[,ii,drop=FALSE],2,quantile,.975)),
      probability_below_limit=colMeans(pred[,ii,drop=FALSE]<log(.01)))
    cbh_atomic_csv(do.call(rbind,cv),file.path(out,paste0("country_held_out_predictions",suffix,".csv")))
    rm(f,pred);gc(FALSE)
  }
}
message("Child incidence panel and posterior trajectories saved to ",private)
