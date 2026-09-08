#!/usr/bin/env Rscript
source("R_cbh/load_pipeline.R")
source("R_cbh/hiv/model.R")
args<-commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% c("--cv","--force")))
private<-"data/derived_cbh/hiv_incidence"
out<-"results/cbh/hiv_incidence"
dir.create(private,recursive=TRUE,showWarnings=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_file<-"data/HIV_Epidemiology_Children_Adolescents_2025.xlsx"
files<-c(source_file,"R_cbh/hiv/incidence.stan","R_cbh/hiv/model.R","R_cbh/hiv/01_fit_incidence.R")
signature<-cbh_hash(list(vapply(files,cbh_file_hash,character(1)),as.character(packageVersion("rstan"))))
d<-hiv_incidence_data(source_file)
sd<-hiv_stan_data(d)
cbh_atomic_csv(d,file.path(out,"source_country_year.csv"))
message("Joint incidence model: ",sd$N," adolescent entries, ",sd$M," child entries, ",sd$C," countries.")
model<-rstan::stan_model(file="R_cbh/hiv/incidence.stan",auto_write=FALSE)
run_fit<-function(tag,held_out=character(),chains=4L,iter=2000L,seed=20260908L) {
  path<-file.path(private,paste0(tag,".rds"))
  ss<-hiv_stan_data(d,held_out)
  cache<-if(file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
  if(is.null(cache) || !identical(cache$signature,signature)) {
    f<-hiv_sample(model,ss,seed,chains,iter)
    cache<-list(fit=f,signature=signature,held_out=held_out,stan_data=ss,session_info=sessionInfo())
    cbh_atomic_rds(cache,path)
  }
  diag<-hiv_fit_diagnostics(cache$fit)
  cbh_atomic_csv(diag,file.path(out,paste0(tag,"_diagnostics.csv")))
  print(diag)
  if(diag$max_rhat>1.05 || diag$min_ess<100 || diag$divergences>0 || diag$treedepth_hits>0)
    stop("Numerical diagnostics need resolution for ",tag,"; fitted object retained.")
  cache
}
main<-run_fit("fit")
draws<-hiv_child_draws(main$fit,sd,d,n_draws=200)
panel<-d[c("iso3","country_name","region","year","child_rate","child_source_value","adolescent_source_value")]
panel$hiv_incidence_per1000<-exp(apply(draws,2,median))
panel$lower_95<-exp(apply(draws,2,quantile,.025))
panel$upper_95<-exp(apply(draws,2,quantile,.975))
panel$hiv_incidence_status<-ifelse(is.finite(d$child_rate),"reported_numeric",
  ifelse(d$child_censored,"censored_posterior","imputed_child_series"))
panel$imputation_signature<-signature
stopifnot(max(abs(panel$hiv_incidence_per1000[is.finite(d$child_rate)]-d$child_rate[is.finite(d$child_rate)]))<1e-8)
cbh_atomic_csv(panel,file.path(private,"child_incidence_country_year.csv"))
cbh_atomic_csv(panel,file.path(out,"child_incidence_country_year.csv"))
cbh_atomic_rds(list(log_incidence=draws,keys=d[c("iso3","year")],signature=signature),file.path(private,"child_incidence_draws.rds"))
sm<-rstan::summary(main$fit,pars=c("beta_between","beta_within","country_c_sd","sigma_c","rho_c"))$summary
cbh_atomic_csv(data.frame(parameter=rownames(sm),sm,row.names=NULL),file.path(out,"model_parameters.csv"))
if("--cv" %in% args) {
  # Five folds stratified by region; whole child series withheld, including censoring.
  countries<-unique(d[is.finite(d$child_rate) | d$child_censored,c("iso3","region")])
  set.seed(20260908)
  countries$fold<-0L
  for(r in unique(countries$region)) {
    ii<-which(countries$region==r)
    countries$fold[ii]<-sample(rep(1:5,length.out=length(ii)))
  }
  cbh_atomic_csv(countries,file.path(out,"validation_folds.csv"))
  cv<-list()
  for(k in 1:5) {
    message("Country-held-out validation fold ",k,"/5")
    held<-countries$iso3[countries$fold==k]
    f<-run_fit(paste0("cv",k),held,chains=2L,iter=2000L,seed=20260908L+k)
    pred<-hiv_child_draws(f$fit,f$stan_data,d,n_draws=400L,seed=20260908L+k)
    ii<-which(d$iso3 %in% held & (is.finite(d$child_rate)|d$child_censored))
    cv[[k]]<-data.frame(d[ii,c("iso3","region","year","child_rate","child_censored")],fold=k,
      predicted=exp(apply(pred[,ii,drop=FALSE],2,median)),
      lower_95=exp(apply(pred[,ii,drop=FALSE],2,quantile,.025)),
      upper_95=exp(apply(pred[,ii,drop=FALSE],2,quantile,.975)),
      probability_below_limit=colMeans(pred[,ii,drop=FALSE]<log(.01)))
    cbh_atomic_csv(do.call(rbind,cv),file.path(out,"country_held_out_predictions.csv"))
    rm(f,pred);gc(FALSE)
  }
}
message("Child incidence panel and posterior trajectories saved to ",private)
