source("R_cbh/load_pipeline.R")
source("R_cbh/hiv/model.R")
source("R_cbh/analysis/model.R")
# Country withholding removes every child observation, but retains adolescence.
d<-data.frame(iso3=rep(c("AAA","BBB","CCC"),each=4),region=rep(c("West","West","East"),each=4),
  year=rep(2000:2003,3),adolescent_rate=rep(c(.2,.3,.4,.5),3),
  adolescent_censored=FALSE,child_rate=rep(c(.1,.2,.3,.4),3),child_censored=FALSE)
d$adolescent_rate[5]<-NA;d$adolescent_censored[5]<-TRUE
d$child_rate[2]<-NA;d$child_censored[2]<-TRUE
full<-hiv_stan_data(d)
held<-hiv_stan_data(d,"AAA")
stopifnot(full$M==12,held$M==8,!any(d$iso3[held$child_row]=="AAA"),
  identical(full$log_a,held$log_a),identical(full$time_basis,held$time_basis),
  full$NC_cens==1,held$NC_cens==0,full$NA_cens==1)
for(ss in list(full,held)) {
  ii<-which(ss$prev_c>0)
  stopifnot(all(d$iso3[ss$child_row[ii]]==d$iso3[ss$child_row[ss$prev_c[ii]]]))
}

# Exact country/entry-year join: no interview-year or prevalence substitution.
bands<-data.frame(country=c("NGA","NGA","STP","COM"),entry_year=c(2010,2011,2010,2010),
  interview_year=2020,hiv_prev_pct=c(9,9,9,9),log_hiv_prev=log(9))
panel<-data.frame(iso3=c("NGA","NGA","COM"),year=c(2010,2011,2010),
  hiv_incidence_per1000=c(.2,.3,.005),
  hiv_incidence_status=c("imputed_child_series","imputed_child_series","censored_posterior"))
j<-cbh_attach_incidence(bands,panel)
stopifnot(nrow(j)==4,identical(j$hiv_incidence_per1000,c(.2,.3,NA_real_,.005)),
  all(j$hiv_prev_pct==9),identical(j$log_hiv_incidence,log(j$hiv_incidence_per1000)),
  j$hiv_incidence_status[3]=="missing_country_year")
stopifnot(inherits(try(cbh_attach_incidence(bands,rbind(panel,panel[1,])),silent=TRUE),"try-error"))
spec<-cbh_trial_spec()
stopifnot("log_hiv_incidence" %in% spec$covariates,!"log_hiv_prev" %in% spec$covariates,
  !grepl("hiv_prev",paste(deparse(cbh_trial_formula(spec)),collapse=" ")),
  "log_hiv_prev" %in% cbh_trial_spec("prevalence")$covariates)

# Stable inverse normal CDF used by the censored innovation parameterization.
inv_log<-function(lp) {
  if(lp > -36) return(qnorm(exp(lp)))
  z<- -sqrt(-2*lp)
  for(k in 1:6) z<-z-(pnorm(z,log.p=TRUE)-lp)/exp(dnorm(z,log=TRUE)-pnorm(z,log.p=TRUE))
  z
}
lp<-c(-.001,-1,-20,-36,-100,-1000,-10000)
stopifnot(max(abs(vapply(lp,inv_log,numeric(1))-qnorm(lp,log.p=TRUE)))<1e-8)
cat("Incidence checks passed: country withholding, censoring, entry-year join, covariate replacement, and log-tail quantiles.\n")
