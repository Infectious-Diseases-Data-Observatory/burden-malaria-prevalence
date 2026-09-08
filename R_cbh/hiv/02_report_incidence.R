#!/usr/bin/env Rscript
source("R_cbh/load_pipeline.R")
library(ggplot2)
out<-"results/cbh/hiv_incidence"
panel<-cbh_read_csv(file.path(out,"child_incidence_country_year.csv"))
source_data<-cbh_read_csv(file.path(out,"source_country_year.csv"))
targets<-panel[panel$iso3 %in% c("NGA","COM"),]
targets$country_name<-factor(targets$country_name,levels=c("Nigeria","Comoros"))
p<-ggplot(targets,aes(year,hiv_incidence_per1000))+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),fill="#176B87",alpha=.18)+
  geom_line(colour="#176B87",linewidth=.8)+facet_wrap(~country_name,scales="free_y")+
  scale_y_log10()+scale_x_continuous(breaks=c(2000,2005,2010,2015,2020,2024),
    expand=expansion(mult=.04))+labs(title="Imputed child HIV incidence",
    subtitle="Ages 0-14 | New infections per 1,000 uninfected population",
    x="Calendar year",y="Child HIV incidence rate (log scale)",
    caption="Median and pointwise 95% posterior predictive intervals. Country effects and annual dependence included.\nNeither country has an observed child series here; Comoros' adolescent predictor is censored below 0.01.\nIntervals condition on numeric source estimates; panel y-axis ranges differ.")+
  theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold"),
    plot.caption=element_text(hjust=0),panel.grid.minor=element_blank(),
    plot.margin=margin(12,22,12,12))
ggsave(file.path(out,"imputed_child_incidence_nga_com.png"),p,device=ragg::agg_png,width=11,height=6,dpi=180,bg="white")
cbh_atomic_csv(targets,file.path(out,"imputed_nga_com.csv"))
fit<-readRDS("data/derived_cbh/hiv_incidence/fit.rds")$fit
sm<-rstan::summary(fit)$summary
# Known numeric source rates are deterministic transformed quantities, not
# sampled parameters; their ESS/R-hat calculations are undefined or misleading.
diagnostic_rows<-sm[,"sd"]>1e-8 & is.finite(sm[,"Rhat"]) & is.finite(sm[,"n_eff"])
extra<-data.frame(max_rhat_all=max(sm[diagnostic_rows,"Rhat"]),
  min_ess_all=min(sm[diagnostic_rows,"n_eff"]))
cbh_atomic_csv(extra,file.path(out,"full_parameter_diagnostics.csv"))
print(extra)
if(file.exists(file.path(out,"country_held_out_predictions.csv"))) {
  cv<-cbh_read_csv(file.path(out,"country_held_out_predictions.csv"))
  stopifnot(setequal(unique(cv$fold),1:5))
  numeric<-cv[is.finite(cv$child_rate),]
  by_country<-do.call(rbind,lapply(split(numeric,numeric$iso3),function(z)
    data.frame(iso3=z$iso3[1],region=z$region[1],n=nrow(z),
      mean_squared_log_error=mean((log(z$predicted)-log(z$child_rate))^2),
      mean_log_error=mean(log(z$predicted)-log(z$child_rate)),
      coverage_95=mean(z$child_rate>=z$lower_95 & z$child_rate<=z$upper_95))))
  cbh_atomic_csv(by_country,file.path(out,"validation_by_country.csv"))
  score<-function(z) data.frame(countries=nrow(z),
    country_weighted_rmse_log=sqrt(mean(z$mean_squared_log_error)),
    country_weighted_bias_log=mean(z$mean_log_error),coverage_95=mean(z$coverage_95))
  metrics<-rbind(cbind(group="All countries",score(by_country)),
    do.call(rbind,lapply(split(by_country,by_country$region),function(z)cbind(group=z$region[1],score(z)))))
  cbh_atomic_csv(metrics,file.path(out,"validation_metrics.csv"))
  censored<-cv[cv$child_censored,]
  cbh_atomic_csv(data.frame(censored_observations=nrow(censored),
    mean_probability_below_reported_limit=mean(censored$probability_below_limit)),
    file.path(out,"validation_censored.csv"))
  latest<-numeric[numeric$year==max(numeric$year),]
  v<-ggplot(latest,aes(child_rate,predicted))+
    geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60")+
    geom_linerange(aes(ymin=lower_95,ymax=upper_95),colour="#176B87",alpha=.3)+
    geom_point(colour="#176B87",size=2)+scale_x_log10()+scale_y_log10()+
    labs(title="Predicting entirely withheld child incidence series",
      subtitle=paste("Five-fold country-held-out validation |",max(numeric$year),"snapshot"),
      x="Reported child incidence per 1,000",y="Predicted child incidence per 1,000",
      caption="Each country's entire child series was withheld; its adolescent series remained available.\nBars: pointwise 95% predictive intervals. Numeric reported child rates only; source uncertainty not included.")+
    theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0),
      panel.grid.minor=element_blank())
  ggsave(file.path(out,"country_held_out_validation.png"),v,device=ragg::agg_png,width=8,height=7,dpi=180,bg="white")
  print(metrics)
  fmt<-function(v)formatC(v,digits=3,format="f")
  rows<-apply(metrics,1,function(z)paste0("| ",z[["group"]]," | ",z[["countries"]],
    " | ",fmt(as.numeric(z[["country_weighted_rmse_log"]]))," | ",
    sprintf("%.1f%%",100*as.numeric(z[["coverage_95"]]))," |"))
  diag<-cbh_read_csv(file.path(out,"fit_diagnostics.csv"))
  writeLines(c("# Child HIV incidence imputation", "",
    "Fitted from the existing UNICEF / UNAIDS 2025 workbook using both-sex rates per 1,000 uninfected population. The joint model includes 3,468 adolescent country-years and 2,121 child country-years, across 139 countries; 85 countries contribute child observations. Numeric child rates are preserved, `<0.01` is treated as left censoring, and missing child trajectories are predicted.","",
    "The child equation separates between-country and within-country adolescent associations. Child and adolescent equations have separate calendar-year spline terms, regional/country effects and AR(1) residuals. See the [model and run instructions](../../../R_cbh/hiv/README.md).", "",
    sprintf("Main fit: four chains, 1,000 post-warmup draws each; maximum nonconstant-quantity R-hat %.4f, minimum effective sample size %.0f, %d divergences and %d maximum-treedepth hits.",extra$max_rhat_all,extra$min_ess_all,diag$divergences,diag$treedepth_hits),"",
    "## Countries without child series", "",
    "Nigeria's numeric adolescent trajectory and Comoros' censored adolescent series inform their predicted child rates. Liberia and São Tomé and Príncipe have no usable adolescent incidence series in this workbook and remain missing. Predictions include country-effect and temporally correlated residual uncertainty.","",
    "![Imputed child incidence](imputed_child_incidence_nga_com.png)","",
    "## Country-held-out validation", "",
    "Five folds withhold each country's entire child series, including censored child values; adolescent observations remain available. The table gives country-weighted log-scale error and pointwise interval coverage for numeric reported child rates. Countries with only censored child values cannot contribute to these numeric-error metrics; censoring predictions are reported separately.","",
    "| Group | Countries with numeric child rates | Log RMSE | 95% interval coverage |",
    "|---|---:|---:|---:|",rows,"",
    "![Country-held-out validation](country_held_out_validation.png)","",
    "Validation measures prediction of published estimated rates, not independently observed infections. Pointwise coverage below 95% indicates undercoverage in that group. The transport assumption for entirely unobserved countries remains an assumption, particularly for countries with only censored adolescent information.","",
    "## Downstream use and limits", "",
    "The active mortality input uses `log_hiv_incidence`, joined by country and band-entry year. It replaces prevalence; adolescent incidence does not enter the mortality formula. Two hundred coherent country-year trajectories are saved, and ten are used for exploratory uncertainty propagation through the age-band mortality model. The median-incidence fit is reported separately from the pooled refits.","",
    "Published source lower/upper estimates are retained but not included in the likelihood. Intervals condition on numeric source rates. Mortality propagation also conditions on fixed smoothing parameters and omits survey-design uncertainty; it is not joint outcome-compatible imputation. Ten imputations have Monte Carlo uncertainty, reported with the pooled contrasts.","",
    "- [Fitted incidence panel](child_incidence_country_year.csv)",
    "- [Validation predictions](country_held_out_predictions.csv), [metrics](validation_metrics.csv), [censored-rate validation](validation_censored.csv)",
    "- [Parameter estimates](model_parameters.csv), [main diagnostics](fit_diagnostics.csv), [all nonconstant-quantity diagnostics](full_parameter_diagnostics.csv)",
    "- [Incidence-adjusted mortality results](../age_band_hiv_incidence_shared_time_v3/REPORT.md)"),file.path(out,"REPORT.md"))
}
