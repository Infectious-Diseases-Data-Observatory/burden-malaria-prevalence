#!/usr/bin/env Rscript
# Recalculate primary contrasts and national burden, using saved national inputs.
source("R_cbh/load_pipeline.R")
library(mgcv)
source("R_cbh/primary/settings.R")
settings <- cbh_primary_settings()
root <- settings$out
out <- file.path(root,"burden")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
component_path <- file.path(settings$private,"pfpr_components.rds")
components <- readRDS(component_path)
manifest <- cbh_read_csv(file.path(root,"fit_manifest.csv"))
diag <- cbh_read_csv(file.path(root,"fit_diagnostics.csv"))
stopifnot(length(components)==7,nrow(manifest)==7,all(diag$converged),all(diag$input_verified))
for(piece in components) stopifnot(identical(piece$model_md5,manifest$md5[match(piece$fit_id,manifest$fit_id)]))
ages <- cbh_config()$age_bands$age_band
years <- settings$years
series <- "map_full"
inputs <- c(component_path,file.path(root,c("fit_manifest.csv","fit_diagnostics.csv")))
all_age <- all_totals <- list()
for(year in years) {
  olddir <- file.path("results/cbh/age_band_separate_v1",paste0("country_burden_",year))
  age_path <- file.path(olddir,paste0("country_age_attributable_",year,".csv"))
  compare_path <- file.path(olddir,paste0("model_vs_ihme_malaria_",year,".csv"))
  nat_path <- file.path(olddir,paste0("national_pfpr_",year,".csv"))
  inputs <- c(inputs,age_path,compare_path,nat_path)
  old <- cbh_read_csv(age_path);previous <- cbh_read_csv(compare_path);nat <- cbh_read_csv(nat_path)
  cbh_unique(old,c("iso3","age_band"),"Original age inputs")
  stopifnot(nrow(old)==45*7,nrow(previous)==45,all(old$exposure_year==year),
    identical(old$pfpr_pct,nat$pfpr_pct[match(old$iso3,nat$iso3)]))
  keep <- c("iso3","age_band","Location","pfpr_pct","ihme_deaths","implied_person_years","ihme_rate_per100000",
    "age_mapping","status","map_population_coverage_within_raster","map_coverage_below_95pct","population_weight_year")
  bands <- old[keep]
  for(s in series) {
    res <- do.call(rbind,lapply(ages,function(a) {
      z <- bands[bands$age_band==a,]
      piece <- components[[paste0(s,"_age_",match(a,ages))]]
      use <- is.finite(z$pfpr_pct)
      nd <- data.frame(pfpr_pct=z$pfpr_pct[use]);zero <- nd;zero$pfpr_pct <- 0
      L <- PredictMat(piece$smooth,zero)-PredictMat(piece$smooth,nd)
      variance <- rowSums((L%*%piece$covariance)*L)
      stopifnot(all(is.finite(variance)),min(variance)>=-1e-9)
      z$log_hr_zero_vs_current <- z$log_hr_se <- NA_real_
      z$log_hr_zero_vs_current[use] <- drop(L%*%piece$coef)
      z$log_hr_se[use] <- sqrt(pmax(variance,0))
      z$hr_zero_vs_current <- exp(z$log_hr_zero_vs_current)
      z$hr_lower_95 <- exp(z$log_hr_zero_vs_current-1.96*z$log_hr_se)
      z$hr_upper_95 <- exp(z$log_hr_zero_vs_current+1.96*z$log_hr_se)
      z$attributable_fraction <- 1-z$hr_zero_vs_current
      z$af_lower_95 <- 1-z$hr_upper_95;z$af_upper_95 <- 1-z$hr_lower_95
      for(unit in c("deaths","rate_per100000")) {
        base <- z[[paste0("ihme_",unit)]]
        z[[paste0("counterfactual_",unit)]] <- base*z$hr_zero_vs_current
        z[[paste0("attributable_",unit)]] <- base*z$attributable_fraction
        z[[paste0("attributable_",unit,"_lower_95")]] <- base*z$af_lower_95
        z[[paste0("attributable_",unit,"_upper_95")]] <- base*z$af_upper_95
        stopifnot(max(abs(base-z[[paste0("counterfactual_",unit)]]-z[[paste0("attributable_",unit)]]),na.rm=TRUE)<1e-7)
      }
      z$zero_below_observed_support <- piece$support[1]>0
      z$current_pfpr_outside_central95 <- z$pfpr_pct<piece$support[2] | z$pfpr_pct>piece$support[3]
      z$current_pfpr_outside_observed_support <- z$pfpr_pct<piece$support[1] | z$pfpr_pct>piece$support[4]
      z$negative_attributable_estimate <- z$attributable_fraction<0
      z$series <- s;z$year <- year;z$exposure_source <- "MAP"
      z$temporal_transport_from_pre2016_fit <- s!="map_full" && year>2015
      z
    }))
    all_age[[paste(s,year)]] <- res
    total <- do.call(rbind,lapply(split(res,res$iso3),function(z) {
      j <- match(z$iso3[1],previous$iso3)
      stopifnot(nrow(z)==7,abs(sum(z$ihme_deaths)-previous$ihme_under5_deaths[j])<1e-7)
      data.frame(iso3=z$iso3[1],country=z$Location[1],year=year,series=s,pfpr_pct=z$pfpr_pct[1],
        exposure_source="MAP",status=z$status[1],ihme_under5_deaths=sum(z$ihme_deaths),
        attributable_under5_deaths=sum(z$attributable_deaths),counterfactual_under5_deaths=sum(z$counterfactual_deaths),
        attributable_fraction=sum(z$attributable_deaths)/sum(z$ihme_deaths),
        ihme_malaria_deaths=previous$ihme_malaria_deaths[j],
        map_population_coverage=z$map_population_coverage_within_raster[1],
        any_current_pfpr_outside_observed_support=any(z$current_pfpr_outside_observed_support),
        any_zero_below_observed_support=any(z$zero_below_observed_support),
        temporal_transport_from_pre2016_fit=any(z$temporal_transport_from_pre2016_fit))
    }))
    all_totals[[paste(s,year)]] <- total
  }
}
r <- do.call(rbind,all_age);totals <- do.call(rbind,all_totals)
cbh_unique(r,c("series","year","iso3","age_band"),"New country age estimates")
cbh_unique(totals,c("series","year","iso3"),"New country totals")
stopifnot(nrow(r)==3*45*7,nrow(totals)==3*45,
  all(table(totals$series[is.finite(totals$attributable_under5_deaths)],totals$year[is.finite(totals$attributable_under5_deaths)])==42))
cbh_atomic_csv(r,file.path(out,"country_age_estimates.csv"))
cbh_atomic_csv(totals,file.path(out,"country_totals.csv"))
# Year totals use the same estimable countries as the cause-specific comparator.
sums <- do.call(rbind,lapply(years,function(y) {
  z <- totals[totals$year==y & is.finite(totals$attributable_under5_deaths),]
  data.frame(year=y,countries=nrow(z),attributable_under5_deaths=sum(z$attributable_under5_deaths),
    counterfactual_under5_deaths=sum(z$counterfactual_under5_deaths),
    ihme_under5_deaths=sum(z$ihme_under5_deaths),ihme_malaria_deaths=sum(z$ihme_malaria_deaths))
}))
cbh_atomic_csv(sums,file.path(out,"year_summary.csv"))
# Preserve the equal person-time assumption for the unresolved IHME 2-4 bin.
for(y in years) for(iso in unique(r$iso3)) {
  z <- r[r$year==y & r$iso3==iso & r$age_band %in% ages[5:7],]
  stopifnot(nrow(z)==3,diff(range(z$ihme_deaths))<1e-7,
    diff(range(z$implied_person_years))<1e-7,diff(range(z$ihme_rate_per100000))<1e-7)
}
# Curve summaries retain the covariance of both evaluation points.
curve_path <- file.path(root,"pfpr_curves.csv")
sm_path <- file.path(root,"smooth_summaries.csv")
d <- cbh_read_csv(curve_path);sm <- cbh_read_csv(sm_path)
cbh_atomic_csv(sm[sm$term=="s(pfpr_pct)",c("series","age_band","edf")],file.path(root,"pfpr_edf.csv"))
hr <- d[d$pfpr_pct==40,c("series","age_band","log_hazard_ratio","standard_error","pfpr_p025","pfpr_p975")]
hr$hazard_ratio_40_to_20 <- exp(-hr$log_hazard_ratio)
hr$lower_95 <- exp(-hr$log_hazard_ratio-1.96*hr$standard_error)
hr$upper_95 <- exp(-hr$log_hazard_ratio+1.96*hr$standard_error)
hr$within_central_support <- hr$pfpr_p025<=20 & hr$pfpr_p975>=40
cbh_atomic_csv(hr,file.path(root,"pfpr_40_to_20_contrasts.csv"))
zero <- d[d$pfpr_pct==0,c("series","age_band","log_hazard_ratio","standard_error","pfpr_min")]
zero$hazard_ratio_20_to_zero <- exp(zero$log_hazard_ratio)
zero$attributable_fraction_at_20pct <- 1-zero$hazard_ratio_20_to_zero
zero$lower_95 <- 1-exp(zero$log_hazard_ratio+1.96*zero$standard_error)
zero$upper_95 <- 1-exp(zero$log_hazard_ratio-1.96*zero$standard_error)
zero$zero_below_observed_support <- zero$pfpr_min>0
cbh_atomic_csv(zero,file.path(root,"pfpr_20_to_zero_contrasts.csv"))
anchors <- do.call(rbind,lapply(components,function(piece) {
  nd <- data.frame(pfpr_pct=c(10,20,30,50));zero <- nd;zero$pfpr_pct <- 0
  L <- PredictMat(piece$smooth,zero)-PredictMat(piece$smooth,nd)
  estimate <- drop(L%*%piece$coef);variance <- rowSums((L%*%piece$covariance)*L)
  stopifnot(all(is.finite(variance)),min(variance)>=-1e-9)
  se <- sqrt(pmax(variance,0))
  data.frame(age_band=piece$age_band,pfpr_pct=nd$pfpr_pct,reference_pct=0,
    attributable_fraction=1-exp(estimate),lower_95=1-exp(estimate+1.96*se),upper_95=1-exp(estimate-1.96*se),
    zero_below_observed_support=piece$support[1]>0,
    current_outside_observed_support=nd$pfpr_pct<piece$support[1] | nd$pfpr_pct>piece$support[4])
}))
cbh_atomic_csv(anchors,file.path(root,"pfpr_attributable_fraction_anchors.csv"))
# Independent numerical comparison with the previous gamma=2 results.
old_path <- file.path(settings$reference,"burden/country_totals.csv")
old <- cbh_read_csv(old_path);old <- old[old$series=="map_full",]
j <- match(cbh_key(totals,c("year","iso3")),cbh_key(old,c("year","iso3")))
stopifnot(!anyNA(j),identical(is.na(totals$attributable_under5_deaths),is.na(old$attributable_under5_deaths[j])))
comparison <- totals[c("year","iso3","country")]
comparison$previous_deaths <- old$attributable_under5_deaths[j]
comparison$rerun_deaths <- totals$attributable_under5_deaths
comparison$difference_deaths <- comparison$rerun_deaths-comparison$previous_deaths
cbh_atomic_csv(comparison,file.path(out,"comparison_with_previous_primary.csv"))
old_curve_path <- file.path(settings$reference,"pfpr_curves.csv")
old <- cbh_read_csv(old_curve_path);old <- old[old$series=="map_full",]
j <- match(cbh_key(d,c("age_band","pfpr_pct")),cbh_key(old,c("age_band","pfpr_pct")))
stopifnot(!anyNA(j))
curve_comparison <- d[c("age_band","pfpr_pct")]
curve_comparison$log_hr_difference <- d$log_hazard_ratio-old$log_hazard_ratio[j]
curve_comparison$se_difference <- d$standard_error-old$standard_error[j]
cbh_atomic_csv(curve_comparison,file.path(root,"comparison_with_previous_curves.csv"))
# Synthetic-cohort life table for DRC. Annual death counts remain separate.
life <- do.call(rbind,lapply(years,function(y) {
  z <- r[r$year==y & r$iso3=="COD",];z <- z[match(ages,z$age_band),]
  ip <- file.path("results/cbh/age_band_separate_v1",paste0("country_burden_",y),"ihme_disjoint_age_inputs.csv")
  inputs <<- c(inputs,ip)
  src <- cbh_read_csv(ip);src <- src[src$iso3=="COD",]
  edges <- c(0,28/365.25*12,6,12,24,36,48,60)
  H <- z$ihme_rate_per100000/1e5*diff(edges)/12
  early <- src$rate_per100000[src$source_age=="0 to 6 days (early neonatal)"]
  late <- src$rate_per100000[src$source_age=="7 to 27 days (late neonatal)"]
  stopifnot(length(early)==1,length(late)==1)
  H[1] <- early/1e5*7/365.25+late/1e5*21/365.25
  H0 <- H*z$hr_zero_vs_current
  q <- -expm1(-H);q0 <- -expm1(-H0)
  S <- exp(-cumsum(H));S0 <- exp(-cumsum(H0))
  stopifnot(max(abs(S-cumprod(1-q)))<1e-12,max(abs(S0-cumprod(1-q0)))<1e-12)
  data.frame(year=y,iso3="COD",age_band=ages,start_month=edges[-8],end_month=edges[-1],
    observed_q=q,zero_pfpr_q=q0,observed_survival=S,zero_pfpr_survival=S0)
}))
cbh_atomic_csv(life,file.path(out,"drc_life_table.csv"))
inputs <- unique(c(inputs,curve_path,sm_path,old_path,old_curve_path,"R_cbh/primary/02_effects.R","R_cbh/primary/settings.R"))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"provenance.csv"))
message("Primary effects and country burden complete; maximum country change = ",
  signif(max(abs(comparison$difference_deaths),na.rm=TRUE),4)," deaths")
