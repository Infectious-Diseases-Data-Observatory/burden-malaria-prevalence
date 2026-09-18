#!/usr/bin/env Rscript
# Secondary subnational burden: Nigerian states, 2024. Applies the saved primary
# age-band PfPR effects to IHME state all-cause deaths by age band, exactly as
# primary/02_effects.R does for countries. No refitting, no downloads.
# Run from the project root after the primary fits and effects stages.
source("R_cbh/load_pipeline.R")
library(mgcv)
source("R_cbh/primary/settings.R")
settings <- cbh_primary_settings(Sys.getenv("CBH_PRIMARY_VERSION","regional"))
root <- settings$out
out <- file.path(root,"nigeria_states")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
YEAR <- 2024L
ages <- cbh_config()$age_bands$age_band

## ---- inputs and their verification -------------------------------------------
ext <- "data/external/nigeria_states"
allcause_path <- list.files(ext,"^Data Explorer.*All causes.*Deaths.*[.]csv$",full.names=TRUE)
stopifnot(length(allcause_path)==1)
malaria_path <- file.path(ext,"GBD NG Malaria mortality incidence prevalence rates 2011-2024.csv")
pfpr_path <- "data/pfpr_admin1_ng_cd_2024.csv"
component_path <- file.path(settings$private,"pfpr_components.rds")
national_path <- file.path(root,"burden/country_age_estimates.csv")
inputs <- c(component_path,file.path(root,c("fit_manifest.csv","fit_diagnostics.csv")),
  allcause_path,malaria_path,pfpr_path,file.path(ext,"source.json"),national_path,
  "R_cbh/burden/05_nigeria_state_burden.R","R_cbh/primary/settings.R")
stopifnot(all(file.exists(inputs)))
components <- readRDS(component_path)
manifest <- cbh_read_csv(file.path(root,"fit_manifest.csv"))
diag <- cbh_read_csv(file.path(root,"fit_diagnostics.csv"))
stopifnot(length(components)==7,nrow(manifest)==7,all(diag$converged),all(diag$input_verified),all(diag$gamma==2))
for(piece in components) stopifnot(identical(piece$model_md5,manifest$md5[match(piece$fit_id,manifest$fit_id)]))

# State names differ across sources: IHME writes "FCT (Abuja)" and "Nasarawa",
# the MAP admin-1 extraction "Abuja" and "Nassarawa".
state_key <- function(x) {
  k <- gsub("[^a-z0-9]","",tolower(iconv(as.character(x),"","ASCII//TRANSLIT")))
  k[k %in% c("fctabuja","abuja")] <- "fct"; k[k=="nassarawa"] <- "nasarawa"; k
}

## ---- IHME all-cause deaths and implied person-years by state and age group ----
g <- cbh_read_csv(allcause_path)
g <- g[g$Year==YEAR & g$Location!="Global" & g$Measure=="Deaths" & g$Sex=="Both",]
stopifnot(all(g$Condition=="All causes"),all(is.finite(g$Value)),all(g$Value>=0))
num <- g[g$Unit=="Number",];rate <- g[g$Unit=="Rate (per 100,000)",]
j <- match(cbh_key(num,c("Location","Age")),cbh_key(rate,c("Location","Age")))
stopifnot(!anyNA(j))
src <- data.frame(state=num$Location,key=state_key(num$Location),source_age=num$Age,
  deaths=num$Value,deaths_lower=num$Lower,deaths_upper=num$Upper,rate_per100000=rate$Value[j])
src$implied_person_years <- src$deaths/src$rate_per100000*1e5
cbh_unique(src,c("key","source_age"),"IHME state age inputs")
states <- sort(unique(src$key))
stopifnot(length(states)==37)
pick <- function(age,col) {z <- src[src$source_age==age,];stopifnot(nrow(z)==37);setNames(z[[col]],z$key)[states]}
# GBD age group -> model band, the documented national convention: 0-27 days to
# the <1 completed-month effect; 1-5, 6-11 direct; 12-23 as "1 to 4" minus
# "2 to 4"; the "2 to 4" group shares one rate and divides deaths/person-time
# equally across 24-35, 36-47 and 48-59 months.
build <- function(col) {
  en <- pick("0 to 6 days (early neonatal)",col);ln <- pick("7 to 27 days (late neonatal)",col)
  y14 <- pick("1 to 4",col);y24 <- pick("2 to 4",col)
  cbind(`<1`=en+ln,`1-5`=pick("1 to 5 months",col),`6-11`=pick("6-11 months",col),
    `12-23`=y14-y24,`24-35`=y24/3,`36-47`=y24/3,`48-59`=y24/3)
}
band_deaths <- build("deaths");band_py <- build("implied_person_years")
mapping <- c(`<1`="0-27 days mapped to model <1 completed month",`1-5`="direct age match",
  `6-11`="direct age match",`12-23`="IHME 1 to 4 minus 2 to 4",
  `24-35`="IHME 2-4 rate shared; deaths and person-time divided equally",
  `36-47`="IHME 2-4 rate shared; deaths and person-time divided equally",
  `48-59`="IHME 2-4 rate shared; deaths and person-time divided equally")
stopifnot(identical(colnames(band_deaths),ages))
u5_deaths <- pick("Under 5","deaths");u5_py <- pick("Under 5","implied_person_years")
stopifnot(max(abs(rowSums(band_deaths)-u5_deaths)/u5_deaths)<1e-9,all(band_deaths>=0),all(band_py>0))

# The state export must be the same IHME vintage as the national input used by
# 02_effects.R: state bands summed over Nigeria must equal the national bands.
nat <- cbh_read_csv(national_path);nat <- nat[nat$iso3=="NGA" & nat$year==YEAR,]
stopifnot(nrow(nat)==7)
nat_band <- setNames(nat$ihme_deaths,nat$age_band)[ages]
stopifnot(max(abs(colSums(band_deaths)-nat_band)/nat_band)<1e-6)

## ---- MAP PfPR2-10 by state, 2024 ----------------------------------------------
pf <- cbh_read_csv(pfpr_path);pf <- pf[pf$iso3=="NGA",];pf$key <- state_key(pf$area)
stopifnot(setequal(pf$key,states),!anyDuplicated(pf$key),all(is.finite(pf$pfpr_pct)))
pfpr <- setNames(pf$pfpr_pct,pf$key)[states]
state_name <- setNames(src$state[match(states,src$key)],states)

## ---- apply the saved age-band effects -----------------------------------------
res <- do.call(rbind,lapply(seq_along(ages),function(i) {
  a <- ages[i];piece <- components[[paste0("map_full_age_",i)]]
  stopifnot(identical(piece$age_band,a))
  nd <- data.frame(pfpr_pct=unname(pfpr));zero <- nd;zero$pfpr_pct <- 0
  L <- PredictMat(piece$smooth,zero)-PredictMat(piece$smooth,nd)
  variance <- rowSums((L%*%piece$covariance)*L)
  stopifnot(all(is.finite(variance)),min(variance)>=-1e-9)
  z <- data.frame(iso3="NGA",state=unname(state_name),key=states,age_band=a,year=YEAR,series="map_full",
    exposure_source="MAP",pfpr_pct=unname(pfpr),
    ihme_deaths=unname(band_deaths[,a]),implied_person_years=unname(band_py[,a]),
    age_mapping=unname(mapping[a]))
  z$ihme_rate_per100000 <- z$ihme_deaths/z$implied_person_years*1e5
  z$log_hr_zero_vs_current <- drop(L%*%piece$coef)
  z$log_hr_se <- sqrt(pmax(variance,0))
  z$hr_zero_vs_current <- exp(z$log_hr_zero_vs_current)
  z$hr_lower_95 <- exp(z$log_hr_zero_vs_current-1.96*z$log_hr_se)
  z$hr_upper_95 <- exp(z$log_hr_zero_vs_current+1.96*z$log_hr_se)
  z$attributable_fraction <- 1-z$hr_zero_vs_current
  z$af_lower_95 <- 1-z$hr_upper_95;z$af_upper_95 <- 1-z$hr_lower_95
  z$counterfactual_deaths <- z$ihme_deaths*z$hr_zero_vs_current
  z$attributable_deaths <- z$ihme_deaths*z$attributable_fraction
  z$attributable_deaths_lower_95 <- z$ihme_deaths*z$af_lower_95
  z$attributable_deaths_upper_95 <- z$ihme_deaths*z$af_upper_95
  stopifnot(max(abs(z$ihme_deaths-z$counterfactual_deaths-z$attributable_deaths))<1e-7)
  z$zero_below_observed_support <- piece$support[1]>0
  z$current_pfpr_outside_central95 <- z$pfpr_pct<piece$support[2] | z$pfpr_pct>piece$support[3]
  z$current_pfpr_outside_observed_support <- z$pfpr_pct<piece$support[1] | z$pfpr_pct>piece$support[4]
  z$negative_attributable_estimate <- z$attributable_fraction<0
  z
}))
cbh_unique(res,c("key","age_band"),"State age estimates")
stopifnot(nrow(res)==37*7,!any(res$current_pfpr_outside_observed_support))
# Preserve the equal person-time assumption for the unresolved IHME 2-4 bin.
for(k in states) {
  z <- res[res$key==k & res$age_band %in% ages[5:7],]
  stopifnot(nrow(z)==3,diff(range(z$ihme_deaths))<1e-7,diff(range(z$implied_person_years))<1e-7)
}
cbh_atomic_csv(res[order(res$state,match(res$age_band,ages)),],file.path(out,"state_age_estimates_2024.csv"))

## ---- state totals and the cause-specific comparator ----------------------------
m <- cbh_read_csv(malaria_path)
m <- m[m$Year==YEAR & m$Measure=="Deaths" & m$Unit=="Rate (per 100,000)" & m$Location!="Nigeria",]
stopifnot(all(m$Condition=="Malaria"),all(m$Age=="Under 5"),all(m$Sex=="Both"))
m$key <- state_key(m$Location);stopifnot(setequal(m$key,states),!anyDuplicated(m$key))
totals <- do.call(rbind,lapply(states,function(k) {
  z <- res[res$key==k,];stopifnot(nrow(z)==7)
  mr <- m[m$key==k,]
  data.frame(iso3="NGA",state=state_name[[k]],key=k,year=YEAR,series="map_full",exposure_source="MAP",
    pfpr_pct=z$pfpr_pct[1],under5_person_years=u5_py[[k]],ihme_under5_deaths=u5_deaths[[k]],
    attributable_under5_deaths=sum(z$attributable_deaths),
    counterfactual_under5_deaths=sum(z$counterfactual_deaths),
    attributable_fraction=sum(z$attributable_deaths)/sum(z$ihme_deaths),
    # Person-years for the comparator come from the all-cause count/rate pair,
    # as in the annual national comparison; the IHME malaria rate is applied to it.
    ihme_malaria_rate_per100000=mr$Value,ihme_malaria_rate_lower=mr$Lower,ihme_malaria_rate_upper=mr$Upper,
    ihme_malaria_deaths=mr$Value/1e5*u5_py[[k]],
    any_current_pfpr_outside_central95=any(z$current_pfpr_outside_central95),
    any_zero_below_observed_support=any(z$zero_below_observed_support),
    any_negative_attributable_band=any(z$negative_attributable_estimate),row.names=NULL)
}))
totals$model_over_ihme_malaria <- totals$attributable_under5_deaths/totals$ihme_malaria_deaths
totals$attributable_rate_per100000 <- totals$attributable_under5_deaths/totals$under5_person_years*1e5
totals$ihme_malaria_share_of_allcause <- totals$ihme_malaria_deaths/totals$ihme_under5_deaths
cbh_unique(totals,"key","State totals")
stopifnot(nrow(totals)==37,all(is.finite(totals$attributable_under5_deaths)),all(totals$ihme_malaria_deaths>0))
cbh_atomic_csv(totals[order(-totals$attributable_under5_deaths),],file.path(out,"state_totals_2024.csv"))

## ---- reconciliation with the national estimate --------------------------------
# The national calculation evaluates each nonlinear curve at the population-
# weighted national mean PfPR (24.70%). Summing state estimates evaluated at
# state prevalence gives a different total; the difference is recorded, not hidden.
national_total <- sum(nat$attributable_deaths)
state_sum <- sum(totals$attributable_under5_deaths)
pop_weighted_pfpr <- sum(totals$pfpr_pct*totals$under5_person_years)/sum(totals$under5_person_years)
recon <- data.frame(year=YEAR,
  national_ihme_under5_deaths=sum(nat_band),state_sum_ihme_under5_deaths=sum(totals$ihme_under5_deaths),
  national_pfpr_pct=nat$pfpr_pct[1],state_person_year_weighted_pfpr_pct=pop_weighted_pfpr,
  national_attributable_deaths=national_total,state_sum_attributable_deaths=state_sum,
  difference_deaths=state_sum-national_total,relative_difference=state_sum/national_total-1,
  national_attributable_fraction=national_total/sum(nat_band),state_sum_attributable_fraction=state_sum/sum(totals$ihme_under5_deaths),
  ihme_malaria_deaths_state_sum=sum(totals$ihme_malaria_deaths),
  model_over_ihme_malaria_national=state_sum/sum(totals$ihme_malaria_deaths),
  note="National PfPR uses GPW 2020 density x cell area with national admin-0 overlap; the state extraction (archive/20_pfpr_admin1_ng_cd.R) used MAP admin-1 boundaries and density-only GPW weights.")
cbh_atomic_csv(recon,file.path(out,"national_reconciliation_2024.csv"))
stopifnot(abs(recon$relative_difference)<0.05)

writeLines(c(paste0("Secondary subnational burden, Nigerian states, ",YEAR),
  "Model: seven separate primary MAP gamma=2 age-band fits (series map_full), verified against fit_manifest.csv.",
  "Exposure: population-weighted MAP PfPR2-10 by state for 2024 from data/pfpr_admin1_ng_cd_2024.csv (inherited extraction; see national_reconciliation_2024.csv note).",
  "Baseline: IHME all-cause under-5 deaths and rates by GBD age group and state, both sexes, 2024 (data/external/nigeria_states).",
  "Contrast: HR_zero_vs_current = exp(f_g(0) - f_g(P_state)); attributable = IHME deaths x (1 - HR). Annual person-time held fixed.",
  "Intervals: conditional on fitted smoothing parameters, exposure and the fixed HIV imputation; state totals are point estimates and band limits are not summed.",
  "Comparator: IHME under-5 malaria death rate per 100,000 by state, 2024, applied to the all-cause implied person-years."),
  file.path(out,"model_specification.txt"))
cbh_atomic_csv(data.frame(file=inputs,md5=vapply(inputs,cbh_file_hash,"")),file.path(out,"input_provenance.csv"))
message("Nigerian state burden ",YEAR,": ",format(round(state_sum),big.mark=",")," attributable deaths across 37 states (national estimate ",
  format(round(national_total),big.mark=","),", difference ",sprintf("%+.2f%%",100*recon$relative_difference),"); IHME malaria ",
  format(round(sum(totals$ihme_malaria_deaths)),big.mark=","))
