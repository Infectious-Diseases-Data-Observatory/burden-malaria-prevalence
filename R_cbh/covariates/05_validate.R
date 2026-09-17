#!/usr/bin/env Rscript
# Focused extraction/selection checks; no model fits or respondent exports.
source("R_cbh/load_pipeline.R")
source("R_cbh/covariates/regional.R")
source("R_cbh/covariates/settings.R")
source("R_cbh/covariates/unicef.R")
library(data.table)
br <- data.frame(caseid=c("a","a","b","c"),v008=rep(1500,4),b3=c(1497,1476,1492,1420),
  v005=c(1,1,3,1)*1e6,b4=c(1,2,2,1),v025=c(1,1,2,2),v190=c(1,1,5,3),
  bord=c(2,1,1,1),v011=c(1150,1150,1100,1100),v133=c(6,6,12,4),b0=c(0,0,0,0),
  b5=c(1,1,1,1),b9=c(0,0,0,0),bidx=c(1,2,1,1),b19=c(3,24,8,80),
  m4=c(95,94,95,94),v409=c(0,0,0,0),v410=c(0,0,0,0),v411=c(0,0,0,0),v414e=c(0,0,0,0))
s <- data.frame(svkey="TEST",iso3="TST")
geo <- list(regkey=rep("one",4))
x <- cbh_regional_recode(br,s,NULL,geo)$values
get <- function(v)x$value[x$variable==v]
stopifnot(abs(get("male_pct")-20)<1e-12,
  abs(get("mean_wealth_quintile")-4)<1e-12,
  abs(get("mean_maternal_education_years")-10.5)<1e-12,
  get("exclusive_breastfeeding_pct")==100,
  x$eligible_n[x$variable=="exclusive_breastfeeding_pct"]==1L)
# A missing feeding response must not become an exclusive-breastfeeding success.
br$v409[1] <- NA_real_
y <- cbh_regional_recode(br,s,NULL,geo)$values
stopifnot(is.na(y$value[y$variable=="exclusive_breastfeeding_pct"]))
# Individual wealth missingness can coexist with a valid regional average.
br$v190[3] <- NA_real_
y <- cbh_regional_recode(br,s,NULL,geo)$values
stopifnot(y$value[y$variable=="mean_wealth_quintile"]==1,
  y$missing_n[y$variable=="mean_wealth_quintile"]==1L)
spec <- cbh_regional_spec(); vars <- all.vars(cbh_regional_formula())
# Missing regions do not contaminate a survey mean, nor borrow across surveys.
w <- data.frame(survey=c("A","A","A","B","B","C","C"),regkey=c("a","b","c","a","b","a","b"),
  example=c(20,NA,60,NA,NA,30,NA))
y <- cbh_regional_mean_fill(w,"example")
stopifnot(identical(cbh_regional_mean_fill(as.data.table(w),"example"),y))
stopifnot(identical(y$example,c(20,40,60,NA_real_,NA_real_,30,30)),
  identical(y$example_regional_mean_imputed,c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE)),
  identical(y$example_donor_regions,c(2L,2L,2L,0L,0L,1L,1L)))
bad <- try(cbh_regional_mean_fill(rbind(w,w[1,]),"example"),silent=TRUE)
stopifnot(inherits(bad,"try-error"))
stopifnot(length(spec$covariates)==22L,all(paste0("z_",spec$covariates) %in% vars),
  !any(c("sex","multiple_birth","birth_order","maternal_age_birth","wealth_quintile","urban") %in% vars))
args <- commandArgs(TRUE)
stopifnot(all(args %in% "--complete-case-only"))
settings <- cbh_covariate_settings("--complete-case-only" %in% args)
out <- settings$out
loss <- fread(file.path(out,"selection_by_survey_region.csv"))
stopifnot(all(loss$complete_case_records<=loss$records),
  all(loss$complete_case_children<=loss$children),
  all(loss$complete_case_deaths<=loss$deaths))
summary <- fread(file.path(out,"complete_case_summary.csv"))
stopifnot(all(summary$retained_regions+summary$lost_regions==summary$regions))
choices <- fread(file.path(out,"published_selections.csv"))
stopifnot(all(choices$consistent))
stopifnot(choices[survey=="CD61FL" & regkey=="equateur",all(source_labels=="Equateur")],
  choices[survey=="CD81FL" & regkey=="equateur",all(source_labels=="..Equateur (>= 2015)")])
writeLines(c("PASS: birth-weighted and distinct-mother-weighted summaries",
  "PASS: youngest eligible infant denominator and unknown feeding responses",
  "PASS: regional means retain records with missing individual values",
  "PASS: available-region means exclude NAs, preserve values, leave all-missing surveys unresolved and reject duplicate regions",
  "PASS: 22 regional/annual confounders; no individual confounding terms",
  "PASS: selection and region-loss accounting",
  "PASS: published source consistency and reviewed DRC boundary versions"),file.path(out,"validation.txt"))
cat("Regional adjustment validation passed.\n")
if(settings$unicef_fallback) {
  panel <- data.frame(country=c("TST","TST","TST"),year=c(2020,2021,2020),
    variable=c("dtp3_pct","dtp3_pct","hib3_pct"),value=c(70,90,0),source="UNICEF_WUENIC_country_year")
  input <- data.frame(country=c("TST","TST","TST"),survey_year=c(2020,2021,2019),
    entry_year=c(2020,2020,2019),dtp3_pct=c(NA,55,NA),hib3_pct=c(NA,30,NA))
  y <- cbh_unicef_fill(input,panel,c("dtp3_pct","hib3_pct"))
  stopifnot(identical(y$dtp3_pct,c(70,55,NA_real_)),identical(y$hib3_pct,c(0,30,NA_real_)),
    identical(y$dtp3_pct_imputed,c(TRUE,FALSE,FALSE)),
    identical(y$hib3_pct_imputed,c(TRUE,FALSE,FALSE)))
  bad <- try(cbh_unicef_fill(input,rbind(panel,panel[1,]),"dtp3_pct"),silent=TRUE)
  stopifnot(inherits(bad,"try-error"))
  w <- cbh_read_csv(file.path(settings$private,"regional_covariates_wide.csv"))
  p <- cbh_read_csv(settings$panel)
  for(v in c("dtp3_pct","measles_pct")) {
    fill <- w[[paste0(v,"_imputed")]]
    original <- w[[paste0(v,"_before_imputation")]]
    stopifnot(all(!fill | !is.finite(original)),identical(w[[v]][!fill],original[!fill]))
    z <- p[p$variable==v,]
    j <- match(paste(w$country,w$survey_year),paste(z$country,z$year))
    stopifnot(all(w[[v]][fill]==z$value[j[fill]]))
  }
  write("PASS: UNICEF exact-year matching, observed-value preservation, reported zeros and duplicate rejection",
    file.path(out,"validation.txt"),append=TRUE)
  cat("UNICEF substitution validation passed.\n")
}
