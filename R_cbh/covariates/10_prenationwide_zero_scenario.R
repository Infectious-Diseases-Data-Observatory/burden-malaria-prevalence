#!/usr/bin/env Rscript
# Hypothetical missingness audit only: never overwrites model inputs or fits.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/covariates/regional.R")
source("R_cbh/covariates/settings.R")
source("R_cbh/covariates/unicef.R")
library(data.table)
cfg <- cbh_config(); spec <- cbh_regional_spec(expanded=TRUE); settings <- cbh_covariate_settings(expanded=TRUE)
args <- commandArgs(TRUE); stopifnot(all(args %in% "--region-mean"))
region_mean <- "--region-mean" %in% args
scenario_id <- if(region_mean)"prenationwide_zero_region_mean_scenario" else "prenationwide_zero_scenario"
out <- file.path(settings$out,scenario_id)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
wide <- fread(file.path(settings$private,"regional_covariates_wide.csv"))
filled_wide <- cbh_regional_mean_fill(wide,spec$regional)
if(region_mean) {
  private <- file.path(settings$private,scenario_id)
  dir.create(private,recursive=TRUE,showWarnings=FALSE)
  cbh_atomic_csv(filled_wide,file.path(private,"regional_covariates_wide.csv"))
}
intro_file <- file.path(settings$out,"vaccine_introduction_history.csv")
intro <- fread(intro_file)
stopifnot(!anyDuplicated(intro[,.(country,antigen)]))
vaccine_antigens <- c(hib3_pct="HIB",pcv3_pct="PNEUMO_CONJ",rotavirus_pct="ROTAVIRUS")
vaccines <- names(cbh_unicef_indicators())
panel <- cbh_read_csv(settings$panel)
hiv <- cbh_read_csv(cbh_trial_spec()$incidence_panel)
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds")); stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
required <- c("death","age_band","pfpr_pct","calendar_year","band_years",
              "survey","country","region",cbh_trial_spec()$covariates)
before_nationwide <- function(year,full_year,last_year,last_status) {
  (is.finite(full_year) & year<full_year) |
    (!is.finite(full_year) & !is.na(last_status) & last_status=="No" &
       is.finite(last_year) & year<=last_year)
}
stopifnot(identical(before_nationwide(c(2001,2002,2025,2026,2000),
  c(2002,2002,NA,NA,NA),c(2025,2025,2025,2025,NA),c("Yes","Yes","No","No",NA)),
  c(TRUE,FALSE,TRUE,FALSE,FALSE)))
pieces <- years <- regions <- scope <- list()
for(i in seq_len(nrow(m))) {
  object <- readRDS(file.path(cfg$output_dir,m$file[i]))
  stopifnot(identical(object$signature,m$signature[i]))
  d <- cbh_attach_incidence(object$data,hiv)
  old <- d$model_ready & complete.cases(d[required])
  for(v in required) if(is.numeric(d[[v]]))old <- old & is.finite(d[[v]])
  eligible <- d$model_ready
  old <- old[eligible]; d <- d[eligible,,drop=FALSE]
  if(!nrow(d))next
  j <- match(paste(d$survey,d$regkey),paste(wide$survey,wide$regkey))
  for(v in spec$regional)d[[v]] <- wide[[v]][j]
  d <- cbh_unicef_fill(d,panel,names(vaccine_antigens))
  bad_before <- vapply(spec$covariates,function(v)!is.finite(d[[v]]),logical(nrow(d)))
  bad_after <- bad_before
  j <- match(paste(d$survey,d$regkey),paste(filled_wide$survey,filled_wide$regkey))
  for(v in spec$regional) {
    donor_n <- filled_wide[[paste0(v,"_donor_regions")]][j]
    if(region_mean)bad_after[,v] <- !is.finite(filled_wide[[v]][j])
    a <- data.table(region=d$region,missing=bad_before[,v],donor_n=donor_n)
    a <- a[,.(records=.N,missing_records=sum(missing),
      recoverable_records=sum(missing & is.finite(donor_n) & donor_n>0),
      no_donor_records=sum(missing & (!is.finite(donor_n)|donor_n==0))),by=region]
    scope[[length(scope)+1L]] <- data.table(survey=m$survey[i],country=m$country[i],variable=v,
      regions=nrow(a),missing_regions=sum(a$missing_records>0),
      recoverable_regions=sum(a$recoverable_records>0),regions_without_donors=sum(a$no_donor_records>0),
      missing_records=sum(a$missing_records),recoverable_records=sum(a$recoverable_records),
      no_donor_records=sum(a$no_donor_records))
  }
  for(v in names(vaccine_antigens)) {
    history <- intro[antigen==vaccine_antigens[[v]]]
    j <- match(d$country,history$country); stopifnot(!anyNA(j))
    zero <- bad_before[,v] & before_nationwide(d$entry_year,
      history$first_reported_nationwide_year[j],history$last_reported_year[j],history$last_reported_status[j])
    bad_after[zero,v] <- FALSE
    stopifnot(all(!zero | bad_before[,v])) # Only missing values can change.
  }
  add_joint <- function(bad)cbind(bad,any_vaccine=rowSums(bad[,vaccines,drop=FALSE])>0,
                                  any_required=rowSums(bad)>0)
  bad_before <- add_joint(bad_before); bad_after <- add_joint(bad_after)
  stopifnot(all(!bad_after | bad_before))
  for(b in c("previous_primary","eligible_MAP")) {
    keep <- if(b=="previous_primary")old else rep(TRUE,nrow(d))
    if(!any(keep))next
    for(v in colnames(bad_before)) {
      before <- bad_before[keep,v]; after <- bad_after[keep,v]
      pieces[[length(pieces)+1L]] <- data.table(baseline=b,survey=m$survey[i],country=m$country[i],
        survey_year=paste(sort(unique(d$survey_year[keep])),collapse="/"),variable=v,
        records=sum(keep),missing_before=sum(before),missing_after=sum(after),resolved_records=sum(before & !after))
      a <- data.table(year=d$entry_year[keep],before,after)
      a <- a[,.(records=.N,missing_before=sum(before),missing_after=sum(after)),by=year]
      a[,`:=`(baseline=b,variable=v)]; years[[length(years)+1L]] <- a
    }
    a <- data.table(region=d$region[keep],child=d$child_id[keep],death=d$death[keep],
      full=!bad_after[keep,"any_required"],vaccine_full=!bad_after[keep,"any_vaccine"])
    a <- a[,.(records=.N,retained_records=sum(full),children=uniqueN(child),
      retained_children=uniqueN(child[full]),deaths=sum(death),retained_deaths=sum(death[full]),
      vaccine_missing_records=sum(!vaccine_full)),by=region]
    a[,`:=`(baseline=b,survey=m$survey[i],country=m$country[i])]
    regions[[length(regions)+1L]] <- a
  }
  if(i%%20L==0L)message("Audited ",i,"/",nrow(m)," surveys")
}
survey <- rbindlist(pieces)
summary <- survey[,.(records=sum(records),missing_before=sum(missing_before),missing_after=sum(missing_after),
  resolved_records=sum(resolved_records),surveys=.N,surveys_with_missing=sum(missing_after>0),
  countries_with_missing=uniqueN(country[missing_after>0])),by=.(baseline,variable)]
summary[,`:=`(missing_after_pct=100*missing_after/records,
  originally_missing_remaining_pct=ifelse(missing_before>0,100*missing_after/missing_before,NA_real_))]
# Reconcile starting values with the established UNICEF-fallback audit.
reference <- fread(file.path(settings$out,"missingness_summary.csv"))
check <- merge(summary[variable %in% spec$covariates],reference,by=c("baseline","variable"),suffixes=c("","_reference"))
stopifnot(nrow(check)==2*length(spec$covariates),all(check$records==check$records_reference),
          all(check$missing_before==check$missing_records))
reference <- fread(file.path(settings$out,"complete_case_summary.csv"))
check <- merge(summary[variable=="any_required"],reference,by="baseline",suffixes=c("","_reference"))
stopifnot(all(check$missing_before==check$records_reference-check$retained_records))
country <- survey[,.(records=sum(records),missing_before=sum(missing_before),missing_after=sum(missing_after),
  surveys=.N,surveys_with_missing=sum(missing_after>0),
  surveys_affected=paste(survey[missing_after>0],collapse="; "),
  survey_years_affected=paste(survey_year[missing_after>0],collapse="; ")),
  by=.(baseline,country,variable)]
country[,country_name:=intro$country_name[match(country,intro$country)]]
annual <- rbindlist(years)[,.(records=sum(records),missing_before=sum(missing_before),missing_after=sum(missing_after)),
  by=.(baseline,variable,year)]
annual[,missing_after_pct:=100*missing_after/records]
region <- rbindlist(regions)
selection <- region[,.(records=sum(records),retained_records=sum(retained_records),
  retained_children=sum(retained_children),retained_deaths=sum(retained_deaths),
  regions=.N,retained_regions=sum(retained_records>0),
  retained_surveys=uniqueN(survey[retained_records>0]),retained_countries=uniqueN(country[retained_records>0])),by=baseline]
for(n in c("survey","summary","country","annual","region","selection"))cbh_atomic_csv(get(n),file.path(out,paste0(n,".csv")))
scope <- rbindlist(scope)
cbh_atomic_csv(scope,file.path(out,"regional_missingness_scope_by_survey.csv"))
scope_summary <- scope[,.(missing_records=sum(missing_records),recoverable_records=sum(recoverable_records),
  no_donor_records=sum(no_donor_records),missing_regions=sum(missing_regions),
  recoverable_regions=sum(recoverable_regions),regions_without_donors=sum(regions_without_donors),
  surveys_with_recoverable_regions=sum(recoverable_regions>0),surveys_without_donors=sum(regions_without_donors>0)),by=variable]
cbh_atomic_csv(scope_summary,file.path(out,"regional_missingness_scope_summary.csv"))
inputs <- c(intro_file,file.path(cfg$output_dir,"manifest.rds"),file.path(settings$private,"manifest.rds"),
  file.path(settings$private,"regional_covariates_wide.csv"),settings$panel,cbh_trial_spec()$incidence_panel)
cbh_atomic_csv(data.frame(file=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,"input_manifest.csv"))
print(summary[variable %in% c(vaccines,"any_vaccine","any_required")])
print(country[variable=="any_vaccine" & missing_after>0])
print(selection)

# Standalone report from the saved aggregate outputs.
report_dir <- file.path(cbh_covariate_settings(expanded=TRUE)$out,scenario_id)
s <- fread(file.path(report_dir,"summary.csv"))
c <- fread(file.path(report_dir,"country.csv"))
q <- fread(file.path(report_dir,"survey.csv"))
fmt <- function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)
lines <- c("# Missingness under the pre-nationwide zero scenario", "",
  "Scenario dated 17 September 2026. Apply zero only to missing Hib3, PCV and rotavirus coverage strictly before the first WHO-reported nationwide introduction year. Preserve all observed coverage. Where introduction is not recorded and the latest status is No, apply zero through that last reported year (2025). This includes missing values during partial rollout by user specification; it is not a claim that true coverage was zero. Introduction-year gaps remain missing. Existing UNICEF DTP3/measles fallback is retained. No base datasets, model inputs or fits are overwritten.", "",
  "Records mean child–age-band rows. The previous-primary denominator matches the earlier missingness tables; the full MAP-eligible denominator is the starting pool for the revised regional adjustment model. These count surveys with at least one missing record, not necessarily surveys completely excluded.", "",
  "| Sample | Records | Remaining vaccine gaps | % of records | % of previously vaccine-missing records | Surveys with vaccine gaps | Countries |",
  "|---|---:|---:|---:|---:|---:|---:|")
if(region_mean)lines <- c("# Missingness with pre-nationwide zeros and available-region means","",
  "Additional scenario: each missing regional covariate receives the arithmetic mean of finite available regional estimates within the same survey, calculated separately for each variable. Every donor region counts once; child-band row counts do not weight this mean. Surveys with no available donor region remain missing. National country-year covariates do not use this fallback. Existing values are preserved and donor counts/imputation flags are retained in the private wide overlay. This scenario does not change the adopted primary data policy.","",lines[-c(1,2)])
for(i in which(s$variable=="any_vaccine"))lines <- c(lines,sprintf("| %s | %s | %s | %.2f%% | %.2f%% | %s | %s |",
  s$baseline[i],fmt(s$records[i]),fmt(s$missing_after[i]),s$missing_after_pct[i],s$originally_missing_remaining_pct[i],
  s$surveys_with_missing[i],s$countries_with_missing[i]))
lines <- c(lines,"","## Surveys with remaining vaccine gaps: full MAP-eligible pool","",
  "Survey years are registry start years; surveys can span two calendar years. Kenya 2003 is additional to the previous-primary sample; the other countries/surveys are present in both samples.","",
  "| Country | Number of surveys | Survey start years | Missing vaccine(s) | Records with vaccine gaps |",
  "|---|---:|---|---|---:|")
labels <- c(hib3_pct="Hib3",pcv3_pct="PCV",rotavirus_pct="Rotavirus")
z <- c[baseline=="eligible_MAP" & variable=="any_vaccine" & missing_after>0][order(country_name)]
for(i in seq_len(nrow(z))) {
  cc <- z$country[i]
  vv <- q[baseline=="eligible_MAP" & country==cc & variable %in% names(labels) & missing_after>0,unique(variable)]
  lines <- c(lines,sprintf("| %s | %s | %s | %s | %s |",z$country_name[i],z$surveys_with_missing[i],
    z$survey_years_affected[i],paste(labels[vv],collapse=", "),fmt(z$missing_after[i])))
}
lines <- c(lines,"","## All 22 required covariates","",
  "| Sample | Records with any remaining gap | % of records | Surveys with any gap | Countries |",
  "|---|---:|---:|---:|---:|")
for(i in which(s$variable=="any_required"))lines <- c(lines,sprintf("| %s | %s | %.2f%% | %s | %s |",
  s$baseline[i],fmt(s$missing_after[i]),s$missing_after_pct[i],s$surveys_with_missing[i],s$countries_with_missing[i]))
lines <- c(lines,"","Countries/surveys with any remaining adjustment-variable gap in the full MAP-eligible pool:","",
  "| Country | Surveys with any gap | Survey start years |","|---|---:|---|")
z <- c[baseline=="eligible_MAP" & variable=="any_required" & missing_after>0][order(country_name)]
for(i in seq_len(nrow(z)))lines <- c(lines,sprintf("| %s | %s | %s |",z$country_name[i],z$surveys_with_missing[i],z$survey_years_affected[i]))
lines <- c(lines,"","[Variable totals](summary.csv), [individual surveys and covariates](survey.csv), [country totals](country.csv), [complete-case selection](selection.csv), [annual missingness](annual.csv).", "",
  "Starting missingness was reconciled with the established UNICEF-fallback audit for both samples and all covariates. Boundary checks enforce strictly before nationwide rollout and preservation of observed values. Source hashes are in input_manifest.csv.", "",
  paste0("Reproduce: `Rscript R_cbh/covariates/10_prenationwide_zero_scenario.R",if(region_mean)" --region-mean" else "","`."))
if(region_mean) {
  z <- fread(file.path(report_dir,"regional_missingness_scope_summary.csv"))[missing_records>0]
  lines <- c(lines,"","## Where the regional gaps occur","",
    "Full MAP-eligible pool, before regional-mean substitution. Counts across covariates overlap.","",
    "| Covariate | Regions fillable from other regions | Surveys with fillable regions | Surveys without any donor region |",
    "|---|---:|---:|---:|")
  for(i in seq_len(nrow(z)))lines <- c(lines,sprintf("| %s | %s | %s | %s |",z$variable[i],
    z$recoverable_regions[i],z$surveys_with_recoverable_regions[i],z$surveys_without_donors[i]))
  lines <- c(lines,"",paste0("Within-region means already exclude missing individual responses and invalid weights. ",
    "A regional mean is missing only when it cannot be estimated from usable source values under the existing extraction/quality rules. ",
    "Available-region imputation does not solve whole-survey absence, questionnaire-review exclusions affecting all donor regions, or missing national country-year inputs."))
}
writeLines(lines,file.path(report_dir,"README.md"))
