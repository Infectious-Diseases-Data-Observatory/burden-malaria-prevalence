#!/usr/bin/env Rscript
# Offline regional covariate assembly and complete-case accounting; no fits.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/covariates/regional.R")
source("R_cbh/covariates/settings.R")
source("R_cbh/covariates/unicef.R")
library(data.table)
args <- commandArgs(TRUE)
stopifnot(all(args %in% c("--complete-case-only","--legacy-expanded")))
settings <- cbh_covariate_settings("--complete-case-only" %in% args,"--legacy-expanded" %in% args)
private <- settings$base
stage_private <- settings$private
out <- settings$out
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(stage_private,recursive=TRUE,showWarnings=FALSE)
cfg <- cbh_config(); spec <- cbh_regional_spec(settings$expanded)
spec$id <- settings$id
unicef <- if(settings$unicef_fallback)cbh_read_csv(settings$panel) else NULL
cbh_atomic_rds(list(complete=FALSE,settings=settings,started=as.character(Sys.time())),
  file.path(stage_private,"manifest.rds"))
reg <- cbh_read_csv(cfg$registry)
bounds <- unique(cbh_read_csv(cfg$boundary_regions)[c("svkey","region","regkey")])
crosswalk <- cbh_read_csv("data/derived_cbh/region_crosswalk.csv")
pub <- cbh_read_csv(file.path(private,"published_indicators.csv"))
pub$survey <- reg$svkey[match(pub$SurveyId,reg$SurveyId)]
pub$clean_label <- trimws(gsub("[[:space:]]+"," ",gsub("(^| )(19|20)[0-9]{2}( |$)"," ",
  gsub("\\([^)]*(19|20)[0-9]{2}[^)]*\\)"," ",gsub("^[. ]+","",pub$CharacteristicLabel)))))
pub$value <- suppressWarnings(as.numeric(pub$Value))
pub$value[!is.finite(pub$value)|pub$value<0|pub$value>100] <- NA_real_
pub$observed_n <- suppressWarnings(as.numeric(pub$DenominatorUnweighted))
pub$weighted_n <- suppressWarnings(as.numeric(pub$DenominatorWeighted))
pub$value[is.finite(pub$observed_n) & pub$observed_n<=0] <- NA_real_
aliases <- cbh_read_csv("R_cbh/covariates/published_region_aliases.csv")
overrides <- cbh_read_csv("R_cbh/covariates/published_region_overrides.csv")
maps <- list()
for(s in unique(pub$survey)) {
  b <- bounds[bounds$svkey==s,,drop=FALSE]
  if(!nrow(b)) next
  labels <- unique(pub$clean_label[pub$survey==s & pub$level=="subnational"])
  x <- cbh_match_regions(labels,b)
  cw <- crosswalk[crosswalk$svkey==s & !is.na(crosswalk$regkey),]
  # Survey-specific recode matches are already reviewed and can resolve labels.
  for(j in which(is.na(x$regkey))) {
    keys <- unique(cw$regkey[cbh_rkey(cw$source_label)==cbh_rkey(x$source_label[j])])
    if(length(keys)==1) {x$regkey[j]<-keys;x$match_method[j]<-"reviewed_recode_crosswalk"}
  }
  ov <- overrides[overrides$survey==s,,drop=FALSE]
  j <- match(x$source_label,ov$source_label)
  good <- !is.na(j)
  if(any(good)) {
    stopifnot(all(ov$regkey[j[good]] %in% b$regkey))
    x$regkey[good] <- ov$regkey[j[good]]
    x$match_method[good] <- "survey_specific_published_override"
  }
  for(j in which(is.na(x$regkey))) {
    target <- aliases$region[aliases$api==cbh_rkey(x$source_label[j])]
    keys <- unique(b$regkey[cbh_rkey(b$region) %in% target | b$regkey %in% target])
    if(length(keys)==1) {x$regkey[j]<-keys;x$match_method[j]<-"curated_published_alias"}
  }
  x$survey <- s; maps[[s]] <- x
}
mapping <- cbh_bind(maps)
cbh_atomic_csv(mapping,file.path(out,"published_region_mapping.csv"))
pub$regkey <- mapping$regkey[match(paste(pub$survey,pub$clean_label),paste(mapping$survey,mapping$source_label))]
pub$match_method <- mapping$match_method[match(paste(pub$survey,pub$clean_label),paste(mapping$survey,mapping$source_label))]
versions <- cbh_read_csv("R_cbh/covariates/published_region_versions.csv")
indicator <- c(dtp3_pct="CH_VACC_C_DP3",measles_pct="CH_VACC_C_MSL",
  facility_delivery_pct="RH_DELP_C_DHF",exclusive_breastfeeding_pct="CN_BFSS_C_EBF",
  improved_water_pct="WS_SRCE_H_IMP",improved_sanitation_pct="WS_TLET_H_IMP",
  electricity_pct="HC_ELEC_H_ELC",interval_7_17_pct="FE_BINT_C_I07",interval_18_23_pct="FE_BINT_C_I18")
selected <- list(); choices <- list()
for(s in unique(pub$survey)) for(v in names(indicator)) {
  x <- pub[pub$survey==s & pub$level=="subnational" & pub$IndicatorId==indicator[[v]],]
  if(!nrow(x)) next
  if(v=="facility_delivery_pct") {
    windows <- c("Five years preceding the survey","Three years preceding the survey","Two years preceding the survey")
    window <- windows[windows %in% x$ByVariableLabel][1]
    x <- x[!is.na(x$ByVariableLabel) & x$ByVariableLabel==window,,drop=FALSE]
  } else x <- x[is.na(x$ByVariableLabel) | x$ByVariableLabel=="",,drop=FALSE]
  x <- x[!is.na(x$regkey),,drop=FALSE]
  for(r in unique(x$regkey)) {
    y <- x[x$regkey==r,,drop=FALSE]
    # A fine recode unit grouped to a broad mortality region must not compete
    # with a directly reported aggregate for that broad region.
    direct <- y$match_method %in% c("exact","canonical","curated_synonym","curated_published_alias","survey_specific_published_override")
    if(any(direct)) y <- y[direct,,drop=FALSE]
    choice <- versions$source_label[versions$survey==s & versions$regkey==r]
    if(length(choice)) {
      stopifnot(length(choice)==1L)
      y <- y[y$CharacteristicLabel==choice,,drop=FALSE]
      if(!nrow(y)) next
    }
    # Duplicate geographic presentations are allowed only if estimates AND
    # denominators agree; no fuzzy matching, arbitrary row, or national fallback.
    same <- nrow(unique(y[c("value","observed_n","weighted_n")]))==1L
    selected[[length(selected)+1L]] <- data.frame(survey=s,regkey=r,variable=v,
      value=if(same) y$value[1] else NA_real_,eligible_n=NA_integer_,
      observed_n=if(same)y$observed_n[1] else NA_real_,missing_n=NA_integer_,
      weighted_n=if(same)y$weighted_n[1] else NA_real_,population=indicator[[v]],
      source=if(same)"DHS_published_regional" else "ambiguous_published_region",
      small_denominator=if(same)is.finite(y$observed_n[1]) && y$observed_n[1]<25 else NA,
      low_precision=if(same)is.finite(y$observed_n[1]) && y$observed_n[1]<50 else NA)
    choices[[length(choices)+1L]] <- data.frame(survey=s,regkey=r,variable=v,
      source_rows=nrow(y),consistent=same,source_labels=paste(unique(y$CharacteristicLabel),collapse=";"),
      recall=paste(unique(y$ByVariableLabel),collapse=";"),data_ids=paste(y$DataId,collapse=";"))
  }
}
p <- as.data.table(cbh_bind(selected))
cbh_atomic_csv(cbh_bind(choices),file.path(out,"published_selections.csv"))
# Sum the two exhaustive published components of the <24-month category.
short <- merge(p[variable=="interval_7_17_pct"],p[variable=="interval_18_23_pct"],by=c("survey","regkey"),suffixes=c("","_second"))
short[, value:=ifelse(observed_n==observed_n_second & weighted_n==weighted_n_second,value+value_second,NA_real_)]
short[,variable:="short_birth_interval_pct"]
short[,population:="non_first_births_last_five_years"]
p <- rbind(p[!grepl("^interval_",variable)],short[,names(p),with=FALSE])
rec <- fread(file.path(private,"recode_regional.csv"))
feeding <- fread(file.path(private,"feeding_extraction.csv"))
national <- pub[pub$level=="national" & pub$IndicatorId=="CN_IYCB_C_EXB",]
stopifnot(!anyDuplicated(national$survey))
feeding[,published_national_pct:=national$value[match(survey,national$survey)]]
feeding[,difference_pp:=national_ebf_pct-published_national_pct]
# A large disagreement is a coding/denominator investigation, NOT a true zero.
# This conservative quality gate prevents demonstrably invalid extracted values
# from silently becoming confounders. Published regional estimates take priority.
feeding[,needs_questionnaire_review:=is.finite(difference_pp) & abs(difference_pp)>5]
rec[variable=="exclusive_breastfeeding_pct" & survey %in% feeding[needs_questionnaire_review==TRUE,survey],
  `:=`(value=NA_real_,source="feeding_questionnaire_review_required")]
cbh_atomic_csv(feeding,file.path(out,"breastfeeding_validation.csv"))
key <- function(d) paste(d$survey,d$regkey,d$variable)
# Official regional BF where available; otherwise retain validated recode estimate.
j <- match(key(rec),key(p)); replace <- !is.na(j) & is.finite(p$value[j])
rec <- rec[!replace]
p[,country:=reg$iso3[match(survey,reg$svkey)]]
regional <- rbindlist(list(rec,p),use.names=TRUE,fill=TRUE)
# Remove missing published BF duplicates when a recode value is the fallback.
regional <- regional[order(!is.finite(value))]
regional <- regional[!duplicated(key(regional))]
invalid <- regional[is.finite(value) & grepl("_pct$",variable) & (value< -1e-10 | value>100+1e-10)]
if(nrow(invalid)) {print(invalid);stop("Invalid regional percentage")}
regional[is.finite(value) & grepl("_pct$",variable),value:=pmin(100,pmax(0,value))]
cbh_unique(as.data.frame(regional),c("survey","regkey","variable"),"Regional covariates")
cbh_atomic_csv(regional,file.path(stage_private,"regional_covariates.csv"))
cbh_atomic_csv(regional,file.path(out,"regional_covariates.csv"))
wide <- dcast(regional,survey+regkey~variable,value.var="value")
for(v in setdiff(spec$regional,names(wide))) wide[,(v):=NA_real_]
if(!settings$expanded)wide <- wide[,c("survey","regkey",spec$regional),with=FALSE]
if(settings$unicef_fallback) {
  wide <- as.data.frame(wide)
  wide$country <- reg$iso3[match(wide$survey,reg$svkey)]
  wide$survey_year <- reg$year[match(wide$survey,reg$svkey)]
  wide <- cbh_unicef_fill(wide,unicef,c("dtp3_pct","measles_pct"))
  cbh_atomic_csv(wide[c("survey","country","survey_year","regkey",
    grep("^(dtp3_pct|measles_pct)",names(wide),value=TRUE))],file.path(out,"regional_vaccine_imputation.csv"))
}
if(settings$regional_mean) {
  wide <- cbh_regional_mean_fill(wide,spec$regional)
  cbh_atomic_csv(wide,file.path(out,"regional_mean_imputation.csv"))
}
cbh_atomic_csv(wide,file.path(stage_private,"regional_covariates_wide.csv"))
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds")); stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
hiv <- cbh_read_csv(cbh_trial_spec()$incidence_panel)
old_required <- c("death","age_band","pfpr_pct","calendar_year","band_years","survey","country","region",cbh_trial_spec()$covariates)
miss <- selections <- ages <- vaccines <- imputations <- list()
for(i in seq_len(nrow(m))) {
  object <- readRDS(file.path(cfg$output_dir,m$file[i]))
  stopifnot(identical(object$signature,m$signature[i]))
  d <- cbh_attach_incidence(object$data,hiv)
  old <- d$model_ready & complete.cases(d[old_required])
  for(v in old_required) if(is.numeric(d[[v]])) old <- old & is.finite(d[[v]])
  j <- match(paste(d$survey,d$regkey),paste(wide$survey,wide$regkey))
  for(v in spec$regional) d[[v]] <- wide[[v]][j]
  if(settings$unicef_fallback) {
    for(v in grep("_(before_imputation|imputed|source_year|imputation_source)$",names(wide),value=TRUE)) d[[v]]<-wide[[v]][j]
    d <- cbh_unicef_fill(d,unicef,intersect(spec$annual,names(cbh_unicef_indicators())))
    for(b in c("previous_primary","eligible_MAP")) for(v in intersect(spec$covariates,names(cbh_unicef_indicators()))) {
      mask <- if(b=="previous_primary")old else d$model_ready
      if(!any(mask))next
      z <- data.table(survey=d$survey[mask],country=d$country[mask],variable=v,baseline=b,
        source=d[[paste0(v,"_imputation_source")]][mask],imputed=d[[paste0(v,"_imputed")]][mask],
        missing=!is.finite(d[[v]][mask]))
      imputations[[length(imputations)+1L]]<-z[,.(records=.N,imputed_records=sum(imputed),missing_records=sum(missing)),
        by=.(baseline,survey,country,variable,source)]
    }
  }
  observed <- vapply(spec$covariates,function(v)is.finite(d[[v]]),logical(nrow(d)))
  full <- d$model_ready & rowSums(!observed)==0L
  regional_old <- spec$regional[1:7]
  original_annual <- spec$annual[1:4]
  baseline <- d$model_ready & rowSums(!observed[,c(regional_old,original_annual),drop=FALSE])==0L
  masks <- list(previous_primary=old,eligible_MAP=d$model_ready,regional_baseline=baseline)
  for(b in names(masks)) {
    keep <- masks[[b]]
    if(!any(keep)) next
    region <- d$region[keep]
    for(v in spec$covariates) {
      bad <- !observed[keep,v]
      a <- data.table(region,missing=bad)
      loss <- a[,.(any_missing=any(missing),all_missing=all(missing)),by=region]
      miss[[length(miss)+1L]] <- data.frame(baseline=b,survey=m$survey[i],country=m$country[i],
        variable=v,records=sum(keep),missing_records=sum(bad),regions=uniqueN(region),
        regions_with_any_missing=sum(loss$any_missing),regions_with_all_missing=sum(loss$all_missing))
    }
    a <- data.table(region,country=d$country[keep],survey=d$survey[keep],
      full=full[keep],death=d$death[keep],child=d$child_id[keep])
    z <- a[,.(records=.N,children=uniqueN(child),deaths=sum(death),
      complete_case_records=sum(full),complete_case_children=uniqueN(child[full]),
      complete_case_deaths=sum(death[full])),by=.(survey,country,region)]
    z[,baseline:=b]; selections[[length(selections)+1L]]<-z
    age <- data.table(age_band=d$age_band[keep],full=full[keep],death=d$death[keep])
    age <- age[,.(records=.N,deaths=sum(death),complete_case_records=sum(full),complete_case_deaths=sum(death[full])),by=age_band]
    age[,`:=`(baseline=b,survey=m$survey[i],country=m$country[i])];ages[[length(ages)+1L]]<-age
  }
  # Source status, including excluded pre-series/no-series zero placeholders.
  for(v in c("hib3_pct","pcv3_pct","rotavirus_pct")) {
    a <- data.table(survey=d$survey[old],country=d$country[old],year=d$entry_year[old],
      status=d[[paste0(v,"_status")]][old],variable=v)
    if(nrow(a)) vaccines[[length(vaccines)+1L]]<-a[,.(records=.N),by=.(survey,country,year,variable,status)]
  }
  if(i%%20L==0L)message("Audited ",i,"/",nrow(m)," surveys")
}
loss <- rbindlist(selections); missing <- rbindlist(miss)
summary <- loss[,.(records=sum(records),children=sum(children),deaths=sum(deaths),regions=.N,
  surveys=uniqueN(survey),countries=uniqueN(country),retained_records=sum(complete_case_records),
  retained_children=sum(complete_case_children),retained_deaths=sum(complete_case_deaths),
  retained_regions=sum(complete_case_records>0),lost_regions=sum(complete_case_records==0),
  partially_reduced_regions=sum(complete_case_records>0 & complete_case_records<records),
  retained_surveys=uniqueN(survey[complete_case_records>0]),retained_countries=uniqueN(country[complete_case_records>0])),by=baseline]
summary[,lost_regions_pct:=100*lost_regions/regions]
summary[,missing_records_pct:=100*(1-retained_records/records)]
stopifnot(summary[baseline=="previous_primary",records]==5885022L,
  summary[baseline=="previous_primary",deaths]==82415L,
  summary[baseline=="previous_primary",children]==1817912L,
  summary[baseline=="previous_primary",regions]==1015L)
tot <- missing[,lapply(.SD,sum),by=.(baseline,variable),.SDcols=c("records","missing_records","regions","regions_with_any_missing","regions_with_all_missing")]
tot[,`:=`(missing_records_pct=100*missing_records/records,
  regions_with_any_missing_pct=100*regions_with_any_missing/regions,
  regions_with_all_missing_pct=100*regions_with_all_missing/regions)]
country <- loss[,.(records=sum(records),retained_records=sum(complete_case_records),regions=.N,
  lost_regions=sum(complete_case_records==0),retained_regions=sum(complete_case_records>0)),by=.(baseline,country)]
cbh_atomic_csv(summary,file.path(out,"complete_case_summary.csv"))
cbh_atomic_csv(tot,file.path(out,"missingness_summary.csv"))
cbh_atomic_csv(missing,file.path(out,"missingness_by_survey.csv"))
cbh_atomic_csv(loss,file.path(out,"selection_by_survey_region.csv"))
cbh_atomic_csv(country,file.path(out,"selection_by_country.csv"))
cbh_atomic_csv(rbindlist(ages),file.path(out,"selection_by_survey_age.csv"))
cbh_atomic_csv(rbindlist(vaccines),file.path(out,"vaccine_source_status.csv"))
if(settings$unicef_fallback)cbh_atomic_csv(rbindlist(imputations),file.path(out,"imputation_by_survey.csv"))
paths <- c(file.path(stage_private,"regional_covariates_wide.csv"),
  file.path(private,c("published_indicators.csv","recode_regional.csv")),
  if(settings$unicef_fallback)settings$panel,
  "data/derived_cbh/manifest.rds",cbh_trial_spec()$incidence_panel,
  list.files("R_cbh/covariates",full.names=TRUE))
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"audit_input_manifest.csv"))
writeLines(paste(deparse(cbh_regional_formula(settings$expanded)),collapse=" "),file.path(out,"planned_formula.txt"))
cbh_atomic_rds(list(complete=TRUE,specification=spec,settings=settings,
  regional_overlay=file.path(stage_private,"regional_covariates_wide.csv"),
  input_manifest=meta$manifest,source_hashes=data.frame(file=paths,md5=vapply(paths,cbh_file_hash,""))),
  file.path(stage_private,"manifest.rds"))
print(summary)
