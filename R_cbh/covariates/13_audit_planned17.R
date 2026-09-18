#!/usr/bin/env Rscript
# Revised 17-variable availability audit; no mortality fitting or input overwrite.
source("R_cbh/load_pipeline.R")
source("R_cbh/covariates/regional.R")
source("R_cbh/analysis/model.R")
library(data.table)
cfg <- cbh_config()
private <- file.path(cfg$output_dir,"regional_adjustment/planned17_audit")
out <- "results/cbh/planned17_covariate_missingness"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
reg <- cbh_read_csv(cfg$registry)
bounds <- unique(cbh_read_csv(cfg$boundary_regions)[c("svkey","region","regkey")])
crosswalk <- cbh_read_csv("data/derived_cbh/region_crosswalk.csv")
pub <- cbh_read_csv(file.path(private,"nutrition_published.csv"));pub$level <- "subnational"
# Matching/selection rules copied from the existing audited regional assembler.
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
indicator <- c(stunting_pct="CN_NUTS_C_HA2",wasting_pct="CN_NUTS_C_WH2")
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
nut <- cbh_bind(selected)
cbh_atomic_csv(nut,file.path(out,"nutrition_selected.csv"))
cbh_atomic_csv(pub[is.na(pub$regkey),c("survey","IndicatorId","CharacteristicLabel","value")],
  file.path(out,"nutrition_unmatched.csv"))
oldwide <- cbh_read_csv("data/derived_cbh/regional_adjustment/reduced_v3/regional_covariates_wide.csv")
age <- cbh_read_csv(file.path(private,"first_birth_age.csv"))
newvars <- c("mean_maternal_age_first_birth","wasting_pct","stunting_pct")
regional <- c("mean_maternal_age_first_birth","mean_maternal_education_years","mean_wealth_quintile",
  "urban_pct","dtp3_pct","measles_pct","facility_delivery_pct","short_birth_interval_pct",
  "improved_water_pct","improved_sanitation_pct","electricity_pct","wasting_pct","stunting_pct")
annual <- c("log_hiv_incidence","log_gdp_pc","log_health_expenditure_pc","political_stability")
variables <- c(regional,annual);stopifnot(length(variables)==17L)
key <- function(d)paste(d$survey,d$regkey)
# Keep the existing overlay's entire donor frame, not only MAP-eligible regions.
wide <- oldwide
wide$mean_maternal_age_first_birth <- age$mean_maternal_age_first_birth[match(key(wide),key(age))]
for(v in c("wasting_pct","stunting_pct")) {
  n <- nut[nut$variable==v,];cbh_unique(n,c("survey","regkey"),"Nutrition summary")
  wide[[v]] <- n$value[match(key(wide),key(n))]
}
wide <- cbh_regional_mean_fill(wide,newvars)
cbh_atomic_csv(wide[c("survey","regkey",regional,
  grep("^(mean_maternal_age_first_birth|wasting_pct|stunting_pct)_",names(wide),value=TRUE))],
  file.path(private,"regional_covariates_wide.csv"))
# All covariates' pre-substitution values, including original vaccination data.
before <- wide
for(v in regional) before[[v]] <- wide[[paste0(v,if(v %in% c("dtp3_pct","measles_pct"))
  "_before_imputation" else "_before_regional_mean")]]
cbh_atomic_csv(age,file.path(out,"first_birth_age_summary.csv"))
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds"));stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
hiv_path <- cbh_trial_spec()$incidence_panel; hiv <- cbh_read_csv(hiv_path)
ledger <- selections <- summaries <- list(); child_counts <- list(); hiv_flags <- list()
for(i in seq_len(nrow(m))) {
  object <- readRDS(file.path(cfg$output_dir,m$file[i]));stopifnot(identical(object$signature,m$signature[i]))
  d <- cbh_attach_incidence(object$data,hiv)
  d <- d[d$model_ready,,drop=FALSE]
  if(!nrow(d))next
  j <- match(key(d),key(wide)); stopifnot(!anyNA(j))
  values <- as.data.frame(d[annual]); original <- values
  hj <- match(paste(d$country,d$entry_year),paste(hiv$iso3,hiv$year))
  child_rate <- hiv$child_rate[hj]
  original$log_hiv_incidence <- ifelse(is.finite(child_rate) & child_rate>0,log(pmax(child_rate,1e-100)),NA_real_)
  for(v in regional) {values[[v]] <- wide[[v]][j];original[[v]] <- before[[v]][j]}
  valid <- vapply(values[variables],is.finite,logical(nrow(d)))
  oldvalid <- vapply(d[annual],is.finite,logical(nrow(d)))
  oldkeep <- rowSums(!oldvalid)==0L
  for(v in cbh_regional_spec()$regional)oldkeep <- oldkeep & is.finite(oldwide[[v]][j])
  keep <- rowSums(!valid)==0L
  for(v in variables) {
    z <- data.table(survey=d$survey,country=d$country,region=d$region,
      before_missing=!is.finite(original[[v]]),after_missing=!is.finite(values[[v]]))
    ledger[[length(ledger)+1L]] <- z[,.(records=.N,before_missing=sum(before_missing),
      after_missing=sum(after_missing)),by=.(survey,country,region)][,variable:=v]
  }
  z <- data.table(survey=d$survey,country=d$country,region=d$region,keep,oldkeep)
  selections[[i]] <- z[,.(records=.N,retained_records=sum(keep),previous_records=sum(oldkeep),
    newly_excluded=sum(oldkeep & !keep),newly_recovered=sum(keep & !oldkeep)),by=.(survey,country,region)]
  child_counts[[i]] <- data.frame(survey=m$survey[i],eligible_children=uniqueN(d$child_id),
    retained_children=uniqueN(d$child_id[keep]),retained_deaths=sum(d$death[keep]),previous_records=sum(oldkeep))
  source_value <- hiv$child_source_value[hj]
  status <- ifelse(is.finite(child_rate) & child_rate>0,"observed_numeric",
    ifelse(!is.na(source_value) & grepl("^<",source_value),"censored","missing_or_invalid"))
  hiv_flags[[i]] <- as.data.table(table(status),keep.rownames=TRUE)
  rm(d,object,values,original,z,valid);gc(FALSE)
  if(i%%20==0)message("Missingness: ",i,"/",nrow(m)," surveys")
}
l <- rbindlist(ledger); sel <- rbindlist(selections); cc <- rbindlist(child_counts)
stopifnot(sum(sel$records)==6357802,nrow(sel)==1113,uniqueN(sel$survey)==120,
  sum(sel$previous_records)==5680117)
summary <- l[,.(records=sum(records),regions=.N,surveys=uniqueN(survey),
  before_missing_records=sum(before_missing),after_missing_records=sum(after_missing),
  before_regions_any_missing=sum(before_missing>0),after_regions_any_missing=sum(after_missing>0),
  before_regions_all_missing=sum(before_missing==records),after_regions_all_missing=sum(after_missing==records),
  before_surveys_any_missing=uniqueN(survey[before_missing>0]),after_surveys_any_missing=uniqueN(survey[after_missing>0])),by=variable]
summary[,`:=`(before_missing_pct=100*before_missing_records/records,after_missing_pct=100*after_missing_records/records)]
summary <- summary[match(variables,variable)]
summary[,level:=ifelse(variable %in% regional,"Survey-region","Country / band-entry year")]
stopifnot(nrow(summary)==17,all(summary$after_missing_records<=summary$before_missing_records))
# Retained unchanged variables must reproduce the existing eligible-MAP audit.
oldsum <- fread("results/cbh/regional_adjustment_reduced_v3/missingness_summary.csv")[baseline=="eligible_MAP"]
for(v in setdiff(variables,newvars))stopifnot(summary[variable==v,after_missing_records]==oldsum[variable==v,missing_records])
selection_summary <- data.frame(eligible_records=sum(sel$records),retained_records=sum(sel$retained_records),
  missing_records=sum(sel$records-sel$retained_records),eligible_regions=nrow(sel),
  retained_regions=sum(sel$retained_records>0),lost_regions=sum(sel$retained_records==0),
  partially_reduced_regions=sum(sel$retained_records>0 & sel$retained_records<sel$records),
  eligible_surveys=uniqueN(sel$survey),retained_surveys=uniqueN(sel[retained_records>0,survey]),
  eligible_countries=uniqueN(sel$country),retained_countries=uniqueN(sel[retained_records>0,country]),
  eligible_children=sum(cc$eligible_children),retained_children=sum(cc$retained_children),
  retained_deaths=sum(cc$retained_deaths),previous_records=sum(sel$previous_records),
  newly_excluded=sum(sel$newly_excluded),newly_recovered=sum(sel$newly_recovered))
cbh_atomic_csv(summary,file.path(out,"covariate_missingness.csv"))
cbh_atomic_csv(l,file.path(out,"missingness_by_survey_region.csv"))
cbh_atomic_csv(sel,file.path(out,"selection_by_survey_region.csv"))
cbh_atomic_csv(selection_summary,file.path(out,"complete_case_summary.csv"))
cbh_atomic_csv(rbindlist(hiv_flags)[,.(N=sum(N)),by=status],file.path(out,"hiv_original_value_status.csv"))
labels <- c("Maternal age at first birth","Maternal education (years)","Household wealth-quintile score",
  "Urban residence (%)","DTP3 coverage (%)","Measles coverage (%)","Facility delivery (%)",
  "Birth interval <24 months (%)","Improved drinking water (%)","Improved sanitation (%)",
  "Household electricity (%)","Wasting prevalence (%)","Stunting prevalence (%)",
  "Child HIV incidence (log)","GDP per capita (log)","Health expenditure per capita (log)","Political stability")
summary[,covariate:=labels];cbh_atomic_csv(summary,file.path(out,"covariate_missingness.csv"))
fmt <- function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)
writeLines(c("# Revised 17-covariate missingness audit","",
  "Denominator: 6,357,802 MAP-eligible child–age-band records, 1,113 survey-regions, 120 surveys and 36 countries, before complete-case selection. Percentages below are the unweighted share of records lacking their assigned regional/annual predictor, not missing individual responses. Regional missingness counts treat each survey-region once; national annual measures may be missing for only some entry years.","",
  "Before filling uses original regional values and observed numeric child HIV incidence. After filling applies the existing HIV posterior estimates, exact survey-year UNICEF DTP3/measles fallback and arithmetic mean of available regions in the same survey. Censored HIV values lack a usable numeric point estimate before filling and are counted separately in hiv_original_value_status.csv. No new HIV model, cross-survey borrowing, national nutrition fallback or mortality refit is performed.","",
  "| Covariate | Records missing before filling | Records missing after filling | Regions entirely missing after filling |",
  "|---|---:|---:|---:|",sprintf("| %s | %.2f%% | %.2f%% | %s |",labels,summary$before_missing_pct,summary$after_missing_pct,fmt(summary$after_regions_all_missing)),"",
  sprintf("Complete cases after substitution retain **%s records, %s children and %s deaths, in %s survey-regions, %s surveys and %s countries**. %s survey-regions have no remaining eligible records; %s additional regions retain only part of their records. Covariate-specific losses overlap and must not be summed.",fmt(selection_summary$retained_records),fmt(selection_summary$retained_children),fmt(selection_summary$retained_deaths),fmt(selection_summary$retained_regions),fmt(selection_summary$retained_surveys),fmt(selection_summary$retained_countries),fmt(selection_summary$lost_regions),fmt(selection_summary$partially_reduced_regions)),"",
  sprintf("Against the fitted 18-variable sample (%s records), the new specification excludes %s records and recovers %s. These are availability results; no 17-variable mortality model has been fitted.",fmt(selection_summary$previous_records),fmt(selection_summary$newly_excluded),fmt(selection_summary$newly_recovered)),"",
  "Age at first birth uses v212, valid completed ages 8–49, weighted by v005 among distinct interviewed mothers with a birth in the last 60 months. It is not age at each recent birth. Wasting and stunting use the official WHO-standard DHS indicators CN_NUTS_C_WH2 and CN_NUTS_C_HA2 (below −2 SD, including severe cases). Published regional denominators and reviewed geography rules are preserved; ambiguous or unmatched estimates remain unavailable before the same-survey fallback. Regional nutritional summaries describe measured surviving children at survey, not anthropometry of children who died.","",
  "[CSV table](covariate_missingness.csv) · [Complete-case counts](complete_case_summary.csv) · [Regional missingness](missingness_by_survey_region.csv) · [Regional selection](selection_by_survey_region.csv) · [Nutrition selections](nutrition_selected.csv) · [Unmatched nutrition labels](nutrition_unmatched.csv) · [Maternal summary denominators](first_birth_age_summary.csv)","",
  "Reproduce: Rscript R_cbh/covariates/11_extract_first_birth_age.R; python3 R_cbh/covariates/12_fetch_nutrition.py; Rscript R_cbh/covariates/13_audit_planned17.R. Public API retrieval is explicit and cached. New audit files are separate from all existing modelling data and results. No TeX files are written."),file.path(out,"README.md"))
paths <- c("R_cbh/covariates/11_extract_first_birth_age.R","R_cbh/covariates/12_fetch_nutrition.py",
  "R_cbh/covariates/13_audit_planned17.R","R_cbh/covariates/regional.R",cfg$registry,cfg$boundary_regions,
  file.path(private,c("nutrition_published.csv","nutrition_source_manifest.csv","first_birth_age.csv","first_birth_age_provenance.csv")),
  "data/derived_cbh/regional_adjustment/reduced_v3/regional_covariates_wide.csv",hiv_path,
  "R_cbh/covariates/published_region_aliases.csv","R_cbh/covariates/published_region_overrides.csv",
  "R_cbh/covariates/published_region_versions.csv",file.path(cfg$output_dir,"manifest.rds"))
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"provenance.csv"))
print(selection_summary);message("All 17 missingness columns verified; previous audit reproduced for retained variables.")
