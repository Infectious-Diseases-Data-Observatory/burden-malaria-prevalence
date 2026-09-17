#!/usr/bin/env Rscript
source("R_cbh/load_pipeline.R")
source("R_cbh/covariates/settings.R")
library(data.table)
args <- commandArgs(TRUE)
stopifnot(all(args %in% "--complete-case-only"))
settings <- cbh_covariate_settings("--complete-case-only" %in% args)
out <- settings$out
m <- fread(file.path(out,"missingness_summary.csv"))[baseline=="previous_primary"]
s <- fread(file.path(out,"complete_case_summary.csv"))
f <- fread(file.path(out,"breastfeeding_validation.csv"))
regions <- fread(file.path(out,"selection_by_survey_region.csv"))[baseline=="previous_primary"]
regions[,regkey:=sub("^[^:]*:[^:]*:","",region)]
r <- fread(file.path(out,"regional_covariates.csv"))
b <- merge(regions,r[variable=="exclusive_breastfeeding_pct"],by=c("survey","regkey"))
labels <- c(male_pct="Male births",multiple_birth_pct="Multiple births",mean_birth_order="Mean birth order",
  mean_maternal_age_birth="Mean maternal age at birth",mean_maternal_education_years="Mean maternal education",
  mean_wealth_quintile="Mean wealth-quintile score",urban_pct="Urban residence",
  dtp3_pct="DTP3 coverage",measles_pct="Measles coverage",facility_delivery_pct="Facility delivery",
  exclusive_breastfeeding_pct="Exclusive breastfeeding",short_birth_interval_pct="Birth interval <24 months",
  improved_water_pct="Improved drinking water",improved_sanitation_pct="Improved sanitation",
  electricity_pct="Electricity",log_hiv_incidence="Child HIV incidence (fixed imputation)",
  log_gdp_pc="GDP per capita",log_health_expenditure_pc="Health expenditure",political_stability="Political stability",
  hib3_pct="Hib3 coverage",pcv3_pct="PCV coverage",rotavirus_pct="Rotavirus coverage")
comma <- function(x)format(x,big.mark=",",trim=TRUE,scientific=FALSE)
old <- s[baseline=="previous_primary"]; new <- s[baseline=="eligible_MAP"]
table <- sprintf("| %s | %.2f%% | %s |",labels[m$variable],m$missing_records_pct,comma(m$regions_with_all_missing))
overview <- sprintf("| %s | %s | %s | %s | %.1f%% |",s$baseline,comma(s$records),comma(s$regions),comma(s$lost_regions),s$lost_regions_pct)
source <- b[!is.finite(value),.(regions=.N,records=sum(records)),by=source]
small <- b[complete_case_records>0,sum(small_denominator,na.rm=TRUE)]
low <- b[complete_case_records>0,sum(low_precision,na.rm=TRUE)]
lines <- c("# Regional adjustment availability audit — 17 September 2026","",
  sprintf("**Complete cases retain %s of the previous %s survey-regions: %s are lost (%.1f%%).** A further %s retained regions lose some records. Region counts distinguish surveys/boundary versions.",
    comma(old$retained_regions),comma(old$regions),comma(old$lost_regions),old$lost_regions_pct,comma(old$partially_reduced_regions)),"",
  sprintf("Starting from all eligible MAP records, the revised regional adjustment set yields **%s child-band records, %s distinct children and %s deaths, from %s surveys in %s countries and %s survey-regions**. This is an availability result, not a fitted model.",
    comma(new$retained_records),comma(new$retained_children),comma(new$retained_deaths),comma(new$retained_surveys),comma(new$retained_countries),comma(new$retained_regions)),"",
  sprintf("Within the previous fitted sample, %s records remain and %.2f%% are excluded. Starting before its individual-level complete-case filter recovers %s additional records (%s children) using regional means.",
    comma(old$retained_records),old$missing_records_pct,comma(new$retained_records-old$retained_records),comma(new$retained_children-old$retained_children)),"",
  "## Missingness relative to the previous fitted sample","",
  "Denominator: **5,885,022 child-band records and 1,015 survey-regions**. Percentages are the unweighted share of mortality records lacking the assigned regional or annual covariate, not the percentage of respondents with an unanswered question. A regional mean may remain available despite some missing responses. The final column counts regions with no usable value for that covariate in any of their previously included records. Annual covariates can also remove part of a region's records.","",
  "| Covariate | Records missing | Entirely unavailable survey-regions |","|---|---:|---:|",table,"",
  "Losses overlap; neither percentages nor regional losses should be added across covariates. All eleven earlier adjustment concepts have 0% missingness within the previous fitted sample by its original selection; the regional versions also have complete coverage on that denominator.","",
  "## Alternative denominators","",
  "| Starting sample | Records | Survey-regions | Regions lost | Regions lost (%) |","|---|---:|---:|---:|---:|",overview,"",
  "`eligible_MAP` is every valid complete-band record with geography and exact entry-year MAP exposure. `regional_baseline` requires the regional versions of the previous adjustment concepts and the four existing annual predictors, before adding vaccines/delivery/feeding/interval/WASH/electricity. These definitions explain why regional averaging can increase initial record availability without adding a region to the former fitted sample.","",
  "## Extraction findings and unresolved issues","",
  "- All five vaccine measures are selected. National Hib3/PCV/rotavirus missingness is predominantly pre-series placeholders previously labelled assumed not introduced, plus countries with no series. These are not independently verified zero-coverage estimates. They remain missing under the existing pipeline rule. Resolving introduction dates/source coverage could substantially change these losses.",
  "- Published DHS regional vaccination, delivery, water, sanitation and electricity replace fragile birth-recode label matching. Historical boundary alternatives and nested source regions are resolved explicitly; individual household coverage is not estimated using women's weights.",
  "- Birth interval is a regional percentage among non-first births. First births are not excluded from the mortality model because B11 is structurally unavailable.",
  sprintf("- Breastfeeding remains unavailable for %s previously included regions: %s have no usable recode estimate and %s have recode estimates held out for questionnaire review. Published regional estimates take precedence where available. Ten surveys trigger the >5-percentage-point national discrepancy screen; their source labels, active items and comparison values are in `breastfeeding_validation.csv`. This is a conservative implementation screen, not a calibration or proof that the remaining questionnaires are equivalent.",
    comma(sum(source$regions)),comma(source[source=="local_BR_weighted_summary",regions]),comma(source[source=="feeding_questionnaire_review_required",regions])),
  sprintf("- Among retained regions, %s breastfeeding estimates have fewer than 25 observations and %s have fewer than 50. They remain included in this availability audit; small-denominator uncertainty and a declared sensitivity need assessment before final inference.",comma(small),comma(low)),
  "- Some absent regional values still need source/geography review. The full mapping ledger includes unmatched alternative subregions outside the model geography; these are not all lost analysis regions. Use `selection_by_survey_region.csv` and the wide covariate table to identify the actual losses.",
  "- Facility-delivery recall can be five, three or two years. DHS covariates are survey-time summaries; national annual covariates retain band-entry-year assignment. Regional averages change the adjustment interpretation and do not control individual-level confounding in the same way as the earlier model.","",
  "## Files and reproduction","",
  "- [Complete-case summary](complete_case_summary.csv), [missingness by variable](missingness_summary.csv), [survey-region ledger](selection_by_survey_region.csv), [country losses](selection_by_country.csv), and [age/survey selection](selection_by_survey_age.csv).",
  "- [Regional values and denominators](regional_covariates.csv), [breastfeeding validation](breastfeeding_validation.csv), [published regional mapping](published_region_mapping.csv), [source choices](published_selections.csv), and [vaccine source statuses](vaccine_source_status.csv).",
  "- [Pipeline and definitions](../../../R_cbh/covariates/README.md). The wide overlay is `data/derived_cbh/regional_adjustment/regional_covariates_wide.csv`. Raw recodes and respondent identifiers remain under ignored `data/`.",
  "- The audit reproduces all old sample totals exactly and checks unique joins, equal denominators for interval components, valid percentages and non-overlapping complete-case counts. No mortality fits, existing figures or TeX files were changed.","",
  "Sources: [DHS API](https://api.dhsprogram.com/), the [DHS reference feeding code](https://github.com/DHSProgram/DHS-Indicators-R/blob/main/Chap11_NT/NT_IYCF.R), and [DHS youngest-child selection](https://github.com/DHSProgram/DHS-Indicators-R/blob/main/Chap11_NT/!NTmain.R). Public source URLs and hashes are retained in the data snapshot's `source_manifest.csv`; local inputs are fingerprinted in `audit_input_manifest.csv`.")
writeLines(lines,file.path(out,"README.md"))
if(settings$unicef_fallback) {
  imp <- fread(file.path(out,"imputation_by_survey.csv"))[baseline=="previous_primary",
    .(imputed_records=sum(imputed_records),remaining_missing=sum(missing_records)),by=variable]
  cbh_atomic_csv(imp,file.path(out,"imputation_summary.csv"))
  intro <- c("**Primary policy: national UNICEF vaccination fallback.** Preserve observed regional DTP3/measles values; substitute the exact country/survey-year WUENIC estimate only when the regional value is missing. The three existing national vaccines use band-entry year. Pre-series/no-series assumptions remain missing unless independently verified; no interpolation, nearest-year carry or automatic zero filling is applied.","",
    "| Vaccine | Previously included records imputed | Records still missing |",
    "|---|---:|---:|",sprintf("| %s | %s | %s |",labels[imp$variable],comma(imp$imputed_records),comma(imp$remaining_missing)),"",
    "This is deterministic national substitution, not a statistical multiple-imputation model. Regional heterogeneity and uncertainty in the substituted estimates are not propagated. Preserve the [pre-imputation audit](../regional_adjustment_v1/README.md) as the comparison. The regional covariate long table describes source estimates before fallback; `regional_vaccine_imputation.csv` and the versioned wide overlay contain the final vaccination values and flags.","")
  lines <- append(lines,intro,after=2L)
  lines <- sub("## Missingness relative to the previous fitted sample",
    "## Remaining missingness after UNICEF substitution",lines,fixed=TRUE)
  lines <- sub("requires complete cases","requires complete cases after substitution",lines,fixed=TRUE)
  writeLines(lines,file.path(out,"README.md"))
}
