#!/usr/bin/env Rscript
# Attribute every MAP-eligible child-band record excluded by complete-case
# selection to the FIRST missing covariate in a declared order, so the reasons
# are disjoint and sum to the total. Aggregate counts only; no model change.
source("R_cbh/load_pipeline.R")
source("R_cbh/covariates/regional.R")
source("R_cbh/analysis/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
library(data.table)
# --liberia (24 September 2026): the same attribution with the HIV panel that adds Liberia
# (R_cbh/hiv/03_add_liberia_aidsinfo.R), for the DHS part of the v7 primary; written to its own folder.
args <- commandArgs(trailingOnly=TRUE); stopifnot(all(args %in% "--liberia")); liberia <- "--liberia" %in% args
cfg <- cbh_config(); primary <- cbh_primary_settings("regional"); stopifnot(isTRUE(primary$nutrition))
private <- file.path(cfg$output_dir,"regional_adjustment/planned17_audit")
audit <- "results/cbh/planned17_covariate_missingness"
out <- if(liberia) paste0(audit,"_lbr") else audit
dir.create(out,recursive=TRUE,showWarnings=FALSE)
spec <- cbh_primary_regional_spec(primary)
# National annual series first (country-year gaps), then the regional summaries
# ordered by how often they are missing in the availability audit.
order_vars <- c("log_hiv_incidence","log_health_expenditure_pc","political_stability","log_gdp_pc",
  "wasting_pct","stunting_pct","facility_delivery_pct","electricity_pct","mean_wealth_quintile",
  "dtp3_pct","measles_pct","short_birth_interval_pct","improved_water_pct","improved_sanitation_pct",
  "urban_pct","mean_maternal_education_years","mean_maternal_age_first_birth")
stopifnot(setequal(order_vars,spec$covariates))
wide <- cbh_read_csv(file.path(private,"regional_covariates_wide.csv"))
cbh_unique(wide,c("survey","regkey"),"Regional overlay")
hiv <- cbh_read_csv(if(liberia) cbh_primary_settings("regional_mics")$hiv_panel else cbh_trial_spec()$incidence_panel)
meta <- readRDS(file.path(cfg$output_dir,"manifest.rds")); stopifnot(meta$complete)
m <- meta$manifest[meta$manifest$status %in% c("built","cached"),]
key <- function(d) paste(d$survey,d$regkey)
tallies <- list(); eligible <- retained <- 0
for(i in seq_len(nrow(m))) {
  object <- readRDS(file.path(cfg$output_dir,m$file[i])); stopifnot(identical(object$signature,m$signature[i]))
  d <- cbh_attach_incidence(object$data,hiv)
  d <- d[d$model_ready,,drop=FALSE]
  if(!nrow(d)) next
  j <- match(key(d),key(wide)); stopifnot(!anyNA(j))
  for(v in spec$regional) d[[v]] <- wide[[v]][j]
  valid <- vapply(order_vars,function(v) is.finite(d[[v]]),logical(nrow(d)))
  first <- apply(!valid,1,function(r) if(any(r)) order_vars[which(r)[1]] else "retained")
  eligible <- eligible+nrow(d); retained <- retained+sum(first=="retained")
  tallies[[i]] <- data.table(survey=m$survey[i],country=m$country[i],reason=first,death=d$death)[
    ,.(records=.N,deaths=sum(death)),by=.(survey,country,reason)]
  rm(d,object,valid,first); gc(FALSE)
  if(i%%20==0) message("Attribution: ",i,"/",nrow(m)," surveys")
}
t <- rbindlist(tallies)
summary <- t[reason!="retained",.(records=sum(records),deaths=sum(deaths),surveys=uniqueN(survey),countries=uniqueN(country)),by=reason]
summary <- summary[match(intersect(order_vars,summary$reason),reason)]
summary[,share_of_excluded_pct:=100*records/sum(records)]
cc <- cbh_read_csv(file.path(audit,"complete_case_summary.csv"))
if(liberia) {
  # Eligibility is unchanged; the retained DHS records must equal the v7 DHS part.
  cc$retained_records <- cbh_primary_settings("regional_mics")$expected_dhs_records
  cc$missing_records <- cc$eligible_records-cc$retained_records
}
stopifnot(eligible==cc$eligible_records,retained==cc$retained_records,sum(summary$records)==cc$missing_records)
labels <- c(log_hiv_incidence="Child HIV incidence (national series unavailable)",
  log_health_expenditure_pc="Health expenditure per capita",political_stability="Political stability",
  log_gdp_pc="GDP per capita",wasting_pct="Wasting prevalence",stunting_pct="Stunting prevalence",
  facility_delivery_pct="Facility delivery",electricity_pct="Household electricity",
  mean_wealth_quintile="Wealth-quintile score",dtp3_pct="DTP3 coverage",measles_pct="Measles coverage",
  short_birth_interval_pct="Short birth interval",improved_water_pct="Improved water",
  improved_sanitation_pct="Improved sanitation",urban_pct="Urban residence",
  mean_maternal_education_years="Maternal education",mean_maternal_age_first_birth="Maternal age at first birth")
summary[,label:=labels[reason]]
cbh_atomic_csv(as.data.frame(summary),file.path(out,"exclusion_attribution.csv"))
cbh_atomic_csv(as.data.frame(t[reason!="retained"]),file.path(out,"exclusion_attribution_by_survey.csv"))
writeLines(c("# Disjoint attribution of covariate exclusions","",
  sprintf("Of %s MAP-eligible child-band records, %s were excluded because at least one of the 17 covariates was unavailable after the HIV, UNICEF vaccination and available-region substitutions. Each excluded record is attributed here to the first missing covariate in the order: national annual series (child HIV incidence, health expenditure, political stability, GDP), then regional summaries (wasting, stunting, facility delivery, electricity, wealth, vaccination, birth interval, water, sanitation, urban, education, age at first birth). Counts are therefore disjoint and sum to the total; a record missing several covariates is counted once, under the first.",
    format(cc$eligible_records,big.mark=","),format(cc$missing_records,big.mark=",")),"",
  "| Attributed reason | Records | Deaths | Surveys affected | Share of excluded |","|---|---:|---:|---:|---:|",
  sprintf("| %s | %s | %s | %d | %.1f%% |",summary$label,format(summary$records,big.mark=",",trim=TRUE),
    format(summary$deaths,big.mark=",",trim=TRUE),summary$surveys,summary$share_of_excluded_pct),"",
  paste0(if(liberia) "HIV panel with Liberia's child incidence derived from UNAIDS counts (DHS part of the v7 primary). " else "",
    "Overlapping per-variable counts are in ",if(liberia) "../planned17_covariate_missingness/" else "","covariate_missingness.csv; per-survey attribution in exclusion_attribution_by_survey.csv. Reproduce: `Rscript R_cbh/covariates/14_exclusion_attribution.R",if(liberia) " --liberia" else "","`.")),
  file.path(out,"EXCLUSION_ATTRIBUTION.md"))
print(as.data.frame(summary[,.(label,records,deaths,surveys,share_of_excluded_pct)]))
message("Exclusion attribution complete")
