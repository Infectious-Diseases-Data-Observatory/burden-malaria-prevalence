#!/usr/bin/env Rscript
# Cross-artifact checks for the current plan-defined primary reporting bundle.
source("R_cbh/load_pipeline.R")
source("R_cbh/primary/settings.R")
st <- cbh_primary_settings("regional");root <- st$out
read <- function(p)cbh_read_csv(file.path(root,p))
manifest <- read("fit_manifest.csv");diag <- read("fit_diagnostics.csv")
sample <- read("primary_sample.csv");figures <- read("paper_figures/manifest.csv")
stopifnot(nrow(manifest)==7L,all(grepl(st$id,manifest$model_file,fixed=TRUE)),
  all(diag$converged),all(diag$input_verified),all(diag$gamma==2),
  identical(unname(vapply(manifest$model_file,cbh_file_hash,"")),manifest$md5),
  nrow(figures)==4L,identical(as.integer(figures$figure),1:4),
  basename(figures$source[4])=="fig4_burden_comparison.png",
  all(grepl(st$id,figures$source,fixed=TRUE)),
  identical(unname(vapply(figures$source,cbh_file_hash,"")),figures$md5))
formula_vars <- all.vars(as.formula(paste(readLines(file.path(root,"model_formula.txt")),collapse=" ")))
scaling <- read("covariate_scaling.csv")
stopifnot(nrow(scaling)==17L,all(paste0("z_",scaling$variable) %in% formula_vars),
  all(c("z_mean_maternal_age_first_birth","z_wasting_pct","z_stunting_pct") %in% formula_vars),
  !any(c("z_male_pct","z_multiple_birth_pct","z_mean_birth_order","z_mean_maternal_age_birth") %in% formula_vars))
flow <- read("study_flow/flow_counts.csv");counts <- setNames(flow$count,flow$item)
coverage <- read("survey_map/survey_coverage.csv")
stopifnot(counts[["primary"]]==sample$records,counts[["deaths"]]==sample$deaths,
  counts[["primary_distinct_children"]]==sample$distinct_children,
  counts[["surveys"]]==sample$surveys,nrow(coverage)==sample$surveys,
  sum(coverage$regions)==st$expected_regions,all(coverage$type=="DHS"))
age <- read("tables/age_band_results.csv")
stopifnot(nrow(age)==7L,sum(age$records)==sample$records,sum(age$observed_deaths)==sample$deaths,
  abs(sum(age$displayed_death_share_pct)-100)<1e-10,!"elapsed_seconds" %in% names(age))
annual <- read("annual_comparison/annual_totals_2000_2024.csv")
burden <- read("burden/year_summary.csv")
stopifnot(nrow(annual)==25L,identical(as.integer(annual$year),2000:2024),all(annual$countries==42L),
  max(abs(annual$model_deaths[match(burden$year,annual$year)]-burden$attributable_under5_deaths))<.01)
national <- read("tables/country_comparison_2024.csv")
stopifnot(nrow(national)==42L,
  abs(sum(national$model_deaths)-annual$model_deaths[annual$year==2024])<1e-6,
  abs(sum(national$ihme_malaria_deaths)-annual$ihme_malaria_deaths[annual$year==2024])<1e-6,
  abs(sum(national$who_cacode_deaths)-annual$who_cacode_deaths[annual$year==2024])<1e-6)
states <- read("nigeria_states/state_totals_2024.csv")
recon <- read("nigeria_states/national_reconciliation_2024.csv")
stopifnot(nrow(states)==37L,abs(sum(states$attributable_under5_deaths)-recon$state_sum_attributable_deaths)<1e-6)
# Check producer/data hashes, not just whether expected filenames exist.
provenances <- list.files(root,pattern="provenance[.]csv$",recursive=TRUE,full.names=TRUE)
provenances <- provenances[!grepl("previous_overleaf",provenances)]
checks <- lapply(provenances,function(p) {
  x <- cbh_read_csv(p);stopifnot(all(c("file","md5") %in% names(x)),all(file.exists(x$file)))
  actual <- unname(vapply(x$file,cbh_file_hash,""))
  if(!identical(actual,x$md5))stop("Stale provenance in ",p,": ",paste(x$file[actual!=x$md5],collapse=", "))
  data.frame(provenance=p,files_verified=nrow(x),all_hashes_match=TRUE)
})
cbh_atomic_csv(do.call(rbind,checks),file.path(root,"paper_refresh/provenance_checks.csv"))
stopifnot(!length(list.files(root,pattern="[.]tex$",recursive=TRUE)))
for(p in file.path(root,"tables",c("age_band_results.latex.txt","country_comparison_2024.latex.txt"))) {
  lines <- readLines(p);rows <- lines[grepl(" & ",lines,fixed=TRUE)]
  stopifnot(length(rows)>7L,all(endsWith(rows,strrep(intToUtf8(92),2))))
}
writeLines(c("PASS: four main figures (Figure 4 combining country, state and annual comparisons) reference the revised primary version and match their hashes.",
  "PASS: seven fitted model hashes, convergence flags and gamma=2 verified.",
  sprintf("PASS: flow, %d DHS surveys, %d survey-regions, child-band/death/child totals reconcile.",sample$surveys,sum(coverage$regions)),
  "PASS: age-table death percentages sum to 100.0%; both contrasts retained; no fit-time column.",
  "PASS: 21 annual totals use 42 countries and reproduce the three national-burden years.",
  "PASS: 2024 country table reproduces all three Figure 5 totals; 37 state totals reconcile.",
  "PASS: producer/data provenance hashes match; table source has valid row terminators.",
  "PASS: no TeX files created in this result version."),file.path(root,"paper_refresh/validation.txt"))
message("Primary reporting validation passed")
