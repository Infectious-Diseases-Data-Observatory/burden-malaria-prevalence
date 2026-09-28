#!/usr/bin/env Rscript
# Local-only revised covariate pipeline. Public DHS retrieval is separate.
args <- commandArgs(TRUE)
stopifnot(all(args %in% c("--complete-case-only","--legacy-expanded")))
# 17 replaces excluded published regional rows (published_survey_exclusions.csv) before 03.
scripts <- c("02_extract_regional.R","17_recode_replacements.R",if(!"--complete-case-only" %in% args)"07_prepare_unicef.R",
  "03_audit_missingness.R","04_report.R","05_validate.R")
for(script in scripts) {
  flags <- if(script %in% c("03_audit_missingness.R","04_report.R","05_validate.R")) args else character()
  status <- system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path("R_cbh/covariates",script)),flags))
  if(status!=0L)stop("Covariate stage failed: ",script)
}
