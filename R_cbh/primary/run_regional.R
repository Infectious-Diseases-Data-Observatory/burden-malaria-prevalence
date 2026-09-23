#!/usr/bin/env Rscript
# Current 17-variable primary fit and comparison; no TeX or source extraction.
source("R_cbh/primary/settings.R")
args <- commandArgs(trailingOnly=TRUE)
imputed <- "--imputed" %in% args; args <- setdiff(args,"--imputed")
mics <- "--mics" %in% args; args <- setdiff(args,"--mics")
dhs_only <- "--dhs-only" %in% args; args <- setdiff(args,"--dhs-only")
# Since 23 September 2026 the primary combines DHS and MICS surveys; --dhs-only runs the DHS-only v3.
if(!imputed && !dhs_only) mics <- TRUE
if(length(args)>1L || (length(args) && !args %in% c("--resume","--report-only")))
  stop("Use no arguments (fresh fits), --resume (verified caches), or --report-only (saved fits); add --imputed for the imputed-covariate sensitivity version.")
# --imputed fits the imputed-covariate sensitivity version (plan section 2.4) and
# runs only the model stages: it never touches the study flow, survey map, paper
# manifest, results index or validation, which describe the primary analysis.
# The DHS-plus-MICS version (primary_map_regional17_dhsmics_gamma2_v5) is the default primary and runs every
# stage; the study flow and survey map combine the DHS and MICS ledgers.
Sys.setenv(CBH_PRIMARY_VERSION=if(imputed) "regional_imputed" else if(mics) "regional_mics" else "regional")
settings <- cbh_primary_settings()
stages <- c(prepare="R_cbh/primary/00_prepare_regional.R",fit="R_cbh/primary/01_fit.R",
  effects="R_cbh/primary/02_effects.R",diagnostics="R_cbh/primary/04_diagnostics.R",
  comparison="R_cbh/primary/05_compare_regional.R",
  report="R_cbh/primary/03_report.R",flow="R_cbh/reporting/01_study_flow.R",
  survey_map="R_cbh/reporting/03_survey_map.R",age_fractions="R_cbh/reporting/06_attributable_fraction_by_age.R",
  annual_burden="R_cbh/burden/04_annual_comparison.R",annual_plot="R_cbh/reporting/07_annual_mortality_comparison.R",
  nigeria_burden="R_cbh/burden/05_nigeria_state_burden.R",nigeria_plot="R_cbh/reporting/08_nigeria_state_comparison.R",
  burden_figure="R_cbh/reporting/11_burden_comparison_figure.R",
  burden_tables="R_cbh/reporting/09_burden_tables.R",source_table="R_cbh/reporting/12_source_comparison_table.R",
  paper_manifest="R_cbh/reporting/04_paper_figures.R",index="R_cbh/reporting/02_results_index.R",
  validation="R_cbh/reporting/10_validate_reporting.R")
if(imputed) stages <- stages[names(stages) %in% c("prepare","fit","effects","diagnostics","comparison","report")]
if("--report-only" %in% args)stages <- stages[!names(stages) %in% c("prepare","fit")]
dir.create(settings$out,recursive=TRUE,showWarnings=FALSE)
log <- list()
for(stage in names(stages)) {
  started <- Sys.time();message("Revised primary stage: ",stage)
  trailing <- if(stage=="fit" && !"--resume" %in% args) "--force" else character()
  status <- system2(file.path(R.home("bin"),"Rscript"),c(stages[[stage]],trailing))
  log[[stage]] <- data.frame(stage=stage,started_at=as.character(started),finished_at=as.character(Sys.time()),
    elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),exit_status=status)
  write.csv(do.call(rbind,log),file.path(settings$out,"pipeline_stages.csv"),row.names=FALSE)
  if(status!=0L)stop("Revised primary stage failed: ",stage)
}
files <- list.files(settings$out,recursive=TRUE,full.names=TRUE)
files <- files[!grepl("output_manifest[.]csv$",files)]
write.csv(data.frame(file=files,bytes=file.info(files)$size,md5=unname(tools::md5sum(files))),
  file.path(settings$out,"output_manifest.csv"),row.names=FALSE)
message(if(imputed) "Imputed-covariate sensitivity pipeline complete: " else "Revised primary pipeline complete: ",settings$out,"/REPORT.md")
