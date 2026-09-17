#!/usr/bin/env Rscript
# Run from the project root. Data setup and supplementary fits are excluded.
source("R_cbh/primary/settings.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>1L || (length(args) && !args %in% c("--report-only","--resume")))
  stop("Use no arguments (fresh fit), --resume (valid fit caches), or --report-only (saved fits).")
settings <- cbh_primary_settings()
stages <- c(fit="R_cbh/primary/01_fit.R",effects="R_cbh/primary/02_effects.R",
  diagnostics="R_cbh/primary/04_diagnostics.R",report="R_cbh/primary/03_report.R",
  flow="R_cbh/reporting/01_study_flow.R",survey_map="R_cbh/reporting/03_survey_map.R",
  age_fractions="R_cbh/reporting/06_attributable_fraction_by_age.R",
  annual_burden="R_cbh/burden/04_annual_comparison.R",
  annual_comparison="R_cbh/reporting/07_annual_mortality_comparison.R",
  paper_figures="R_cbh/reporting/04_paper_figures.R",index="R_cbh/reporting/02_results_index.R")
if("--report-only" %in% args) stages <- stages[names(stages)!="fit"]
dir.create(settings$out,recursive=TRUE,showWarnings=FALSE)
log <- list()
for(stage in names(stages)) {
  started <- Sys.time()
  message("Primary pipeline stage: ",stage)
  trailing <- if(stage=="fit" && !"--resume" %in% args) "--force" else character()
  status <- system2(file.path(R.home("bin"),"Rscript"),c(stages[[stage]],trailing))
  log[[stage]] <- data.frame(stage=stage,started_at=as.character(started),finished_at=as.character(Sys.time()),
    elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),exit_status=status)
  write.csv(do.call(rbind,log),file.path(settings$out,"pipeline_stages.csv"),row.names=FALSE)
  if(status!=0L) stop("Primary stage failed: ",stage,". Outputs from earlier stages are retained for inspection.")
}
files <- list.files(settings$out,recursive=TRUE,full.names=TRUE)
files <- files[!grepl("output_manifest[.]csv$",files)]
write.csv(data.frame(file=files,bytes=file.info(files)$size,md5=unname(tools::md5sum(files))),
  file.path(settings$out,"output_manifest.csv"),row.names=FALSE)
message("Primary pipeline complete: ",settings$out,"/REPORT.md")
