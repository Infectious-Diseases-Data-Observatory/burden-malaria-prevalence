#!/usr/bin/env Rscript
# Current 18-variable primary fit and comparison; no TeX or source extraction.
source("R_cbh/primary/settings.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>1L || (length(args) && !args %in% c("--resume","--report-only")))
  stop("Use no arguments (fresh fits), --resume (verified caches), or --report-only (saved fits).")
Sys.setenv(CBH_PRIMARY_VERSION="regional")
settings <- cbh_primary_settings()
stages <- c(prepare="R_cbh/primary/00_prepare_regional.R",fit="R_cbh/primary/01_fit.R",
  effects="R_cbh/primary/02_effects.R",diagnostics="R_cbh/primary/04_diagnostics.R",
  comparison="R_cbh/primary/05_compare_regional.R")
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
message("Revised primary pipeline complete: ",settings$out,"/REPORT.md")
