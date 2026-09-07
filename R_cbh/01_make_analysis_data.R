#!/usr/bin/env Rscript
# Run from the project root. No network requests, installations or model fits.
source("R_cbh/load_pipeline.R")
args <- commandArgs(trailingOnly = TRUE)
allowed <- args %in% c("--check-inputs", "--force") | grepl("^--(surveys|output)=.+$", args)
if (any(!allowed)) stop("Arguments: --check-inputs --force --surveys=NG7BFL,ET8AFL --output=data/derived_cbh")
value <- function(key) {
  hit <- args[startsWith(args, paste0("--", key, "="))]
  if (length(hit) > 1L) stop("Repeated argument: ", key)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit) else NULL
}
cfg <- cbh_config()
if (!is.null(value("output"))) cfg$output_dir <- value("output")
surveys <- if (is.null(value("surveys"))) NULL else strsplit(value("surveys"), ",", fixed = TRUE)[[1]]
if ("--check-inputs" %in% args) {
  checks <- cbh_check_inputs(cfg, surveys)$checks
  print(checks, row.names = FALSE)
} else {
  m <- cbh_build(cfg, surveys, force = "--force" %in% args)
  message(sprintf("Finished: %s eligible band rows; %s with PfPR. Output: %s", sum(m$rows), sum(m$model_ready_rows), cfg$output_dir))
  if (any(m$status == "failed")) stop("Some surveys failed; inspect survey_manifest.csv.")
}
