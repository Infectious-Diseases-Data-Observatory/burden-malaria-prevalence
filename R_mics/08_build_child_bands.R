#!/usr/bin/env Rscript
# Run the existing child age-band builder (R_cbh/R/build.R) on the MICS surveys, writing to
# data/derived_mics/cbh_build so the DHS manifest and shards are untouched. Paths are
# overridden in memory; no shared builder file is edited.
source("R_cbh/load_pipeline.R")
cfg <- cbh_config()
d <- file.path(cfg$root, "data/derived_mics")
cfg$output_dir <- file.path(d, "cbh_build")
cfg$registry <- file.path(d, "survey_registry_mics.csv")
cfg$survey_rules <- file.path(d, "survey_rules_mics.csv")
cfg$region_overrides <- file.path(d, "region_overrides_mics.csv")
cfg$boundary_regions <- file.path(d, "map_pfpr_by_survey_region_mics.csv")
cfg$annual_map <- file.path(d, "map_pfpr_window_years_mics.csv")
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
res <- cbh_build(cfg, surveys = NULL, force = "--force" %in% args)
m <- readRDS(file.path(cfg$output_dir, "manifest.rds"))
cat("MICS build complete:", sum(m$status %in% c("built", "cached")), "built of", nrow(m), "\n")
print(table(m$status))
