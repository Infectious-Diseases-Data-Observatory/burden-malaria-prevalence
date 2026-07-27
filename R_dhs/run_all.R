# =============================================================================
# run_all.R
# Orchestrate the rebuilt DHS/MIS-only analysis.
#
# Full rebuild from existing authorised recodes and cached/downloadable MAP:
#   Rscript R_dhs/run_all.R
#
# Skip data-access and MAP extraction when their outputs already exist:
#   Rscript R_dhs/run_all.R --analysis-only
#
# Validate models against the pre-existing aggregate panel without reading
# individual DHS records:
#   Rscript R_dhs/run_all.R --from-legacy-aggregate
# =============================================================================

source("R_dhs/00_config.R")

args <- commandArgs(trailingOnly = TRUE)
legacy_mode <- "--from-legacy-aggregate" %in% args
analysis_only <- "--analysis-only" %in% args

run_script <- function(script, trailing = character(0)) {
  path <- file.path("R_dhs", script)
  message("\n=== Running ", path, " ", paste(trailing, collapse = " "), " ===")
  status <- system2("Rscript", c(path, trailing))
  if (!identical(status, 0L)) stop(path, " failed with status ", status)
}

if (legacy_mode) {
  run_script("03_build_analysis_dataset.R", "--from-legacy-aggregate")
} else {
  if (!analysis_only) {
    run_script("01_access_dhs_data.R")
    run_script("02_extract_map_pfpr.R")
  }
  run_script("03_build_analysis_dataset.R")
}

run_script("04_fit_main_models.R")
run_script("05_make_main_plots.R")
run_script("06_sensitivity_samples.R")
run_script("07_sensitivity_model_structure.R")
run_script("08_sensitivity_likelihood.R")

legacy_validation_inputs <- c(
  file.path(REPO_ROOT, "results", "penalized_models.rds"),
  file.path(REPO_ROOT, "results", "paper_fig2_negcontrol.csv")
)
if (all(file.exists(legacy_validation_inputs))) {
  run_script("09_validate_reproduction.R")
} else {
  message(
    "\nLegacy fitted-model artifacts are absent; skipping the optional ",
    "migration validation."
  )
}

message("\nRebuilt DHS/MIS analysis complete. Outputs: ", RESULTS_DIR)
