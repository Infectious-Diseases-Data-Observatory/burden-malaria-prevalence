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
# Validate models against the pre-existing aggregate panel, reading local DHS
# recodes only to reconstruct the three direct regional vaccine measures:
#   Rscript R_dhs/run_all.R --from-legacy-aggregate
# =============================================================================

source("R_dhs/00_config.R")

args <- commandArgs(trailingOnly = TRUE)
legacy_mode <- "--from-legacy-aggregate" %in% args
analysis_only <- "--analysis-only" %in% args

run_script <- function(script, trailing = character(0), fatal = TRUE) {
  path <- file.path("R_dhs", script)
  message("\n=== Running ", path, " ", paste(trailing, collapse = " "), " ===")
  status <- system2("Rscript", c(path, trailing))
  if (identical(status, 0L)) return(invisible(TRUE))
  if (fatal) stop(path, " failed with status ", status)
  message(
    "\n", path, " reported differences (status ", status, "). Continuing: ",
    "this is a comparison against a superseded baseline, not a gate on the ",
    "current build. Review its table before relying on the legacy numbers."
  )
  invisible(FALSE)
}

if (file.exists(UNICEF_GLOBAL_CSV)) {
  panel_outdated <- !file.exists(UNICEF_IMMUNISATION_CSV) ||
    file.info(UNICEF_GLOBAL_CSV)$mtime >
      file.info(UNICEF_IMMUNISATION_CSV)$mtime
  if (panel_outdated) {
    run_script("02b_extract_unicef_immunisation.R")
  }
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
run_script("15_specification_forest.R")
run_script("16_sensitivity_random_effects.R")
# Timing sensitivity: needs the cached annual MAP rasters, DHS boundaries and raw
# Births Recodes to build its two intermediates; both are cached per survey and
# the script exits cleanly when those inputs are unavailable.
run_script("17_sensitivity_timing.R")
run_script("18_age_windows_5to14.R")
# Period-stratified refits; needs the intermediates built by script 17.
run_script("19_period_stratified.R")
run_script("20_subregion_stratified.R")
run_script("21_map_vs_measured_prevalence.R")
run_script("22_seasonality_of_fieldwork.R")
run_script("23_deprivation_proxy.R")
run_script("24_intervention_targeting.R")
run_script("25_national_burden_lagged.R")
run_script("26_audit_survey_coverage.R")
# Mortality horizon x prevalence lag. Needs the lagged PfPR panel that script 17
# builds, so it runs after it.
run_script("27_horizon_lag_selection.R")
run_script("28_survey_map.R")
run_script("29_subgroup_fits.R")
# Bayesian refit of the tensor surface. Caches its sampled fit and refits only
# when the panel changes, so this costs about two minutes on a warm cache.
run_script("30_brms_tensor_model.R")
run_script("31_brms_subgroup_fits.R")
run_script("32_brms_three_way.R")

legacy_validation_inputs <- c(
  file.path(REPO_ROOT, "results", "penalized_models.rds"),
  file.path(REPO_ROOT, "results", "paper_fig2_negcontrol.csv")
)
if (all(file.exists(legacy_validation_inputs))) {
  # Informational only. Script 09 compares the rebuild against the legacy
  # aggregate's fitted models, and the raw rebuild deliberately analyses a
  # larger panel, so exceeding its tolerances is expected rather than a fault.
  # Run it directly for a hard pass/fail gate.
  run_script("09_validate_reproduction.R", fatal = FALSE)
} else {
  message(
    "\nLegacy fitted-model artifacts are absent; skipping the optional ",
    "migration validation."
  )
}

# National burden / time-surface comparison. It is downstream of the DHS/MIS
# model rebuild but produces all of the burden, temporal-trend and IHME/WHO
# numbers, so it is run here whenever its (local) inputs are present. Guarded
# like the script-09 migration check so a partial checkout skips gracefully.
burden_inputs <- c(
  file.path(DATA_DIR, "wb_mortality_timeseries.csv"),
  file.path(DATA_DIR, "pfpr_by_country_year.csv"),
  file.path(DATA_DIR, "pfpr_by_country_2024.csv"),
  file.path(DATA_DIR, "igme_mortality_by_country.csv"),
  file.path(DATA_DIR, "wb_livebirths_by_country.csv"),
  file.path(DATA_DIR, "ihme_malaria_u5_deaths_by_country.csv"),
  file.path(REPO_ROOT, "results", "ihme_u5_deaths_ssa_timeseries.csv"),
  file.path(REPO_ROOT, "results", "who_wmr2025_africa_deaths.csv")
)
if (all(file.exists(burden_inputs))) {
  run_script("10_time_surface_and_burden.R")
} else {
  message(
    "\nNational-burden inputs are absent; skipping script 10. Missing: ",
    paste(basename(burden_inputs[!file.exists(burden_inputs)]), collapse = ", ")
  )
}

# Manuscript figures kept in R_dhs: the study-flow diagram (11, self-contained),
# the RCT triangulation of the prevalence->mortality link (12, uses the fitted
# model bundle), and the Method-1 IHME attributable-share comparison (13, which
# skips itself if the Method-1 country panel from R/02 is absent).
run_script("11_study_flow.R")
if (file.exists(MODEL_BUNDLE_RDS)) {
  run_script("12_triangulation_rct.R")
} else {
  message("\nModel bundle absent; skipping script 12 (RCT triangulation).")
}
run_script("13_ihme_share.R")

message("\nRebuilt DHS/MIS analysis complete. Outputs: ", RESULTS_DIR)
