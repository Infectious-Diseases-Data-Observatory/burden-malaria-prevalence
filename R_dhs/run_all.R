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
# The survey-grouped 10-fold cross-validation (script 35) refits three brms
# models ten times each, about an hour on a 10-core machine, so it only runs
# when asked for.
with_kfold <- "--kfold" %in% args
# The Stan refit of the person-time models (script 43) takes about an hour.
with_person_time_brms <- "--person-time-brms" %in% args

run_script <- function(script, trailing = character(0), fatal = TRUE,
                       env = character(0)) {
  path <- file.path("R_dhs", script)
  message("\n=== Running ", paste(c(env, path, trailing), collapse = " "),
          " ===")
  status <- system2("Rscript", c(path, trailing), env = env)
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
# The same two programs refitted to neonatal mortality, the negative control.
# Their outputs carry a "_neonatal" suffix; see brms_outcome() in 00_config.R.
run_script("30_brms_tensor_model.R", env = "BRMS_OUTCOME=neonatal")
run_script("31_brms_subgroup_fits.R", env = "BRMS_OUTCOME=neonatal")
# The model ladder (additive / linear-in-time / full surface) with PSIS-LOO,
# the neonatal-covariate check and the burden trend with credible bands. All
# cache their fits and refit only when the panel changes.
run_script("34_brms_model_ladder.R")
if (with_kfold) run_script("35_brms_ladder_kfold.R")
run_script("36_neonatal_as_covariate.R")
run_script("37_burden_trend_brms.R")
# The person-time design: deaths and person-months by region, 12-month window
# and age segment, each window paired with the MAP prevalence of its own year,
# one negative-binomial model per age group. Script 40 reads every recode
# (about ten minutes cold, seconds warm), 41 re-extracts MAP over the window
# years, 42 fits.
run_script("40_build_person_time.R")
run_script("41_extract_map_window_years.R")
run_script("42_person_time_models.R")
if (with_person_time_brms) run_script("43_person_time_brms.R")
run_script("44_deaths_by_age.R")
run_script("45_person_time_data_and_curves.R")
run_script("46_person_time_covariate_forest.R")

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
# Main-text manuscript figures, drawn from the saved bundle and script 10's
# burden inputs.
run_script("14_manuscript_figures.R")
# Copy the curated key figures into "Key results/" with a README. Runs last so
# it always collects the figures this pass produced.
run_script("33_key_results.R")

message("\nRebuilt DHS/MIS analysis complete. Outputs: ", RESULTS_DIR)
