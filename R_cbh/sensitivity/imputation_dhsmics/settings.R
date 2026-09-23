# Imputed-covariate sensitivity on the DHS and MICS sample (v6).
# The DHS-only version (primary_map_regional17_imputed_gamma2_v4) runs through the
# primary scripts; this one is self-contained so that the committed DHS+MICS primary
# (primary_map_regional17_dhsmics_gamma2_v5), whose code provenance covers
# R_cbh/primary/settings.R and 00_prepare_regional.R, stays untouched.
# Same 17 covariates, formula, reference knots and gamma as v5; every MAP-eligible
# DHS and MICS record is retained after imputing all remaining covariate gaps.
cbh_imputed_dhsmics_settings <- function() {
  source("R_cbh/primary/settings.R", local = TRUE)
  base <- cbh_primary_settings("regional_mics")
  st <- base[c("knots", "years", "gamma", "seed", "nthreads", "regional", "nutrition",
               "mics_output_dir", "mics_overlay")]
  st$id <- "primary_map_regional17_dhsmics_imputed_gamma2_v6"
  st$out <- file.path("results/cbh", st$id)
  st$private <- file.path("data/derived_cbh/models", st$id)
  st$data <- file.path(st$private, "complete_case_dataset.rds")
  st$reference <- base$out
  st$reference_data <- base$data
  st$imputed <- TRUE; st$mics <- TRUE
  st$imputation_private <- "data/derived_cbh/regional_adjustment/imputed_dhsmics_v6"
  st$imputation_results <- "results/cbh/covariate_imputation_dhsmics_v6"
  st$imputed_overlay <- file.path(st$imputation_private, "regional_covariates_wide.csv")
  st$imputed_national <- file.path(st$imputation_private, "national_covariates_imputed.csv")
  st$dhs_overlay <- "data/derived_cbh/regional_adjustment/planned17_audit/regional_covariates_wide.csv"
  st$dhs_manifest <- "data/derived_cbh/survey_manifest.csv"
  st$mics_manifest <- file.path(st$mics_output_dir, "survey_manifest.csv")
  st$registries <- c("data/derived_dhs/survey_registry.csv", "data/derived_mics/survey_registry_mics.csv")
  st$annual_maps <- c("data/derived_dhs/map_pfpr_window_years.csv", "data/derived_mics/map_pfpr_window_years_mics.csv")
  st$hiv_panel <- "data/derived_cbh/hiv_incidence/child_incidence_country_year_extended.csv"
  st$hiv_draws <- "data/derived_cbh/hiv_incidence/child_incidence_draws_extended.rds"
  st$comparison_labels <- c("Complete-case covariates, DHS and MICS (v5)", "Imputed covariates, DHS and MICS (v6)")
  # Every MAP-eligible record (model_ready totals of the two survey manifests). The DHS
  # part equals the DHS-only imputed version v4.
  st$expected_dhs_records <- 6357802L; st$expected_dhs_deaths <- 90938L; st$expected_dhs_regions <- 1113L
  st$expected_records <- 8797963L; st$expected_deaths <- 123419L; st$expected_regions <- 1457L
  st$expected_surveys <- 166L; st$expected_countries <- 40L
  st$m <- 10L; st$maxit <- 20L; st$imputation_seed <- 20260918L
  st
}
