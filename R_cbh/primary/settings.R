# Current primary specification. Prepared data and national inputs are held fixed.
cbh_primary_settings <- function(version=Sys.getenv("CBH_PRIMARY_VERSION","legacy")) {
  stopifnot(version %in% c("legacy","regional","regional18","regional_imputed","regional_mics"))
  settings <- list(
  id = "primary_map_gamma2_v1",
  out = "results/cbh/primary_map_gamma2_v1",
  private = "data/derived_cbh/models/primary_map_gamma2_v1",
  data = "data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds",
  # Historical directory name: its map_full rows are the original primary fits.
  reference = "results/cbh/map_snow_gamma2_v1",
  knots = "R_cbh/primary/reference_knots.csv",
  years = c(2005L, 2015L, 2024L), gamma = 2, seed = 20260907L,
  nthreads = 2L,regional=FALSE,nutrition=FALSE,imputed=FALSE,mics=FALSE,
  hiv_panel = "data/derived_cbh/hiv_incidence/child_incidence_country_year.csv",
  comparison_labels = c("Previous 18-variable adjustment","Revised 17-variable adjustment"),
  expected_records=5885022L,expected_deaths=82415L,expected_regions=1015L)
  if(version %in% c("regional","regional18","regional_imputed","regional_mics")) {
    settings$id <- "primary_map_regional18_gamma2_v2"
    settings$out <- file.path("results/cbh",settings$id)
    settings$private <- file.path("data/derived_cbh/models",settings$id)
    settings$data <- file.path(settings$private,"complete_case_dataset.rds")
    settings$reference <- "results/cbh/primary_map_gamma2_v1"
    settings$regional <- TRUE
    settings$expected_records <- 5680117L;settings$expected_deaths <- 78634L;settings$expected_regions <- 973L
  }
  if(version %in% c("regional","regional_imputed","regional_mics")) {
    settings$id <- "primary_map_regional17_gamma2_v3"
    settings$out <- file.path("results/cbh",settings$id)
    settings$private <- file.path("data/derived_cbh/models",settings$id)
    settings$data <- file.path(settings$private,"complete_case_dataset.rds")
    settings$reference <- "results/cbh/primary_map_regional18_gamma2_v2"
    settings$nutrition <- TRUE
    settings$expected_records <- 5465305L;settings$expected_deaths <- 75726L;settings$expected_regions <- 916L
  }
  if(version=="regional_imputed") {
    # Sensitivity version: same 17 covariates, with the remaining covariate gaps
    # imputed (plan section 2.4) so every MAP-eligible record is retained.
    settings$id <- "primary_map_regional17_imputed_gamma2_v4"
    settings$out <- file.path("results/cbh",settings$id)
    settings$private <- file.path("data/derived_cbh/models",settings$id)
    settings$data <- file.path(settings$private,"complete_case_dataset.rds")
    settings$reference <- "results/cbh/primary_map_regional17_gamma2_v3"
    settings$imputed <- TRUE
    settings$imputed_overlay <- "data/derived_cbh/regional_adjustment/imputed_v4/regional_covariates_wide.csv"
    settings$imputed_national <- "data/derived_cbh/regional_adjustment/imputed_v4/national_covariates_imputed.csv"
    settings$imputation_results <- "results/cbh/covariate_imputation_v4"
    settings$hiv_panel <- "data/derived_cbh/hiv_incidence/child_incidence_country_year_extended.csv"
    settings$comparison_labels <- c("Complete-case covariates (v3)","Imputed covariates (v4)")
    # Every MAP-eligible record: totals from the survey manifest (model_ready rows/deaths).
    settings$expected_records <- 6357802L;settings$expected_deaths <- 90938L;settings$expected_regions <- 1113L
  }
  if(version=="regional_mics") {
    # DHS and MICS surveys with complete birth histories (decided 23 September 2026). The DHS part
    # is the complete-case sample of primary_map_regional17_gamma2_v3 and must reproduce it exactly;
    # MICS shards come from R_mics/08_build_child_bands.R and covariates from R_mics/09_regional_covariates.R.
    # Complete-case selection keeps only MICS surveys with anthropometry (the user's rule).
    settings$id <- "primary_map_regional17_dhsmics_gamma2_v5"
    settings$out <- file.path("results/cbh",settings$id)
    settings$private <- file.path("data/derived_cbh/models",settings$id)
    settings$data <- file.path(settings$private,"complete_case_dataset.rds")
    settings$reference <- "results/cbh/primary_map_regional17_gamma2_v3"
    settings$mics <- TRUE
    settings$mics_output_dir <- "data/derived_mics/cbh_build"
    settings$mics_overlay <- "data/derived_mics/regional_covariates_wide_mics.csv"
    settings$comparison_labels <- c("DHS surveys only (v3)","DHS and MICS surveys (v5)")
    settings$expected_dhs_records <- 5465305L;settings$expected_dhs_deaths <- 75726L;settings$expected_dhs_regions <- 916L
    # Combined DHS + MICS complete-case sample from the first preparation (23 September 2026).
    settings$expected_records <- 7498459L;settings$expected_deaths <- 102282L;settings$expected_regions <- 1211L
  }
  settings
}
