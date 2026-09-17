# Current primary specification. Prepared data and national inputs are held fixed.
cbh_primary_settings <- function(version=Sys.getenv("CBH_PRIMARY_VERSION","legacy")) {
  stopifnot(version %in% c("legacy","regional"))
  settings <- list(
  id = "primary_map_gamma2_v1",
  out = "results/cbh/primary_map_gamma2_v1",
  private = "data/derived_cbh/models/primary_map_gamma2_v1",
  data = "data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds",
  reference = "results/cbh/map_snow_gamma2_v1",
  knots = "R_cbh/primary/reference_knots.csv",
  years = c(2005L, 2015L, 2024L), gamma = 2, seed = 20260907L,
  nthreads = 2L,regional=FALSE)
  if(version=="regional") {
    settings$id <- "primary_map_regional18_gamma2_v2"
    settings$out <- file.path("results/cbh",settings$id)
    settings$private <- file.path("data/derived_cbh/models",settings$id)
    settings$data <- file.path(settings$private,"complete_case_dataset.rds")
    settings$reference <- "results/cbh/primary_map_gamma2_v1"
    settings$regional <- TRUE
  }
  settings
}
