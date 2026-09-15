# Current primary specification. Prepared data and national inputs are held fixed.
cbh_primary_settings <- function() list(
  id = "primary_map_gamma2_v1",
  out = "results/cbh/primary_map_gamma2_v1",
  private = "data/derived_cbh/models/primary_map_gamma2_v1",
  data = "data/derived_cbh/models/age_band_hiv_incidence_shared_time_v3/complete_case_dataset.rds",
  reference = "results/cbh/map_snow_gamma2_v1",
  knots = "R_cbh/primary/reference_knots.csv",
  years = c(2005L, 2015L, 2024L), gamma = 2, seed = 20260907L,
  nthreads = 2L)
