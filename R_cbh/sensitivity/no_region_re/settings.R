# Primary model without the survey-region random intercept, on the v7 primary sample.
cbh_no_region_settings <- function() {
  base <- cbh_primary_settings("regional_mics")
  st <- list(base = base, id = "no_region_re_dhsmics_map_gamma2_v2")
  st$out <- file.path("results/cbh", st$id)
  st$private <- file.path("data/derived_cbh/models", st$id)
  st
}
