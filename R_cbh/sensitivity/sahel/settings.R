# Sensitivity of the primary PfPR models; distinct from the 11N monthly-rate analysis.
cbh_sahel_settings <- function() list(
  out = "results/cbh/sahel_map_gamma2_v1",
  private = "data/derived_cbh/models/sahel_map_gamma2_v1",
  latitude_min = 12, longitude_max = 36,
  exclude_countries = c("ETH", "ERI", "SOM", "DJI"),
  centroids = "data/derived_dhs/dhs_region_centroids.csv",
  registry = "data/derived_dhs/survey_registry.csv")
