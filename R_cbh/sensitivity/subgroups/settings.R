# Subgroup refits of the current primary PfPR models (plan section 3.4.1).
# Four subsets of the fitted 17-variable sample, each refitted separately by age band:
#   sahel        survey regions with boundary centroid >= 12N, < 36E, outside ETH/ERI/SOM/DJI
#   east_africa  UN M49 Eastern Africa countries
#   early        surveys with survey year <= median survey year (one vote per survey)
#   late         surveys with survey year >  median survey year
cbh_subgroup_settings <- function() list(
  id = "subgroups_map_gamma2_v1",
  out = "results/cbh/subgroups_map_gamma2_v1",
  private = "data/derived_cbh/models/subgroups_map_gamma2_v1",
  centroids = "data/derived_dhs/dhs_region_centroids.csv",
  registry = "data/derived_dhs/survey_registry.csv",
  sahel = list(latitude_min = 12, longitude_max = 36, exclude_countries = c("ETH", "ERI", "SOM", "DJI")),
  # UN M49 Eastern Africa (https://unstats.un.org/unsd/methodology/m49/), restricted to
  # countries that can appear in the DHS/MIS sample.
  east_africa = c("BDI", "COM", "DJI", "ERI", "ETH", "KEN", "MDG", "MOZ", "MWI", "RWA",
                  "SOM", "SSD", "TZA", "UGA", "ZMB", "ZWE"),
  labels = c(map_full = "Full analysis", sahel = "Sahel (≥12°N)",
             east_africa = "Eastern Africa", early = "Surveys up to median year",
             late = "Surveys after median year"),
  colours = c(map_full = "#222222", sahel = "#C34D26", east_africa = "#009E73",
              early = "#7B4FA3", late = "#E69F00"))
