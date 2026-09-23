# Subgroup refits of the current primary PfPR models (plan section 3.4.1).
# Four subsets of the fitted 17-variable sample, each refitted separately by age band:
#   sahel        survey regions with boundary centroid >= 12N, < 36E, outside ETH/ERI/SOM/DJI
#   east_africa  UN M49 Eastern Africa countries
#   early        surveys with survey year <= median survey year (one vote per survey)
#   late         surveys with survey year >  median survey year
# sample = "dhsmics" (default) refits the DHS and MICS primary (v5); sample = "dhs"
# reproduces the earlier DHS-only refits of v3 (subgroups_map_gamma2_v1).
cbh_subgroup_settings <- function(sample = c("dhsmics", "dhs")) {
  sample <- match.arg(sample)
  st <- list(
    sample = sample,
    id = "subgroups_map_gamma2_v1",
    primary_version = "regional",
    centroids = "data/derived_dhs/dhs_region_centroids.csv",
    registry = "data/derived_dhs/survey_registry.csv",
    sahel = list(latitude_min = 12, longitude_max = 36, exclude_countries = c("ETH", "ERI", "SOM", "DJI")),
    # UN M49 Eastern Africa (https://unstats.un.org/unsd/methodology/m49/), restricted to
    # countries that can appear in the DHS/MIS/MICS sample.
    east_africa = c("BDI", "COM", "DJI", "ERI", "ETH", "KEN", "MDG", "MOZ", "MWI", "RWA",
                    "SOM", "SSD", "TZA", "UGA", "ZMB", "ZWE"),
    labels = c(map_full = "Full analysis", sahel = "Sahel (≥12°N)",
               east_africa = "Eastern Africa", early = "Surveys up to median year",
               late = "Surveys after median year"),
    colours = c(map_full = "#222222", sahel = "#C34D26", east_africa = "#009E73",
                early = "#7B4FA3", late = "#E69F00"))
  if (sample == "dhsmics") {
    st$id <- "subgroups_dhsmics_map_gamma2_v2"
    st$primary_version <- "regional_mics"
    st$centroids <- c(st$centroids, "data/derived_mics/mics_region_centroids.csv")
    st$registry <- c(st$registry, "data/derived_mics/survey_registry_mics.csv")
  }
  st$out <- file.path("results/cbh", st$id)
  st$private <- file.path("data/derived_cbh/models", st$id)
  st
}
cbh_subgroup_sample <- function(args) if ("--dhs-only" %in% args) "dhs" else "dhsmics"
