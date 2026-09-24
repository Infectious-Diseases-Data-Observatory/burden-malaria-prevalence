# SMC before/after analysis (requested 23 September 2026): the primary DHS+MICS
# age-band models (primary_map_regional17_dhsmics_gamma2_v5) refitted in the
# countries that have introduced seasonal malaria chemoprevention, with one extra
# term: an admin-1 indicator equal to 1 when the band entry year is at or after the
# admin-1 unit's SMC switch-on year. Everything else is the primary model: formula,
# 17 covariates with the full-sample scaling, reference knots, gamma, random effects.
#
# Source: data/SMC_rollout/smc_by_dhs_cluster.csv (GiveWell DataWell, DHS/MIS
# clusters in 10 countries; licensed DHS cluster identifiers, so cluster-level
# data never leave the ignored data folder). smc_first_year is fixed per area.
# Admin-1 unit: the analysis region, except Nigeria, whose analysis regions are the
# six zones; there the unit is the state (from the recode's state variable by
# cluster). Unit switch-on year: pooled over every round's clusters, the first year
# by which at least half of the unit's clusters had SMC (unmatched clusters dropped;
# 'none' = never).
#
# Status basis (codebook): district- or region-level records confirm a switch-on;
# national-scope records only say the country ran SMC that year ("not a treatment
# flag"). Main analysis (decided 23 September 2026): confirmed status only.
#   main      binary indicator from confirmed switch-ons; records of national-basis
#             units are dropped from the unit's first recorded campaign year onward
#   no_smc    the main sample without the SMC term (reference for the PfPR curves)
#   placebo   main sample, pre-switch-on records only, fake switch-on 3 years early
#   national  sensitivity: every SMC country, national-scope years taken as switch-on
#   coverage  sensitivity: main sample, share of the unit's confirmed clusters with SMC
#             by the entry year (0-1) in place of the binary indicator
cbh_smc_settings <- function() list(
  # v2 (24 September 2026) refits on the v7 primary (with Liberia); v1 was fitted on v5.
  id = "smc_dhsmics_map_gamma2_v2",
  out = "results/cbh/smc_dhsmics_map_gamma2_v2",
  private = "data/derived_cbh/models/smc_dhsmics_map_gamma2_v2",
  primary_version = "regional_mics",
  source = "data/SMC_rollout/smc_by_dhs_cluster.csv",
  countries = c(burkina_faso = "BFA", cameroon = "CMR", chad = "TCD", cote_divoire = "CIV", drc = "COD",
                mali = "MLI", niger = "NER", nigeria = "NGA", uganda = "UGA", zambia = "ZMB"),
  majority = 0.5,
  placebo_lead = 3L,
  aliases = data.frame(
    iso3 = c("CMR", "CMR", "CMR", "TCD", "NGA", "NGA"),
    source_key = c("adamaoua", "centrewithoutyaounde", "littoralwithoutdouala", "ennedi", "fctabuja", "fct"),
    regkey = c("adamawa", "centre", "littoral", "ennediestennediouest", "federalcapitalterritory", "federalcapitalterritory")),
  # Nigerian surveys whose band entry years reach the first SMC year (2013): state by cluster.
  nigeria_state_sources = data.frame(
    survey = c("NG6AFL", "NG7BFL", "NG8BFL", "MC_NGA2016"),
    file = c("registry", "registry", "registry", "data/MICS_extracted/NGA_2016_MICS5_v01_M/hh.sav"),
    cluster = c("v001", "v001", "v001", "HH1"), state = c("sstate", "sstate", "sstate", "HH7")),
  # v5 countries with SMC programmes that are not in the cluster file (the codebook names
  # Mozambique; the others are to be confirmed with the data provider).
  omitted_smc_countries = c("SEN", "GMB", "GIN", "GNB", "GHA", "TGO", "BEN", "MRT", "MOZ"),
  variants = c("main", "no_smc", "placebo", "national", "coverage"),
  labels = c(main = "SMC indicator, confirmed status (main)", no_smc = "No SMC term (main sample)",
             placebo = "Placebo: switch-on 3 years early", national = "SMC indicator, national years as recorded",
             coverage = "SMC coverage share, confirmed (0-1)"))
