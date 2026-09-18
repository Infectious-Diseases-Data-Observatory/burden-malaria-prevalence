# Supplementary results

Reporting designation updated 18 September 2026. The primary analysis uses **full-sample MAP prevalence, seven separate age-band models and gamma=2**; see the [analysis plan](../docs/ANALYSIS_PLAN.md). The analyses below are supplementary sensitivity checks of that primary model.

## Removed analyses

The comparison of MAP with the annual Snow prevalence estimates (matched MAP–Snow fits, Snow penalisation trials and Snow-based national mortality scenarios) was removed from the analysis plan on 18 September 2026 and is not part of the manuscript. Its code is in [`archive/2026-09-18-snow-comparison/`](../archive/2026-09-18-snow-comparison/README.md). The saved outputs in `results/cbh/map_snow_gamma2_v1/` and `results/cbh/age_band_snow_2000_2015_v1/` are retained unchanged as historical records; `results/cbh/map_snow_gamma2_v1/` also holds the original primary `map_full` fits and is still read by the saved-primary audit.

## Sahelian child mortality by calendar month

Added 16 September 2026. A descriptive analysis of the all-cause death rate of children under two years by month of the year in Sahelian survey regions (boundary centroid at or above 11°N, west of the Horn of Africa; 247 regions in 41 surveys), by age band (<1, 1–5, 6–11, 12–23 completed months). Month of death is derived from month of birth and reported age at death; see [R_cbh/seasonality/README.md](../R_cbh/seasonality/README.md) for the derivation and its limits. It fits no prevalence model.

- [Figure: pooled rates by month and the negative-binomial month effects](../results/cbh/seasonality_sahel_v1/sahel_mortality_by_calendar_month.png)
- [Report with tables](../results/cbh/seasonality_sahel_v1/REPORT.md)
- [Region list with centroids and mean MAP prevalence](../results/cbh/seasonality_sahel_v1/sahel_regions.csv)

## Sahel-only PfPR sensitivity at gamma=2

Added 16 September 2026. Refit the seven primary age-band models in survey regions with centroids ≥12°N under the existing western/Horn exclusion rule. This uses MAP exposure, the primary adjustment set and fixed HIV incidence imputation. The subset contains 378,680 children and 15,999 deaths from 25 surveys in seven countries. Whole regions are selected by centroid; the cutoff is an approximate seasonal-area proxy. In particular, Nigeria's broad northern zones fall below the cutoff. This is separate from the 11°N descriptive monthly-rate analysis above.

- [Hazard-ratio overlay against the full primary analysis](../results/cbh/sahel_map_gamma2_v1/sahel_vs_full_pfpr_hazard_ratios.png)
- [Log-hazard-ratio version](../results/cbh/sahel_map_gamma2_v1/sahel_vs_full_pfpr_log_hazard_ratios.png)
- [Report, model details and contrasts](../results/cbh/sahel_map_gamma2_v1/REPORT.md)
- [Region selection audit](../results/cbh/sahel_map_gamma2_v1/region_selection.csv)
- [Reproduction instructions](../R_cbh/sensitivity/sahel/README.md)
