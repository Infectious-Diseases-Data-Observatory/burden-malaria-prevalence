# Supplementary results

Reporting designation updated 18 September 2026. The primary analysis uses **full-sample MAP prevalence, seven separate age-band models and gamma=2**; see the [analysis plan](../docs/ANALYSIS_PLAN.md). The analyses below are supplementary sensitivity checks of that primary model.

## Supplementary figures for the manuscript

- [Count version of main Figure 4](../results/cbh/primary_map_regional17_gamma2_v3/burden_comparison/sfig_burden_comparison_counts.png): the 2024 country and Nigerian-state comparisons with IHME and the 2000–2024 annual trend as numbers of deaths (log axes), with its [caption](../results/cbh/primary_map_regional17_gamma2_v3/burden_comparison/CAPTION_counts.md). The main figure shows the same comparisons per 1,000 under-five child-years.
- [Inclusion flow](../results/cbh/primary_map_regional17_gamma2_v3/study_flow/study_flow_diagram.png) and [caption](../results/cbh/primary_map_regional17_gamma2_v3/study_flow/CAPTION.md).
- [PfPR splines by subgroup](../results/cbh/subgroups_map_gamma2_v1/sfig_pfpr_splines_by_subgroup.png) and [caption](../results/cbh/subgroups_map_gamma2_v1/CAPTION.md): the primary curves overlaid with refits in the Sahel, Eastern Africa and surveys before/after the median survey year (see below).

## Imputed-covariate sensitivity version

Added 18 September 2026. The 17-variable primary specification refitted after imputing every remaining covariate gap (whole-survey gaps in wasting, stunting, facility delivery, electricity and wealth by chained-equation multiple imputation; the 2001 WGI round by interpolation; health expenditure for Zimbabwe 2000–2009 and all countries' 2024 from a GAM; child HIV incidence for Liberia and São Tomé from the extended incidence model), so that all 6,357,802 MAP-eligible records from 120 surveys in 36 countries are retained. The complete-case fit remains the primary; this is `primary_map_regional17_imputed_gamma2_v4`. See plan sections 2.4 and 3.4.3.

- [Comparison with the complete-case primary](../results/cbh/primary_map_regional17_imputed_gamma2_v4/REPORT.md)
- [Multiple-imputation check: pooled curves and contrasts](../results/cbh/primary_map_regional17_imputed_gamma2_v4/multiple_imputation/REPORT.md)
- [Supplementary manuscript figure: complete-case versus imputed-covariate PfPR curves](../results/cbh/primary_map_regional17_imputed_gamma2_v4/sfig_pfpr_splines_imputed_covariates.png) and [caption](../results/cbh/primary_map_regional17_imputed_gamma2_v4/CAPTION_sfig_imputed_covariates.md), exported to the manuscript's `Supplementary Figures/` folder
- [Imputation audits: regional](../results/cbh/covariate_imputation_v4/REGIONAL_IMPUTATION.md) and [national](../results/cbh/covariate_imputation_v4/NATIONAL_IMPUTATION.md)
- [Extended HIV incidence panel](../results/cbh/hiv_incidence/imputed_no_adolescent_series.csv) and [latent-adolescent validation](../results/cbh/hiv_incidence/latent_adolescent_validation.csv)

## Subgroup refits of the primary model at gamma=2

Added 18 September 2026. The seven 17-variable primary age-band models refitted separately in four subsets of the fitted sample: Sahel survey regions (centroid ≥12°N, west of 36°E, outside the Horn of Africa; 23 surveys in seven countries), UN M49 Eastern Africa (38 surveys in 12 countries), surveys up to the median survey year of 2013 (50 surveys) and surveys after it (45 surveys). Same covariates and scaling, fixed HIV imputation, gamma=2 and basis dimensions; knots at the subset's own predictor quantiles. Whole regions and whole surveys enter. The subsets overlap the full sample, so differences are descriptive.

- [Overlay figure](../results/cbh/subgroups_map_gamma2_v1/sfig_pfpr_splines_by_subgroup.png) and [40%→20% contrasts by subgroup](../results/cbh/subgroups_map_gamma2_v1/pfpr_40_to_20_by_subgroup.png)
- [Report with contrast table and subset sizes](../results/cbh/subgroups_map_gamma2_v1/REPORT.md)
- [Region and survey selection audits](../results/cbh/subgroups_map_gamma2_v1/region_selection.csv)
- [Reproduction instructions](../R_cbh/sensitivity/README.md)

## Removed analyses

The comparison of MAP with the annual Snow prevalence estimates (matched MAP–Snow fits, Snow penalisation trials and Snow-based national mortality scenarios) was removed from the analysis plan on 18 September 2026 and is not part of the manuscript. Its code is in [`archive/2026-09-18-snow-comparison/`](../archive/2026-09-18-snow-comparison/README.md). The saved outputs in `results/cbh/map_snow_gamma2_v1/` and `results/cbh/age_band_snow_2000_2015_v1/` are retained unchanged as historical records; `results/cbh/map_snow_gamma2_v1/` also holds the original primary `map_full` fits and is still read by the saved-primary audit.

## Sahelian child mortality by calendar month

Added 16 September 2026. A descriptive analysis of the all-cause death rate of children under two years by month of the year in Sahelian survey regions (boundary centroid at or above 11°N, west of the Horn of Africa; 247 regions in 41 surveys), by age band (<1, 1–5, 6–11, 12–23 completed months). Month of death is derived from month of birth and reported age at death; see [R_cbh/seasonality/README.md](../R_cbh/seasonality/README.md) for the derivation and its limits. It fits no prevalence model.

- [Figure: pooled rates by month and the negative-binomial month effects](../results/cbh/seasonality_sahel_v1/sahel_mortality_by_calendar_month.png)
- [Report with tables](../results/cbh/seasonality_sahel_v1/REPORT.md)
- [Region list with centroids and mean MAP prevalence](../results/cbh/seasonality_sahel_v1/sahel_regions.csv)

## Sahel-only PfPR sensitivity at gamma=2 (historical, 11-variable fit)

Added 16 September 2026; superseded by the subgroup refits above, which repeat the Sahel restriction on the current 17-variable sample. Refit the seven primary age-band models in survey regions with centroids ≥12°N under the existing western/Horn exclusion rule. This uses MAP exposure, the primary adjustment set and fixed HIV incidence imputation. The subset contains 378,680 children and 15,999 deaths from 25 surveys in seven countries. Whole regions are selected by centroid; the cutoff is an approximate seasonal-area proxy. In particular, Nigeria's broad northern zones fall below the cutoff. This is separate from the 11°N descriptive monthly-rate analysis above.

- [Hazard-ratio overlay against the full primary analysis](../results/cbh/sahel_map_gamma2_v1/sahel_vs_full_pfpr_hazard_ratios.png)
- [Log-hazard-ratio version](../results/cbh/sahel_map_gamma2_v1/sahel_vs_full_pfpr_log_hazard_ratios.png)
- [Report, model details and contrasts](../results/cbh/sahel_map_gamma2_v1/REPORT.md)
- [Region selection audit](../results/cbh/sahel_map_gamma2_v1/region_selection.csv)
- [Reproduction instructions](../R_cbh/sensitivity/sahel/README.md)

## Adjustment without wasting and stunting

Added 22 September 2026. The seven primary age-band models refitted on the identical complete-case sample with the two nutritional covariates removed (15 covariates instead of 17). Wasting and stunting are measured in surviving children at the time of the survey and may lie on the causal pathway from malaria to death, so adjusting for them risks removing part of the association of interest. Only the adjustment set changes: sample, exposure, reference knots, basis dimensions, gamma=2 and random effects are those of the primary. See plan section 3.4.4.

- [Supplementary figure](../results/cbh/nutrition_adjustment_map_gamma2_v1/sfig_pfpr_splines_without_nutrition.png) and [caption](../results/cbh/nutrition_adjustment_map_gamma2_v1/CAPTION.md)
- [Hazard-ratio comparison](../results/cbh/nutrition_adjustment_map_gamma2_v1/CONTRASTS.md) and [contrast table](../results/cbh/nutrition_adjustment_map_gamma2_v1/comparison_contrasts.csv)
- [Fit diagnostics](../results/cbh/nutrition_adjustment_map_gamma2_v1/fit_diagnostics.csv)
