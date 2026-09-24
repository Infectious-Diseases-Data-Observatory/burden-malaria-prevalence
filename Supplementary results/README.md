# Supplementary results

Reporting designation updated 24 September 2026. The primary analysis uses **full-sample MAP prevalence, seven separate age-band models and gamma=2** on DHS and UNICEF MICS surveys, `primary_map_regional17_dhsmics_gamma2_v7` (135 surveys: 98 DHS and 37 MICS; 37 countries; 7,607,122 child-band records, 103,987 deaths); see the [analysis plan](../docs/ANALYSIS_PLAN.md). v7 is the 23 September primary `primary_map_regional17_dhsmics_gamma2_v5` plus Liberia, whose child HIV incidence is now derived from UNAIDS counts. The analyses below are supplementary sensitivity checks of that primary model. The subgroup (v3), no-nutrition (v3), imputed-covariate (v8) and SMC (v2) analyses are refitted on the v7 sample. Their versions on v5 are kept under [DHS and MICS v5 history](#dhs-and-mics-v5-history), and the earlier refits on the DHS-only primary (`primary_map_regional17_gamma2_v3`) under [DHS-only history](#dhs-only-history).

This page is maintained by hand; no script writes it.

## Supplementary figures for the manuscript

- [Count version of main Figure 3](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden_comparison/sfig_burden_comparison_counts.png): the 2024 country and Nigerian-state comparisons with IHME as numbers of deaths (log axes), and the 2000–2024 annual trend with IHME and UN IGME per 100,000 child-years, with its [caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden_comparison/CAPTION_counts.md). The main figure shows the same comparisons per 1,000 under-five child-years.
- [Inclusion flow](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/study_flow/study_flow_diagram.png) and [caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/study_flow/CAPTION.md).
- [PfPR and mortality by age band](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/pfpr_splines.png) (formerly main Figure 2), captioned in the [paper-figure captions](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/paper_figures/CAPTIONS.md).
- [PfPR splines by subgroup](../results/cbh/subgroups_dhsmics_map_gamma2_v3/sfig_pfpr_splines_by_subgroup.png) and [caption](../results/cbh/subgroups_dhsmics_map_gamma2_v3/CAPTION.md): the primary curves overlaid with refits in the Sahel, Eastern Africa and surveys before/after the median survey year (see below).
- [PfPR splines without wasting and stunting](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/sfig_pfpr_splines_without_nutrition.png) and [caption](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/CAPTION.md) (see below).
- Complete-case versus imputed-covariate PfPR curves (v8): written by `R_cbh/sensitivity/imputation_dhsmics/06_report.R` as `sfig_pfpr_splines_imputed_covariates.png` and `CAPTION_sfig_imputed_covariates.md` in `results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/`.

Main Figure 4 (added 24 September 2026), the probability of dying before age 5 in 2024 by country, is a primary result: `results/cbh/primary_map_regional17_dhsmics_gamma2_v7/under5_probability/` ([caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/under5_probability/CAPTION.md)).

## Subgroup refits of the primary model at gamma=2

Added 23 September 2026 on the DHS and MICS sample; refitted on v7 on 24 September 2026 as `subgroups_dhsmics_map_gamma2_v3`. The seven 17-variable primary age-band models refitted separately in four subsets of the v7 sample. Same covariates and full-sample scaling, fixed HIV imputation, gamma=2 and basis dimensions; knots at the subset's own predictor quantiles. Whole regions and whole surveys enter. MICS region centroids come from the analysis-region polygons (`R_mics/11_region_centroids.R`). The subsets overlap the full sample and each other, so differences are descriptive.

Subset sizes: Sahel (centroid ≥12°N, west of 36°E, outside the Horn of Africa) 31 surveys, 8 countries, 18,919 deaths; UN M49 Eastern Africa 51 surveys, 12 countries, 31,589 deaths (Liberia is in neither subset, so both are unchanged from v2); surveys up to the median survey year, 2014: 72 surveys; after 2014: 63 surveys.

- [Overlay figure](../results/cbh/subgroups_dhsmics_map_gamma2_v3/sfig_pfpr_splines_by_subgroup.png) and [40%→20% contrasts by subgroup](../results/cbh/subgroups_dhsmics_map_gamma2_v3/pfpr_40_to_20_by_subgroup.png)
- [Report with contrast table and subset sizes](../results/cbh/subgroups_dhsmics_map_gamma2_v3/REPORT.md) and [sample summary](../results/cbh/subgroups_dhsmics_map_gamma2_v3/sample_summary.csv)
- [Region selection](../results/cbh/subgroups_dhsmics_map_gamma2_v3/region_selection.csv) and [survey groups](../results/cbh/subgroups_dhsmics_map_gamma2_v3/survey_groups.csv)
- Reproduction: `Rscript R_cbh/sensitivity/subgroups/01_fit.R` then `02_report.R` ([instructions](../R_cbh/sensitivity/README.md))

## Adjustment without wasting and stunting

Added 23 September 2026 on the DHS and MICS sample; refitted on v7 on 24 September 2026 as `nutrition_adjustment_dhsmics_map_gamma2_v3`. The seven primary age-band models refitted on the identical v7 complete-case sample with the two nutritional covariates removed (15 covariates instead of 17). Wasting and stunting are measured in surviving children at the time of the survey and may lie on the causal pathway from malaria to death, so adjusting for them risks removing part of the association of interest. Only the adjustment set changes: sample, exposure, reference knots, basis dimensions, gamma=2 and random effects are those of the primary. Results: see the [hazard-ratio comparison](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/CONTRASTS.md). See plan section 3.4.4.

- [Supplementary figure](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/sfig_pfpr_splines_without_nutrition.png) and [caption](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/CAPTION.md)
- [Hazard-ratio comparison](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/CONTRASTS.md) and [contrast table](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/comparison_contrasts.csv)
- [Fit diagnostics](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/fit_diagnostics.csv) and [restart record](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/restarts.csv)
- Reproduction: `Rscript R_cbh/sensitivity/nutrition/01_fit.R` then `02_report.R`

## Imputed-covariate sensitivity version

Refitted on 24 September 2026 as `primary_map_regional17_dhsmics_imputed_gamma2_v8`, compared with v7; it replaces v6 (below), which was compared with v5. The 17-variable primary specification refitted after imputing every remaining covariate gap, so that every MAP-eligible record is retained: the same 8,797,963 records, 2,656,472 children, 123,419 deaths, 1,457 survey-regions, 166 surveys and 40 countries as v6. The only input change from v6 is Liberia's child HIV incidence, which now comes from the UNAIDS counts and is held fixed across the 10 imputations; São Tomé and Príncipe keeps the latent-adolescent imputation. The complete-case fit (v7) remains the primary. The covariate imputations are those of v6. Results: see `REPORT.md` and `multiple_imputation/REPORT.md` in `results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/`, written by `06_report.R` and `05_propagate.R`. See plan sections 2.4 and 3.4.3 and the [code instructions](../R_cbh/sensitivity/imputation_dhsmics/README.md).

- [Sample](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/prepared_sample.csv), [records with imputed values](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/imputation_record_counts.csv) and [by survey](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/imputation_record_counts_by_survey.csv)
- [Fit diagnostics](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/fit_diagnostics.csv)
- [Imputation audits: regional](../results/cbh/covariate_imputation_dhsmics_v6/REGIONAL_IMPUTATION.md) and [national](../results/cbh/covariate_imputation_dhsmics_v6/NATIONAL_IMPUTATION.md), with the [Somalia, South Sudan and Zimbabwe series](../results/cbh/covariate_imputation_dhsmics_v6/somalia_south_sudan_zimbabwe_series.csv)
- [Liberia's UNAIDS-derived HIV incidence](../results/cbh/hiv_incidence/liberia_aidsinfo/REPORT.md); [extended HIV incidence panel](../results/cbh/hiv_incidence/imputed_no_adolescent_series.csv) and [latent-adolescent validation](../results/cbh/hiv_incidence/latent_adolescent_validation.csv)

## SMC introduction before/after (DHS and MICS)

Added 23 September 2026; a supplementary result that is not part of the paper (24 September 2026). Refitted on v7 on 24 September 2026 as `smc_dhsmics_map_gamma2_v2` (v1, on v5, is below). The seven primary age-band models refitted in the eight countries with SMC in the DataWell cluster file, with an admin-1 indicator for band entry at or after the unit's SMC switch-on (Nigeria by state). Main analysis: confirmed (district- or region-level) switch-ons only; national-scope years, a coverage share and a 3-year-early placebo are sensitivities. Plan Section 3.4.5.

- [Report](../results/cbh/smc_dhsmics_map_gamma2_v2/REPORT.md) and [hazard ratios by age](../results/cbh/smc_dhsmics_map_gamma2_v2/smc_hazard_ratios_by_age.png)
- [Serial coverage maps, 2015–2022](../results/cbh/smc_dhsmics_map_gamma2_v2/smc_coverage_maps_2015_2022.png) and [caption](../results/cbh/smc_dhsmics_map_gamma2_v2/CAPTION_smc_coverage_maps.md)
- [Admin-1 switch-on years](../results/cbh/smc_dhsmics_map_gamma2_v2/admin1_smc_years.csv) and [sample by country](../results/cbh/smc_dhsmics_map_gamma2_v2/country_summary.csv)

## Sahelian child mortality by calendar month

Added 16 September 2026. A descriptive analysis of the all-cause death rate of children under two years by month of the year in Sahelian survey regions (boundary centroid at or above 11°N, west of the Horn of Africa; 247 regions in 41 surveys), by age band (<1, 1–5, 6–11, 12–23 completed months). Month of death is derived from month of birth and reported age at death; see [R_cbh/seasonality/README.md](../R_cbh/seasonality/README.md) for the derivation and its limits. It fits no prevalence model and uses DHS and MIS surveys only (no MICS).

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

## DHS and MICS v5 history

The four sensitivities below were fitted on the 23 September DHS+MICS primary `primary_map_regional17_dhsmics_gamma2_v5` (132 surveys: 95 DHS and 37 MICS; 36 countries; 7,498,459 child-band records, 102,282 deaths) before Liberia entered the primary. They are superseded by the v7 refits above and kept unchanged. The v5 primary outputs, including its [count version of Figure 3](../results/cbh/primary_map_regional17_dhsmics_gamma2_v5/burden_comparison/sfig_burden_comparison_counts.png), [inclusion flow](../results/cbh/primary_map_regional17_dhsmics_gamma2_v5/study_flow/study_flow_diagram.png) and [paper-figure captions](../results/cbh/primary_map_regional17_dhsmics_gamma2_v5/paper_figures/CAPTIONS.md), are preserved.

### Subgroup refits on v5 (v2)

Added 23 September 2026.

| Subset | Children | Records | Deaths | Survey-regions | Surveys | Countries |
|---|---:|---:|---:|---:|---:|---:|
| Sahel (centroid ≥12°N, west of 36°E, outside the Horn of Africa) | 468,275 | 1,538,691 | 18,919 | 239 | 31 | 8 |
| UN M49 Eastern Africa | 828,170 | 2,680,402 | 31,589 | 447 | 51 | 12 |
| Surveys up to the median survey year, 2014 | 1,076,607 | 3,399,915 | 55,805 | 549 | 70 | 33 |
| Surveys after 2014 | 1,215,482 | 4,098,544 | 46,477 | 662 | 62 | 34 |

The Sahel subset covers Burkina Faso, the Gambia, Guinea-Bissau, Mali, Mauritania, Niger, Senegal and Chad. All 28 subgroup fits converged; three (Sahel and Eastern Africa at 48–59 months, surveys after 2014 at 12–23 months) used the tighter-tolerance restart.

- [Overlay figure](../results/cbh/subgroups_dhsmics_map_gamma2_v2/sfig_pfpr_splines_by_subgroup.png) and [40%→20% contrasts by subgroup](../results/cbh/subgroups_dhsmics_map_gamma2_v2/pfpr_40_to_20_by_subgroup.png)
- [Report with contrast table and subset sizes](../results/cbh/subgroups_dhsmics_map_gamma2_v2/REPORT.md) and [sample summary](../results/cbh/subgroups_dhsmics_map_gamma2_v2/sample_summary.csv)
- [Region selection](../results/cbh/subgroups_dhsmics_map_gamma2_v2/region_selection.csv) and [survey groups](../results/cbh/subgroups_dhsmics_map_gamma2_v2/survey_groups.csv)

### Adjustment without wasting and stunting on v5 (v2)

Added 23 September 2026. The seven primary age-band models refitted on the identical v5 complete-case sample (7,498,459 records, 102,282 deaths, 1,211 survey-regions, 132 surveys, 36 countries) without wasting and stunting. The largest absolute change in the 40%→20% hazard ratio is 0.008, at 1–5 months (0.915 with all 17 covariates, 0.907 without). The 48–59-month fit used the tighter-tolerance restart.

- [Supplementary figure](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/sfig_pfpr_splines_without_nutrition.png) and [caption](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/CAPTION.md)
- [Hazard-ratio comparison](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/CONTRASTS.md) and [contrast table](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/comparison_contrasts.csv)
- [Fit diagnostics](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/fit_diagnostics.csv) and [restart record](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2/restarts.csv)

### Imputed-covariate version on v5 (v6)

Added 23 September 2026: `primary_map_regional17_dhsmics_imputed_gamma2_v6`, compared with v5, on the same 8,797,963 records as v8.

- Regional gaps: chained-equation multiple imputation (`mice`, predictive mean matching, m = 10) on the combined DHS and MICS survey-region overlay, imputed jointly. Imputed survey-regions: wasting 111, stunting 89, wealth 66, facility delivery 36, electricity 12. The MICS gaps are the whole-survey anthropometry gaps of Nigeria 2021, Madagascar 2012 (South), and Somalia North-East and Somaliland 2011. 795,848 records (558,906 DHS, 236,942 MICS) carry at least one imputed regional value.
- National gaps: the 2001 WGI round by interpolation; South Sudan's political stability before independence (2000–2010) and its GDP outside 2008–2015 from GAMs; health expenditure for Zimbabwe 2000–2009, Somalia 2000–2012, South Sudan outside 2017–2023 and every country in 2024 from a GAM. Records affected: political stability 247,862; GDP 33,840 (MICS South Sudan 2010 only); health expenditure 260,385.
- Child HIV incidence for countries without an adolescent series (Liberia; São Tomé and Príncipe, including MICS 2014 and 2019) from the extended incidence model: 165,149 records. (v8 takes Liberia's rate from UNAIDS counts instead.)
- The South Sudan pre-independence values are a modelling choice. Alternatives are whole-Sudan values, carrying the 2011 values back, or omitting MICS South Sudan 2010 (57,509 records, 0.65% of the sample).
- All seven point fits converged without restarts. Point-fit hazard ratios for PfPR 40%→20% by age band (<1 to 48–59 months): 0.940, 0.913, 0.825, 0.843, 0.864, 0.858, 0.947, against 0.944, 0.915, 0.830, 0.823, 0.850, 0.843, 0.907 in v5.
- Multiple-imputation check (pooled curves and contrasts): [report](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/multiple_imputation/REPORT.md) and [pooled contrasts](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/multiple_imputation/pooled_contrasts.csv); 10 imputed datasets changed the hazard ratios by at most 0.002, with a between-imputation variance share of at most 2.3%
- Comparison report and supplementary figure: [report](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/REPORT.md), [figure](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/sfig_pfpr_splines_imputed_covariates.png) and [caption](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/CAPTION_sfig_imputed_covariates.md); largest change in the 40%→20% hazard ratio from v5 0.040 (48–59 months)
- [Sample](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/prepared_sample.csv), [records with imputed values](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/imputation_record_counts.csv) and [by survey](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/imputation_record_counts_by_survey.csv); [fit diagnostics](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/fit_diagnostics.csv) and [point-fit curves](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6/pfpr_curves.csv)

### SMC before/after on v5 (v1)

Added 23 September 2026. No age band showed a statistically detectable association (hazard ratios 0.83–1.15, all intervals include 1); placebo hazard ratios 0.92–1.07.

- [Report](../results/cbh/smc_dhsmics_map_gamma2_v1/REPORT.md) and [hazard ratios by age](../results/cbh/smc_dhsmics_map_gamma2_v1/smc_hazard_ratios_by_age.png)
- [Serial coverage maps, 2015–2022](../results/cbh/smc_dhsmics_map_gamma2_v1/smc_coverage_maps_2015_2022.png) and [caption](../results/cbh/smc_dhsmics_map_gamma2_v1/CAPTION_smc_coverage_maps.md)
- [Admin-1 switch-on years](../results/cbh/smc_dhsmics_map_gamma2_v1/admin1_smc_years.csv) and [sample by country](../results/cbh/smc_dhsmics_map_gamma2_v1/country_summary.csv)

## DHS-only history

The three sensitivities below were fitted on the DHS-only primary `primary_map_regional17_gamma2_v3` (5,465,305 records, 1,686,004 children, 75,726 deaths, 916 survey-regions, 95 surveys, 34 countries) before MICS surveys entered the primary. They are superseded by the DHS and MICS refits above and kept unchanged. `--dhs-only` on the subgroup and nutrition scripts reproduces them; the imputed version runs through `R_cbh/primary/run_regional.R --imputed`.

### Subgroup refits on the DHS-only sample (v1)

Added 18 September 2026. The seven 17-variable primary age-band models refitted separately in four subsets of the DHS-only fitted sample: Sahel survey regions (centroid ≥12°N, west of 36°E, outside the Horn of Africa; 23 surveys in seven countries), UN M49 Eastern Africa (38 surveys in 12 countries), surveys up to the median survey year of 2013 (50 surveys) and surveys after it (45 surveys). Same covariates and scaling, fixed HIV imputation, gamma=2 and basis dimensions; knots at the subset's own predictor quantiles. Whole regions and whole surveys enter. The subsets overlap the full sample, so differences are descriptive.

- [Overlay figure](../results/cbh/subgroups_map_gamma2_v1/sfig_pfpr_splines_by_subgroup.png) and [40%→20% contrasts by subgroup](../results/cbh/subgroups_map_gamma2_v1/pfpr_40_to_20_by_subgroup.png)
- [Report with contrast table and subset sizes](../results/cbh/subgroups_map_gamma2_v1/REPORT.md)
- [Region and survey selection audits](../results/cbh/subgroups_map_gamma2_v1/region_selection.csv)

### Imputed-covariate version on the DHS-only data (v4)

Added 18 September 2026. The 17-variable primary specification refitted after imputing every remaining covariate gap (whole-survey gaps in wasting, stunting, facility delivery, electricity and wealth by chained-equation multiple imputation; the 2001 WGI round by interpolation; health expenditure for Zimbabwe 2000–2009 and all countries' 2024 from a GAM; child HIV incidence for Liberia and São Tomé from the extended incidence model), so that all 6,357,802 MAP-eligible DHS and MIS records from 120 surveys in 36 countries are retained. This is `primary_map_regional17_imputed_gamma2_v4`, compared with v3.

- [Comparison with the complete-case primary](../results/cbh/primary_map_regional17_imputed_gamma2_v4/REPORT.md)
- [Multiple-imputation check: pooled curves and contrasts](../results/cbh/primary_map_regional17_imputed_gamma2_v4/multiple_imputation/REPORT.md)
- [Supplementary figure: complete-case versus imputed-covariate PfPR curves](../results/cbh/primary_map_regional17_imputed_gamma2_v4/sfig_pfpr_splines_imputed_covariates.png) and [caption](../results/cbh/primary_map_regional17_imputed_gamma2_v4/CAPTION_sfig_imputed_covariates.md)
- [Imputation audits: regional](../results/cbh/covariate_imputation_v4/REGIONAL_IMPUTATION.md) and [national](../results/cbh/covariate_imputation_v4/NATIONAL_IMPUTATION.md)

### Adjustment without wasting and stunting on the DHS-only sample (v1)

Added 22 September 2026. The seven primary age-band models refitted on the identical DHS-only complete-case sample with wasting and stunting removed (15 covariates instead of 17); only the adjustment set changes.

- [Supplementary figure](../results/cbh/nutrition_adjustment_map_gamma2_v1/sfig_pfpr_splines_without_nutrition.png) and [caption](../results/cbh/nutrition_adjustment_map_gamma2_v1/CAPTION.md)
- [Hazard-ratio comparison](../results/cbh/nutrition_adjustment_map_gamma2_v1/CONTRASTS.md) and [contrast table](../results/cbh/nutrition_adjustment_map_gamma2_v1/comparison_contrasts.csv)
- [Fit diagnostics](../results/cbh/nutrition_adjustment_map_gamma2_v1/fit_diagnostics.csv)

## Removed analyses

The comparison of MAP with the annual Snow prevalence estimates (matched MAP–Snow fits, Snow penalisation trials and Snow-based national mortality scenarios) was removed from the analysis plan on 18 September 2026 and is not part of the manuscript. Its code is in [`archive/2026-09-18-snow-comparison/`](../archive/2026-09-18-snow-comparison/README.md). The saved outputs in `results/cbh/map_snow_gamma2_v1/` and `results/cbh/age_band_snow_2000_2015_v1/` are retained unchanged as historical records; `results/cbh/map_snow_gamma2_v1/` also holds the original primary `map_full` fits and is still read by the saved-primary audit.
