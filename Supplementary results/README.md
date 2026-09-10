# Supplementary results

Reporting designation updated 10 September 2026. The primary analysis uses **full-sample MAP prevalence, seven separate age-band models and gamma=2**; see the [analysis plan](../docs/ANALYSIS_PLAN.md). The Snow prevalence analyses below are supplementary. These links retain the original files and provenance; no fits or results have been moved or recomputed.

## Snow versus MAP at gamma=2

The principal supplementary comparison uses identical pre-2016 records. Its matched MAP fit is a supplementary comparator, not the full-sample primary MAP model.

- [Matched age-band spline figure](../results/cbh/map_snow_gamma2_v1/matched_pfpr_splines.png)
- [Full-sample comparison figure](../results/cbh/map_snow_gamma2_v1/full_sample_pfpr_splines.png) — changes fitting period/sample as well as exposure source
- [Model comparison report, EDF, contrasts and diagnostics](../results/cbh/map_snow_gamma2_v1/REPORT.md)
- [Country mortality comparison report](../results/cbh/map_snow_gamma2_v1/burden/REPORT.md)
- [All-country mortality figure](../results/cbh/map_snow_gamma2_v1/burden/country_deaths_all_years.png)
- [Country comparison table](../results/cbh/map_snow_gamma2_v1/burden/country_comparisons.csv)

The country comparison evaluates both sets of curves at the **same national MAP exposure** in 2005, 2015 and 2024. Snow has no 2024 prevalence estimate; its 2024 mortality result transports the pre-2016 fitted relationship to MAP exposure. Preserve that label and the source/transport limitations. Primary MAP rows appearing as references in these supplementary comparisons remain identifiable as `map_full` (or `map_full_gamma2` in year totals).

## Snow exposure audit and smoothing analyses

- [Annual Snow extraction, joins and original comparison](../results/cbh/age_band_snow_2000_2015_v1/REPORT.md)
- [Separate gamma=1.4 and PfPR cs sensitivities](../results/cbh/age_band_snow_2000_2015_v1/penalty_sensitivity/REPORT.md)
- [Snow gamma=1, 1.4 and 2 comparison](../results/cbh/age_band_snow_2000_2015_v1/penalty_sensitivity/gamma2/comparison/REPORT.md)

The Snow cs fit used gamma=1; it changed only the PfPR basis. All Snow analyses retain their actual settings. MAP/gamma=2/cr is the primary specification.
