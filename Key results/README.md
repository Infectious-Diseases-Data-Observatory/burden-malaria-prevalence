# Key results — primary MAP analysis

Seven separate age-band models; MAP PfPR at band entry; gamma=2; one median child HIV-incidence imputation. The primary sample contains 5,885,022 child-band records and 82,415 deaths from 105 surveys in 34 countries.

| Result | Authoritative output / selection |
|---|---|
| Inclusion flow | [Figure](../results/cbh/map_snow_gamma2_v1/study_flow/study_flow_diagram.png) and [caption](../results/cbh/map_snow_gamma2_v1/study_flow/CAPTION.md) |
| Seven PfPR curves | [Saved curve estimates](../results/cbh/map_snow_gamma2_v1/pfpr_curves.csv), `series == map_full` |
| PfPR 40% to 20% effects | [Contrasts](../results/cbh/map_snow_gamma2_v1/pfpr_40_to_20_contrasts.csv), `series == map_full` |
| Numerical diagnostics | [Diagnostics](../results/cbh/map_snow_gamma2_v1/fit_diagnostics.csv), `series == map_full` |
| Country burden, 2005/2015/2024 | [Country-age estimates](../results/cbh/map_snow_gamma2_v1/burden/country_age_estimates.csv) and [country totals](../results/cbh/map_snow_gamma2_v1/burden/country_totals.csv), `series == map_full`; [year totals](../results/cbh/map_snow_gamma2_v1/burden/year_summary.csv), `map_full_gamma2` |

A dedicated MAP-only curve and burden plotting stage is still needed. The combined MAP–Snow figures are supplementary. See the [supplementary index](<../Supplementary results/README.md>) and the [component-by-component audit](../docs/CODE_AUDIT.md).

## Historical figures

Existing older PNGs in this directory are retained for traceability, not designated primary results. Their [original index](../archive/2026-09-15-code-audit/Key%20results/README.md) is preserved. The gamma=1 [geography/period/structure figure](pfpr_spline_sensitivity_by_age_geography_period.png) remains a key exploratory result; its planned gamma=2 replacements have not been run.

Regenerate the inclusion figure with `Rscript R_cbh/reporting/01_study_flow.R` and this index with `Rscript R_cbh/reporting/02_results_index.R`. The older `R_dhs/11_study_flow.R` and `R_dhs/33_key_results.R` commands delegate to these scripts.
