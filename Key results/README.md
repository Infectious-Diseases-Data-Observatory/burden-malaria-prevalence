# Key results — primary MAP analysis

The primary analysis uses seven separate age-band models, MAP prevalence, cr PfPR splines and gamma=2. The full sample contains 1,817,912 children, 5,885,022 child-band records and 82,415 deaths from 105 surveys in 34 countries.

[Primary analysis report and rerun verification](../results/cbh/primary_map_gamma2_v1/REPORT.md)

| Result | Output |
|---|---|
| Figure 1: survey map and timing | [Figure](../results/cbh/primary_map_gamma2_v1/survey_map/survey_map_and_timing.png) and [caption](../results/cbh/primary_map_gamma2_v1/survey_map/CAPTION.md) |
| Inclusion flow | [Figure](../results/cbh/primary_map_gamma2_v1/study_flow/study_flow_diagram.png) and [caption](../results/cbh/primary_map_gamma2_v1/study_flow/CAPTION.md) |
| Figure 2: seven PfPR curves | [Figure](../results/cbh/primary_map_gamma2_v1/pfpr_splines.png) and [estimates](../results/cbh/primary_map_gamma2_v1/pfpr_curves.csv) |
| Figure 3: malaria-attributable share by age | [Figure](../results/cbh/primary_map_gamma2_v1/attributable_fraction_by_age.png) and [estimates](../results/cbh/primary_map_gamma2_v1/attributable_fraction_by_age.csv) |
| Age-band results table | [Table](../results/cbh/primary_map_gamma2_v1/tables/age_band_results.md) and [LaTeX fragment](../results/cbh/primary_map_gamma2_v1/tables/age_band_results.tex) |
| PfPR 40% to 20% effects | [Figure](../results/cbh/primary_map_gamma2_v1/pfpr_40_to_20.png) and [contrasts](../results/cbh/primary_map_gamma2_v1/pfpr_40_to_20_contrasts.csv) |
| Diagnostics | [Numerical checks](../results/cbh/primary_map_gamma2_v1/fit_diagnostics.csv) and [fitted outcomes](../results/cbh/primary_map_gamma2_v1/fitted_outcome_checks.csv) |
| Mortality, 2005/2015/2024 | [Country figure](../results/cbh/primary_map_gamma2_v1/burden/country_deaths_all_years.png) and [country totals](../results/cbh/primary_map_gamma2_v1/burden/country_totals.csv) |
| Figure 4: comparison with IHME | [Scatter plots](../results/cbh/primary_map_gamma2_v1/burden/country_vs_ihme.png) |
| Figure 5: annual mortality, 2004–2024 | [Figure](../results/cbh/primary_map_gamma2_v1/annual_comparison/fig5_annual_malaria_mortality.png) and [totals and input audit](../results/cbh/primary_map_gamma2_v1/annual_comparison/README.md) |
| Age contributions | [Figure](../results/cbh/primary_map_gamma2_v1/burden/deaths_by_age.png) and [country-age estimates](../results/cbh/primary_map_gamma2_v1/burden/country_age_estimates.csv) |
| DRC synthetic-cohort survival | [Figure](../results/cbh/primary_map_gamma2_v1/burden/drc_survival.png) |

[Paper figure manifest and captions](../results/cbh/primary_map_gamma2_v1/paper_figures/CAPTIONS.md)

Snow exposure comparisons are [supplementary](<../Supplementary results/README.md>). The previous MAP gamma=2 outputs remain in `results/cbh/map_snow_gamma2_v1/`, with a numerical comparison in the new report.

## Historical figures

Older PNGs in this directory are retained for traceability. The gamma=1 [geography/period/structure figure](pfpr_spline_sensitivity_by_age_geography_period.png) remains a key exploratory result; planned gamma=2 replacements have not been run. See the [original index](../archive/2026-09-15-code-audit/Key%20results/README.md).

Run `Rscript run_all.R --primary` to refit and regenerate the primary results from prepared inputs. Use `Rscript R_cbh/primary/run.R --report-only` for saved fits. The standalone reporting scripts and old `R_dhs/11_study_flow.R` and `R_dhs/33_key_results.R` wrappers use the same primary selection.
