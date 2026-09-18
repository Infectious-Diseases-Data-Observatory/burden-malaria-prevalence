# Key results — primary MAP analysis

The primary analysis uses seven separate age-band MAP models, cr PfPR splines, gamma=2 and 17 regional/annual covariates. The sample contains 1,686,004 children, 5,465,305 child-band records and 75,726 deaths from 95 surveys in 34 countries.

[Current primary figures and tables](../results/cbh/primary_map_regional17_gamma2_v3/RESULTS.md)

[Comparison with the previous iteration](../results/cbh/primary_map_regional17_gamma2_v3/REPORT.md)

| Result | Output |
|---|---|
| Figure 1: survey map and timing | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/survey_map/survey_map_and_timing.png) and [caption](../results/cbh/primary_map_regional17_gamma2_v3/survey_map/CAPTION.md) |
| Inclusion flow | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/study_flow/study_flow_diagram.png) and [caption](../results/cbh/primary_map_regional17_gamma2_v3/study_flow/CAPTION.md) |
| Figure 2: seven PfPR curves | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/pfpr_splines.png) and [estimates](../results/cbh/primary_map_regional17_gamma2_v3/pfpr_curves.csv) |
| Figure 3: malaria-attributable share by age | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/attributable_fraction_by_age.png) and [estimates](../results/cbh/primary_map_regional17_gamma2_v3/attributable_fraction_by_age.csv) |
| Age-band results table | [Table](../results/cbh/primary_map_regional17_gamma2_v3/tables/age_band_results.md) and [LaTeX source as text](../results/cbh/primary_map_regional17_gamma2_v3/tables/age_band_results.latex.txt) |
| Country, annual and state tables | [Table index](../results/cbh/primary_map_regional17_gamma2_v3/tables/README.md) |
| PfPR 40% to 20% effects | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/pfpr_40_to_20.png) and [contrasts](../results/cbh/primary_map_regional17_gamma2_v3/pfpr_40_to_20_contrasts.csv) |
| Diagnostics | [Numerical checks](../results/cbh/primary_map_regional17_gamma2_v3/fit_diagnostics.csv) and [fitted outcomes](../results/cbh/primary_map_regional17_gamma2_v3/fitted_outcome_checks.csv) |
| Mortality, 2005/2015/2024 | [Country figure](../results/cbh/primary_map_regional17_gamma2_v3/burden/country_deaths_all_years.png) and [country totals](../results/cbh/primary_map_regional17_gamma2_v3/burden/country_totals.csv) |
| Figure 4: countries (A), Nigerian states (B) and annual trend (C) versus IHME/UN IGME | [Combined figure](../results/cbh/primary_map_regional17_gamma2_v3/burden_comparison/fig4_burden_comparison.png) and [caption](../results/cbh/primary_map_regional17_gamma2_v3/burden_comparison/CAPTION.md) |
| Figure 4A source: country comparison, 2005/2015/2024 | [Three-year scatter](../results/cbh/primary_map_regional17_gamma2_v3/burden/country_vs_ihme.png) |
| Figure 4B source: Nigerian states, 2024 | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/nigeria_states/fig6_nigeria_states_vs_ihme.png) and [state table and audit](../results/cbh/primary_map_regional17_gamma2_v3/nigeria_states/README.md) |
| Figure 4C source: annual mortality, 2004–2024 | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/annual_comparison/fig5_annual_malaria_mortality.png) and [totals and input audit](../results/cbh/primary_map_regional17_gamma2_v3/annual_comparison/README.md) |
| Age contributions | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/burden/deaths_by_age.png) and [country-age estimates](../results/cbh/primary_map_regional17_gamma2_v3/burden/country_age_estimates.csv) |
| DRC synthetic-cohort survival | [Figure](../results/cbh/primary_map_regional17_gamma2_v3/burden/drc_survival.png) |

[Paper figure manifest and captions](../results/cbh/primary_map_regional17_gamma2_v3/paper_figures/CAPTIONS.md)

Sensitivity analyses are [supplementary](<../Supplementary results/README.md>). The previous MAP gamma=2 outputs remain in `results/cbh/primary_map_regional18_gamma2_v2/`, with a comparison in the new report.

## Historical figures

Older PNGs in this directory are retained for traceability. The gamma=1 [geography/period/structure figure](pfpr_spline_sensitivity_by_age_geography_period.png) remains a key exploratory result; planned gamma=2 replacements have not been run. See the [original index](../archive/2026-09-15-code-audit/Key%20results/README.md).

Run `Rscript run_all.R --primary` to refit and regenerate the primary results from prepared inputs. Primary figures and tables follow the current analysis plan; historical sensitivity plots are not relabelled. Use `Rscript R_cbh/primary/run_regional.R --report-only` for saved fits. The standalone reporting scripts and old `R_dhs/11_study_flow.R` and `R_dhs/33_key_results.R` wrappers use the same primary selection.
