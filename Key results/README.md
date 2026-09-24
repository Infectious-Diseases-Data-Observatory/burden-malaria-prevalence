# Key results — primary MAP analysis

The primary analysis uses seven separate age-band MAP models, cr PfPR splines, gamma=2 and 17 regional/annual covariates. The sample contains 2,325,130 children, 7,607,122 child-band records and 103,987 deaths from 135 surveys (98 DHS and 37 MICS) in 37 countries.

[Current primary figures and tables](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/RESULTS.md)

[Comparison with the previous iteration](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/REPORT.md)

| Result | Output |
|---|---|
| Figure 1: survey map and timing | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/survey_map/survey_map_and_timing.png) and [caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/survey_map/CAPTION.md) |
| Figure 2: malaria-attributable share by age | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/attributable_fraction_by_age.png) and [estimates](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/attributable_fraction_by_age.csv) |
| Figure 3: countries (A), Nigerian states (B) and annual trend (C) versus IHME/UN IGME, deaths per 1,000 child-years | [Combined figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden_comparison/fig4_burden_comparison.png) and [caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden_comparison/CAPTION.md) |
| Figure 3A source: country comparison, 2005/2015/2024 | [Three-year scatter](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden/country_vs_ihme.png) |
| Figure 3B source: Nigerian states, 2024 | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/nigeria_states/fig6_nigeria_states_vs_ihme.png) and [state table and audit](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/nigeria_states/README.md) |
| Figure 3C source: annual mortality, 2000–2024 | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/annual_comparison/fig5_annual_malaria_mortality.png) and [totals and input audit](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/annual_comparison/README.md) |
| Figure 4: probability of dying before age 5 by country, all causes and caused by malaria | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/under5_probability/under5_death_probability_2024.png) and [values](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/under5_probability/under5_death_probability.csv) and [caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/under5_probability/CAPTION.md) |
| Supplementary: inclusion flow | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/study_flow/study_flow_diagram.png) and [caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/study_flow/CAPTION.md) |
| Supplementary: seven PfPR curves | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/pfpr_splines.png) and [estimates](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/pfpr_curves.csv) |
| Supplementary: count version of Figure 3 | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden_comparison/sfig_burden_comparison_counts.png) and [caption](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden_comparison/CAPTION_counts.md) |
| Supplementary: subgroup refits | [Report](../results/cbh/subgroups_dhsmics_map_gamma2_v3/REPORT.md) |
| Supplementary: refit without nutrition covariates | [Contrasts](../results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v3/CONTRASTS.md) |
| Supplementary: imputed covariates | [Report](../results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v8/REPORT.md) |
| Age-band results table | [Table](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/tables/age_band_results.md) and [LaTeX source as text](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/tables/age_band_results.latex.txt) |
| Country, annual and state tables | [Table index](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/tables/README.md) |
| PfPR 40% to 20% effects | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/pfpr_40_to_20.png) and [contrasts](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/pfpr_40_to_20_contrasts.csv) |
| Diagnostics | [Numerical checks](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/fit_diagnostics.csv) and [fitted outcomes](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/fitted_outcome_checks.csv) |
| Mortality, 2005/2015/2024 | [Country figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden/country_deaths_all_years.png) and [country totals](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden/country_totals.csv) |
| Age contributions | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden/deaths_by_age.png) and [country-age estimates](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden/country_age_estimates.csv) |
| DRC synthetic-cohort survival | [Figure](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/burden/drc_survival.png) |

[Paper figure manifest and captions](../results/cbh/primary_map_regional17_dhsmics_gamma2_v7/paper_figures/CAPTIONS.md)

Sensitivity analyses are [supplementary](<../Supplementary results/README.md>). The comparator primary remains in `results/cbh/primary_map_regional17_dhsmics_gamma2_v5/`, with a comparison in the new report; this sample adds Liberia's 2007, 2013 and 2019-20 DHS to it (child HIV incidence derived from UNAIDS counts).

## Historical figures

Older PNGs in this directory are retained for traceability. The gamma=1 [geography/period/structure figure](pfpr_spline_sensitivity_by_age_geography_period.png) remains a key exploratory result; planned gamma=2 replacements have not been run. See the [original index](../archive/2026-09-15-code-audit/Key%20results/README.md).

Run `Rscript run_all.R --primary` to refit and regenerate the primary results from prepared inputs. Primary figures and tables follow the current analysis plan. Use `Rscript R_cbh/primary/run_regional.R --report-only` for saved fits. The standalone reporting scripts and old `R_dhs/11_study_flow.R` and `R_dhs/33_key_results.R` wrappers use the same primary selection. The sensitivity refits are run separately from `R_cbh/sensitivity/subgroups/`, `R_cbh/sensitivity/nutrition/` and `R_cbh/sensitivity/imputation_dhsmics/`.
