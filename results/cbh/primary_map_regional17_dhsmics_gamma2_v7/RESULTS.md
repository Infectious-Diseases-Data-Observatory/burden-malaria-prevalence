# Current primary figures and tables

The revised 17-variable PfPR-ACM model uses 2,325,130 children, 7,607,122 child-band records and 103,987 deaths in 135 surveys and 37 countries.

Seven separate MAP gamma=2 age-band fits, including survey-region urban percentage, with fixed PfPR/time knots and the declared HIV, UNICEF and available-region substitutions. The raw-source data and HIV imputation were not refitted. All seven models passed the fitted-input and numerical checks.

![PfPR curves](pfpr_splines.png)

| Completed months | Records | Deaths (% of U5 deaths) | PfPR EDF | HR: 40% to 20% (95% interval) | HR: 20% to 0% (95% interval) |
|---|---:|---:|---:|---:|---:|
| <1 | 1,341,316 | 38,245 (36.8%) | 2.33 | 0.94 (0.92–0.97) | 0.94 (0.90–0.99) |
| 1-5 | 1,206,366 | 15,873 (15.3%) | 2.28 | 0.92 (0.88–0.95) | 0.86 (0.80–0.93) |
| 6-11 | 1,160,005 | 14,384 (13.8%) | 2.97 | 0.83 (0.79–0.88) | 0.69 (0.63–0.76) |
| 12-23 | 1,003,279 | 13,389 (12.9%) | 3.51 | 0.83 (0.78–0.88) | 0.53 (0.47–0.60) |
| 24-35 | 994,653 | 11,488 (11.0%) | 3.66 | 0.86 (0.81–0.92) | 0.47 (0.41–0.54) |
| 36-47 | 963,051 | 6,732 (6.5%) | 3.33 | 0.85 (0.79–0.91) | 0.55 (0.48–0.64) |
| 48-59 | 938,452 | 3,876 (3.7%) | 3.20 | 0.92 (0.84–1.01) | 0.57 (0.49–0.67) |
| **Total** | **7,607,122** | **103,987 (100.0%)** | — | — | — |

Deaths are observed deaths in the primary analysis sample; percentages use all 103,987 under-five deaths as the denominator. Displayed percentages use largest-remainder rounding to one decimal place and sum to 100.0%. Records are child–age-band observations, so a child can contribute multiple records. EDF: effective degrees of freedom of the PfPR spline. Both contrast columns are adjusted mortality hazard ratios for reducing PfPR from the first value to the second. Intervals are conditional on fitted smoothing parameters, exposure, the fixed HIV imputation and filled regional covariates; zero PfPR is below observed exposure support in every band.

[Table CSV](tables/age_band_results.csv) · [LaTeX source as plain text](tables/age_band_results.latex.txt)

[40% to 20% contrast figure](pfpr_40_to_20.png) · [Full contrasts and comparison with the previous adjustment](REPORT.md)

![Country estimates](burden/country_vs_ihme.png)

Country/IHME axes use identical base-10 scales starting at 1,000 deaths. All country estimates, signed age contributions and missing/support flags remain in the CSVs. Age intervals are conditional on smoothing parameters, exposure and filled covariates; marginal intervals are not summed into total intervals.

[Country totals](burden/country_totals.csv) · [Country-age estimates](burden/country_age_estimates.csv) · [Annual comparison](annual_comparison/README.md) · [Nigeria states](nigeria_states/README.md)

[Age contributions](burden/deaths_by_age.png) · [DRC synthetic-cohort survival](burden/drc_survival.png) · [Aggregate outcome check](outcome_residuals_by_pfpr.png)

[Paper figure manifest/captions](paper_figures/CAPTIONS.md) · [Study flow](study_flow/CAPTION.md) · [Fit diagnostics](fit_diagnostics.csv)

These are primary results. Sahel, geographic, period and other sensitivity analyses are separate; they are not refitted by reporting. No TeX files are written. Reproduce with `Rscript R_cbh/primary/run_regional.R --report-only`.
