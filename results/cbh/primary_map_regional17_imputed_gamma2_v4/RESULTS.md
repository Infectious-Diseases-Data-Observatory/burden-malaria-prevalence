# Current primary figures and tables

The imputed-covariate (sensitivity) 17-variable PfPR-ACM model uses 1,932,325 children, 6,357,802 child-band records and 90,938 deaths in 120 surveys and 36 countries.

Seven separate MAP gamma=2 age-band fits, including survey-region urban percentage, with fixed PfPR/time knots and the declared HIV, UNICEF and available-region substitutions. The raw-source data and HIV imputation were not refitted. All seven models passed the fitted-input and numerical checks.

![PfPR curves](pfpr_splines.png)

| Completed months | Records | Deaths (% of U5 deaths) | PfPR EDF | HR: 40% to 20% (95% interval) | HR: 20% to 0% (95% interval) |
|---|---:|---:|---:|---:|---:|
| <1 | 1,131,361 | 32,994 (36.3%) | 2.28 | 0.95 (0.92–0.98) | 0.94 (0.89–0.99) |
| 1-5 | 1,011,107 | 13,525 (14.9%) | 2.04 | 0.92 (0.88–0.96) | 0.90 (0.85–0.97) |
| 6-11 | 967,350 | 12,573 (13.8%) | 2.42 | 0.81 (0.78–0.85) | 0.76 (0.70–0.82) |
| 12-23 | 835,276 | 12,383 (13.6%) | 3.39 | 0.85 (0.80–0.91) | 0.59 (0.52–0.67) |
| 24-35 | 828,103 | 9,999 (11.0%) | 3.53 | 0.85 (0.79–0.92) | 0.51 (0.44–0.58) |
| 36-47 | 802,270 | 5,981 (6.6%) | 3.22 | 0.85 (0.78–0.92) | 0.58 (0.50–0.68) |
| 48-59 | 782,335 | 3,483 (3.8%) | 3.30 | 0.98 (0.89–1.08) | 0.57 (0.48–0.68) |
| **Total** | **6,357,802** | **90,938 (100.0%)** | — | — | — |

Deaths are observed deaths in the primary analysis sample; percentages use all 90,938 under-five deaths as the denominator. Displayed percentages use largest-remainder rounding to one decimal place and sum to 100.0%. Records are child–age-band observations, so a child can contribute multiple records. EDF: effective degrees of freedom of the PfPR spline. Both contrast columns are adjusted mortality hazard ratios for reducing PfPR from the first value to the second. Intervals are conditional on fitted smoothing parameters, exposure, the fixed HIV imputation and filled regional covariates; zero PfPR is below observed exposure support in every band.

[Table CSV](tables/age_band_results.csv) · [LaTeX source as plain text](tables/age_band_results.latex.txt)

[40% to 20% contrast figure](pfpr_40_to_20.png) · [Full contrasts and comparison with the previous adjustment](REPORT.md)

![Country estimates](burden/country_vs_ihme.png)

Country/IHME axes use identical base-10 scales starting at 1,000 deaths. All country estimates, signed age contributions and missing/support flags remain in the CSVs. Age intervals are conditional on smoothing parameters, exposure and filled covariates; marginal intervals are not summed into total intervals.

[Country totals](burden/country_totals.csv) · [Country-age estimates](burden/country_age_estimates.csv) · [Annual comparison](annual_comparison/README.md) · [Nigeria states](nigeria_states/README.md)

[Age contributions](burden/deaths_by_age.png) · [DRC synthetic-cohort survival](burden/drc_survival.png) · [Aggregate outcome check](outcome_residuals_by_pfpr.png)

[Paper figure manifest/captions](paper_figures/CAPTIONS.md) · [Study flow](study_flow/CAPTION.md) · [Fit diagnostics](fit_diagnostics.csv)

These are primary results. Sahel, geographic, period and other sensitivity analyses are separate; they are not refitted by reporting. No TeX files are written. Reproduce with `Rscript R_cbh/primary/run_regional.R --report-only`.
