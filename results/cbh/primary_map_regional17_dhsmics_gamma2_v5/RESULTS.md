# Current primary figures and tables

The revised 17-variable PfPR-ACM model uses 2,292,089 children, 7,498,459 child-band records and 102,282 deaths in 132 surveys and 36 countries.

Seven separate MAP gamma=2 age-band fits, including survey-region urban percentage, with fixed PfPR/time knots and the declared HIV, UNICEF and available-region substitutions. The raw-source data and HIV imputation were not refitted. All seven models passed the fitted-input and numerical checks.

![PfPR curves](pfpr_splines.png)

| Completed months | Records | Deaths (% of U5 deaths) | PfPR EDF | HR: 40% to 20% (95% interval) | HR: 20% to 0% (95% interval) |
|---|---:|---:|---:|---:|---:|
| <1 | 1,322,149 | 37,623 (36.8%) | 2.34 | 0.94 (0.92–0.97) | 0.95 (0.90–1.00) |
| 1-5 | 1,189,209 | 15,586 (15.2%) | 2.25 | 0.92 (0.88–0.95) | 0.86 (0.81–0.93) |
| 6-11 | 1,143,915 | 14,126 (13.8%) | 2.90 | 0.83 (0.79–0.87) | 0.70 (0.63–0.77) |
| 12-23 | 989,211 | 13,138 (12.9%) | 3.46 | 0.82 (0.77–0.88) | 0.54 (0.48–0.62) |
| 24-35 | 980,344 | 11,338 (11.1%) | 3.63 | 0.85 (0.80–0.91) | 0.47 (0.41–0.54) |
| 36-47 | 949,032 | 6,650 (6.5%) | 3.34 | 0.84 (0.78–0.91) | 0.54 (0.47–0.63) |
| 48-59 | 924,599 | 3,821 (3.7%) | 3.16 | 0.91 (0.83–0.99) | 0.58 (0.49–0.68) |
| **Total** | **7,498,459** | **102,282 (100.0%)** | — | — | — |

Deaths are observed deaths in the primary analysis sample; percentages use all 102,282 under-five deaths as the denominator. Displayed percentages use largest-remainder rounding to one decimal place and sum to 100.0%. Records are child–age-band observations, so a child can contribute multiple records. EDF: effective degrees of freedom of the PfPR spline. Both contrast columns are adjusted mortality hazard ratios for reducing PfPR from the first value to the second. Intervals are conditional on fitted smoothing parameters, exposure, the fixed HIV imputation and filled regional covariates; zero PfPR is below observed exposure support in every band.

[Table CSV](tables/age_band_results.csv) · [LaTeX source as plain text](tables/age_band_results.latex.txt)

[40% to 20% contrast figure](pfpr_40_to_20.png) · [Full contrasts and comparison with the previous adjustment](REPORT.md)

![Country estimates](burden/country_vs_ihme.png)

Country/IHME axes use identical base-10 scales starting at 1,000 deaths. All country estimates, signed age contributions and missing/support flags remain in the CSVs. Age intervals are conditional on smoothing parameters, exposure and filled covariates; marginal intervals are not summed into total intervals.

[Country totals](burden/country_totals.csv) · [Country-age estimates](burden/country_age_estimates.csv) · [Annual comparison](annual_comparison/README.md) · [Nigeria states](nigeria_states/README.md)

[Age contributions](burden/deaths_by_age.png) · [DRC synthetic-cohort survival](burden/drc_survival.png) · [Aggregate outcome check](outcome_residuals_by_pfpr.png)

[Paper figure manifest/captions](paper_figures/CAPTIONS.md) · [Study flow](study_flow/CAPTION.md) · [Fit diagnostics](fit_diagnostics.csv)

These are primary results. Sahel, geographic, period and other sensitivity analyses are separate; they are not refitted by reporting. No TeX files are written. Reproduce with `Rscript R_cbh/primary/run_regional.R --report-only`.
