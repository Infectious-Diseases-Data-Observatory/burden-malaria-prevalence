# Current primary figures and tables

The revised 17-variable PfPR-ACM model uses 1,686,004 children, 5,465,305 child-band records and 75,726 deaths in 95 surveys and 34 countries.

Seven separate MAP gamma=2 age-band fits, including survey-region urban percentage, with fixed PfPR/time knots and the declared HIV, UNICEF and available-region substitutions. The raw-source data and HIV imputation were not refitted. All seven models passed the fitted-input and numerical checks.

![PfPR curves](pfpr_splines.png)

| Completed months | Records | Deaths (% of U5 deaths) | PfPR EDF | HR: 40% to 20% (95% interval) | HR: 20% to 0% (95% interval) |
|---|---:|---:|---:|---:|---:|
| <1 | 966,553 | 27,788 (36.7%) | 2.26 | 0.936 (0.906–0.966) | 0.938 (0.889–0.990) |
| 1-5 | 868,070 | 11,146 (14.7%) | 1.77 | 0.924 (0.886–0.963) | 0.908 (0.851–0.970) |
| 6-11 | 831,667 | 10,286 (13.6%) | 2.33 | 0.814 (0.773–0.857) | 0.755 (0.689–0.826) |
| 12-23 | 717,905 | 10,215 (13.5%) | 3.31 | 0.835 (0.780–0.895) | 0.578 (0.507–0.659) |
| 24-35 | 713,239 | 8,383 (11.1%) | 3.38 | 0.824 (0.765–0.887) | 0.529 (0.458–0.610) |
| 36-47 | 692,778 | 4,992 (6.6%) | 3.21 | 0.839 (0.773–0.912) | 0.555 (0.475–0.648) |
| 48-59 | 675,093 | 2,916 (3.8%) | 3.16 | 0.929 (0.838–1.030) | 0.564 (0.470–0.676) |
| **Total** | **5,465,305** | **75,726 (100.0%)** | — | — | — |

Deaths are observed deaths in the primary analysis sample; percentages use all 75,726 under-five deaths as the denominator. Displayed percentages use largest-remainder rounding to one decimal place and sum to 100.0%. Records are child–age-band observations, so a child can contribute multiple records. EDF: effective degrees of freedom of the PfPR spline. Both contrast columns are adjusted mortality hazard ratios for reducing PfPR from the first value to the second. Intervals are conditional on fitted smoothing parameters, exposure, the fixed HIV imputation and filled regional covariates; zero PfPR is below observed exposure support in every band.

[Table CSV](tables/age_band_results.csv) · [LaTeX source as plain text](tables/age_band_results.latex.txt)

[40% to 20% contrast figure](pfpr_40_to_20.png) · [Full contrasts and comparison with the previous adjustment](REPORT.md)

![Country estimates](burden/country_vs_ihme.png)

Country/IHME axes use identical base-10 scales starting at 1,000 deaths. All country estimates, signed age contributions and missing/support flags remain in the CSVs. Age intervals are conditional on smoothing parameters, exposure and filled covariates; marginal intervals are not summed into total intervals.

[Country totals](burden/country_totals.csv) · [Country-age estimates](burden/country_age_estimates.csv) · [Annual comparison](annual_comparison/README.md) · [Nigeria states](nigeria_states/README.md)

[Age contributions](burden/deaths_by_age.png) · [DRC synthetic-cohort survival](burden/drc_survival.png) · [Aggregate outcome check](outcome_residuals_by_pfpr.png)

[Paper figure manifest/captions](paper_figures/CAPTIONS.md) · [Study flow](study_flow/CAPTION.md) · [Fit diagnostics](fit_diagnostics.csv)

These are primary results. Sahel, geographic, period and other sensitivity analyses are separate; they are not refitted by reporting. No TeX files are written. Reproduce with `Rscript R_cbh/primary/run_regional.R --report-only`.
